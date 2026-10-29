# Acceptance harness (ALF GPU port)

Implements [docs/gpu/VERIFICATION.md](../../../docs/gpu/VERIFICATION.md).

## Files

| File | Where | What it does |
|---|---|---|
| `acc_ship.sh HOST TICKET WHO CAND BASE PREREG` | laptop | checks the pre-registration (PREREG strict ancestor of CAND; tools unchanged since BASE except the added `TICKET*.json`, and since PREREG at all; existing tests unmodified unless allowed), then `git archive` of CAND and BASE to `~/alf-gpu/work/t<TICKET>/cand/<WHO>/<sha12>/`, hashes checked on both ends, `changed_files.txt`, `.alf_candidate` |
| `acc_build.sh TICKET WHO SHA12 NAME TOOL ARGS` | machine | extracts the re-hashed tarball, builds under `env -i`, refuses an HDF5 download; `bin/NAME.out` + `NAME.build.json` (with `tar_sha256_built`) |
| `acc_run.sh RUNDIR BIN CONFIG [--allow-debug] [VAR=v]` | machine | one run under `env -i` + whitelist (`OMP_TARGET_OFFLOAD=MANDATORY` for GPU builds), timeout, GPU and user-load check, `provenance.json` with output hashes |
| `acc_ctest.sh TICKET WHO SHA12 OUT TOOL [--base] [--regex R] [--repeat N] [--env V=x] -- CONFIGURE` | machine | unit tests from the verified tarball under `env -i`; `ctest.exit` 0/1, 2 = invalid |
| `acc_gpunone.sh TICKET WHO SHA12 OUT [--plant]` | machine | GPU=none: preprocessed changed Fortran (9 macro sets, own headers per side), counts of added/deleted, build and other changed files, ctest; `--plant` is its positive control |
| `acc_driver.py --ticket N --sha12 S --who mgr\|agent [--part P] [--reason R]` | machine | lint, configs and seeds (fresh `acc*` for mgr, `devN` for agent), builds, runs (GPU one at a time, lock), ctests, mutations, GPU=none, verdict; attempts log |
| `verdict.py ACC RUNROOT --candidate SHA [--seeds-record F] [--attempts F]` | machine | PASS / FAIL / INPUT_ERROR per gate, controls, BLOCKED items; `--selfcheck`; `--lint ACC` |
| `mutate.py IN OUT KIND` | machine | script-level positive controls on run directories |
| `audit.py CLAIMS.json verdict.json` | anywhere | an implementer's reported numbers vs the recomputed ones |
| `refcache.py key/put/get/list/check` | machine | CPU reference cache in `~/alf-gpu/refcache/` (manager writes only) |
| `acc_check_commit.sh ACCEPTED [COMMIT]` | laptop | the commit equals the accepted candidate up to `docs/**` (scripts under `docs/gpu/scripts/` count as code) |
| `NNN.json`, `NNN.<part>.json` | candidate tree | pre-registered acceptance files; `*.draft.json` are drafts, not in force |

## Prerequisites on the machine
- `~/alf-gpu/env/{gnu,nvhpc,intel}.sh`; prebuilt HDF5 per compiler in `~/alf-gpu/hdf5/<compiler_vers>/` (or
  `ALF_HDF5_DIR`); the Python of `~/alf-gpu/venv` (pyALF, h5py, numpy, scipy) for the driver and the verdicts.

## Typical T2 sequence (manager)

```sh
# laptop: CAND = frozen candidate, BASE = accepted commit the ticket starts from, PREREG = pre-registration commit
testsuite/gpu_reference/acceptance/acc_ship.sh nividia 207 mgr <CAND> <BASE> <PREREG>
# machine (unattended; every run <= 15 min)
. ~/alf-gpu/env/gnu.sh
C=~/alf-gpu/work/t207/cand/mgr/<sha12>
python3 $C/src/testsuite/gpu_reference/acceptance/acc_driver.py --ticket 207 --sha12 <sha12> --who mgr
```

Then `audit.py` against the implementer's report, the `acceptance-verifier`, and `acc_check_commit.sh` before the
commit. Raw outputs stay until the owner's gate.

## Acceptance file (schema `alf-acc/1`)
- `ticket`, `candidate_sha` (`"self"`), `base_sha`, `created`, `_status`, `approvals` (owner quotes, dated), `machine`
  (`nividia` | `intel`), `gpu_exclusive`, `timeout_s`, `disk_min_free_GB`, `cpu_parallel`, `false_alarm_budget`,
  `active_lines` (lines every GPU run must print), `blocked` (owner-decided BLOCKED items), `allowed_test_changes`.
- `builds`: name → `tool`, `configure`, `flags_ext`, `gpu`, optional `active_lines` (replaces the file-wide list,
  e.g. for a base build without the new GPU path); names starting with `base_` build the base.
- `configs`: name → `base` (from `configs.py`/`configs_stat.py`), `overrides`, `parameters_sha256` (required when frozen).
- `runs`: id → `build`, `config`, `seeds` (`dev`, `accN`), `env`, `allow_debug`, `refcache`.
- `ctests`: `{id, tool, configure, regex, repeat, env, control, gpu, base}` (acc_ctest.sh); `gpunone`: true.
- `mutations` (`mutate.py` on run directories) and `log_mutations` (planted failures in harness logs).
- `gates`: `{id, type, …, _reason}`; types `selftest`, `faults`, `same_seed`, `replay`, `paired`, `precision`,
  `compare2`, `log_regex`, `files_identical`.
- `controls`: `{id, like: <gate id>, replace: {runs|pairs|cases|a|b|checks}, expect}`. The control keeps its gate's
  criterion (only the inputs change); every gate needs one with `expect: FAIL`; a control counts only while its gate
  passes on the clean runs.
