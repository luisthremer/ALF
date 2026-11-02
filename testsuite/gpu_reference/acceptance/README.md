# Acceptance harness (ALF GPU port)

Implements [docs/gpu/VERIFICATION.md](../../../docs/gpu/VERIFICATION.md).

## Files

| File | Where | What it does |
|---|---|---|
| `acc_ship.sh HOST TICKET WHO CAND BASE PREREG` | laptop | checks the pre-registration (PREREG strict ancestor of CAND; tools unchanged since BASE except the added `TICKET*.json`, and since PREREG at all; existing tests unmodified unless allowed), then `git archive` of CAND and BASE to `~/alf-gpu/work/t<TICKET>/cand/<WHO>/<sha12>/`, hashes checked on both ends, `changed_files.txt`, `.alf_candidate` |
| `acc_build.sh TICKET WHO SHA12 NAME TOOL ARGS` | machine | extracts the re-hashed tarball, builds under `env -i`, refuses an HDF5 download; `bin/NAME.out` + `NAME.build.json` (with `tar_sha256_built`) |
| `acc_run.sh RUNDIR BIN CONFIG [--allow-debug] [VAR=v]` | machine | one run under `env -i` + whitelist (`OMP_TARGET_OFFLOAD=MANDATORY` for GPU builds, `ONEAPI_DEVICE_SELECTOR=level_zero:0` for `GPU=intel` builds; `ALF_GPU_*`, `ALFSH_*`, `ALF_LA_TEST_*` only with `--allow-debug`), timeout, GPU and user-load check, `provenance.json` with output hashes |
| `acc_qbuild.sh TICKET WHO SHA12 def\|log "CONFIGURE" QFLAGS` | machine | (207-5) offline quad-reference program + GNU CPU-oracle objects from the hash-checked tarball (the 207-2 scripts `t207_d5_build.sh`/`t207_d5_qbuild.sh` under `env -i`, work tree `cand/.../qw`); `bin/qref_<path>.qbuild.json` (tarball/candidate/program/script hashes, compiler). `ACC_QBUILD_HOME` only for a machine without `~/alf-gpu/env/gnu.sh` (intel) |
| `acc_qproc.sh TICKET WHO SHA12 RUNID def\|log QS NUL HOSTIN KEEP` | machine | (207-5) after a `acc_run.sh --qref` run: bridge records into `runs/<id>/qproc/` (207-2 layout, generated from the harness provenance), move the dumps, `t207_d5_proc.sh` (strict check, manifest, thinning by the gcd keep rule; interpreter `ACC_PYTHON`, default `~/alf-gpu/venv/bin/python3`), `acc_qproc.json` with the machine record; the driver calls it right after each QREF run, before the next run, then writes the run's `qref_stamp` line into the attempts log |
| `qref_gate.py`, `qref_selfcheck.py` | machine | (207-5) gate types `qref` (fills the instance fields of the pre-registered spec from the harness records, checks env, stamps, machine, build records, bridge records, then runs the approved `t207_d5_gate.py` and maps its exit status) and `h5_identical`; their lint; `machine` / `stamp` CLI; `--selfcheck` cases (part A synthetic; part B against the real gate with `verdict.py --selfcheck --qref-fixture CANDROOT`) |
| `acc_ctest.sh TICKET WHO SHA12 OUT TOOL [--base] [--regex R] [--repeat N] [--env V=x] -- CONFIGURE` | machine | unit tests from the verified tarball under `env -i`, `ctest -V`, one test at a time, build paths and addresses normalised; GPU configurations with MANDATORY offload under the GPU lock (`GPU=intel`: `ONEAPI_DEVICE_SELECTOR=level_zero:0`, stamp line `selector`, effective env in `env_rep_N.txt`; `--env` refuses the offload/device variables); `ctest.exit` 0 passed, 1 a test failed its own check, 2 invalid (crash, timeout, not run, build failure) |
| `acc_gpunone.sh TICKET WHO SHA12 OUT [--plant]` | machine | GPU=none: preprocessed changed Fortran (9 macro sets, own headers per side), counts of added/deleted, build and other changed files, ctest; `--plant` is its positive control |
| `acc_driver.py --ticket N --sha12 S --who mgr\|agent [--part P] [--reason R]` | machine | lint, configs and seeds (fresh `acc*` for mgr, `devN` for agent), builds, runs (GPU one at a time, lock), ctests, mutations, GPU=none, verdict; attempts log |
| `verdict.py ACC RUNROOT --candidate SHA [--seeds-record F] [--attempts F]` | machine | PASS / FAIL / INPUT_ERROR per gate, controls, BLOCKED items; `--selfcheck [--qref-fixture CANDROOT]` (without the fixture of a finished QREF dry run on that machine it ends INCOMPLETE, exit 4); `--lint ACC`; QREF gates need `--attempts` (stamps) and write the filled spec and the 207-2 gate report next to `verdict.json` |
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
  `active_lines` (lines every GPU run must print), `forbidden_lines` (lines no GPU run may print; also per build and
  per run, all three levels apply together; ignored for CPU runs), `blocked` (owner-decided BLOCKED items),
  `allowed_test_changes`.
- `builds`: name → `tool`, `configure`, `flags_ext`, `gpu` (must match a `GPU=` token in `configure`; lint),
  optional `forbidden_lines`, optional `active_lines` (replaces the file-wide list,
  e.g. for a base build without the new GPU path); names starting with `base_` build the base.
- `configs`: name → `base` (from `configs.py`/`configs_stat.py`), `overrides`, `parameters_sha256` (required when frozen).
- `runs`: id → `build`, `config`, `seeds` (`dev`, `accN`), `env`, `allow_debug`, `refcache`, and for debug runs
  only `active_lines` (replaces the build's/file's list, e.g. without the production line "GPU debug env: none").
- `ctests`: `{id, tool, configure, regex, repeat, env, control, gpu, base}` (acc_ctest.sh; `gpu` must match the
  `GPU=` token; `env` may not hold offload or device variables); `gpunone`: true. A `log_regex` check on a ctest
  directory also reads its `ctest.exit` (2 = INVALID) when the check names no `exit_file`.
- `mutations` (`mutate.py` on run directories) and `log_mutations` (planted failures in harness logs).
- `gates`: `{id, type, …, _reason}`; types `selftest`, `faults`, `same_seed`, `replay`, `paired`, `precision`,
  `compare2`, `log_regex`, `files_identical`, `qref`, `h5_identical`.
- QREF (ticket 207-5, docs/gpu/VERIFICATION.md 4b): `qbuilds`: path (`def`|`log`) -> `configure` (CPU, no GPU= token), `flags_ext`,
  `qflags` (`-DQR_PREPIVOT` / `-DSTABLOG`); `qref_env_allowed`: the names the env file of the GPU build sets (no debug, fault or steering
  names; lint); `qref_env_pinned` `{NAME: value}` + `_qref_env_pinned_reason`: a name the env file sets that the allow-list refuses by name, allowed
  only with exactly that value (intel: `{"NEO_CACHE_PERSISTENT": "0"}`; never offload, device-steering, debug, fault or trace names);
  `qref_gate_timeout_s` (hang guard of the gate subprocess, default 7200); a run with a `qref` entry `{path, qs, nul, hostin, keep}` (= the spec's) is a QREF run (needs `allow_debug`, only `ALF_GPU_*`
  in its env, whose role in the spec's `env_roles` must be the role it is listed for); gate `qref`: `spec` (the 207-2 gate spec, instance
  fields `binary_sha256`, `prog_sha`, `seeds` may be `"@harness"`, `parameters_sha256` only in a DRAFT; `"@harness"` seeds outside a DRAFT
  need `acc*` seeds in every run of the gate and its controls), `spec_template_sha256`, `runs` (main),
  `ctl_runs` `{role: [ids]}`, `tools_sha256` (the 207-2 scripts); gate `h5_identical`: `pairs [[on, off]], min_datasets`. A control of a
  `qref` gate replaces only `runs`: `mutations` of kind `qref_swap:dev:var8` of exactly the gate's main runs; it counts as a detection only if the
  device fails a bound and the copies' device equals the clean report's var8 (otherwise INPUT_ERROR "uninformative control"). QREF runs are
  never `base_*` builds; their env file must be the one of the build record. Stamps (`qref_stamp` lines in the attempts log) detect edits made
  after the stamp, not fabrication before it: the attempts log is writable by the same user. Host member: the offline program always prints `nvh`;
  specs use `host_member` `nvh` on both machines ("in-run host member, label fixed by the program").
- `controls`: `{id, like: <gate id>, replace: {runs|pairs|cases|a|b|checks}, expect}`. The control keeps its gate's
  criterion (only the inputs change); every gate needs one with `expect: FAIL`; a control counts only while its gate
  passes on the clean runs.
