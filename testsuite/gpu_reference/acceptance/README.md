# Acceptance harness (ALF GPU port)

Implements [docs/gpu/VERIFICATION.md](../../../docs/gpu/VERIFICATION.md). Proposed 2026-10-29; in force after the
owner approves branch `gpu-port-process-lanes`.

## Files

| File | Where it runs | What it does |
|---|---|---|
| `acc_ship.sh HOST TICKET SHA [BASE]` | laptop | `git archive` of a commit (and its base) to `~/alf-gpu/work/t<TICKET>/cand/<sha12>/{src,base_src}`, read-only, tarball sha256 checked on both ends, `.alf_candidate` |
| `acc_build.sh TICKET SHA12 NAME TOOL ARGS` | machine | builds one variant into `cand/<sha12>/bin/NAME.out` + `NAME.build.json` (provenance); `base_*` names build the base |
| `acc_run.sh RUNDIR BIN CONFIG [--allow-debug] [VAR=v]` | machine | one run under `env -i` + whitelist, timeout 900 s, `provenance.json`; refuses debug/fault variables without `--allow-debug`; never overwrites |
| `acc_driver.py NNN.json --ticket --sha12 --who` | machine | prepares configs (fresh acceptance seeds), builds, runs (GPU one at a time, CPU runs in a pool of ≤ `cpu_parallel`), then `verdict.py` |
| `acc_gpunone.sh TICKET SHA12 OUT` | machine | GPU=none static checks: preprocessed sources of every changed Fortran file (computed list, 6 macro sets, errors gated) and ctest with its exit code |
| `verdict.py NNN.json RUNROOT` | machine | PASS/FAIL per gate and control from raw outputs, fail-closed; `verdict.json`, `verdict.md`; `--selfcheck` (CPU, seconds) must PASS whenever the script changes |
| `mutate.py IN OUT KIND` | machine | script-level positive controls (scale one small entry, shift, NaN, drop a group, rotate seeds, inflate/bias TAUDIAG) |
| `audit.py CLAIMS.json verdict.json` | anywhere | implementer's reported numbers vs recomputed ones, at printed digits |
| `refcache.py key/put/get/list/check` | machine | CPU reference cache in `~/alf-gpu/refcache/` (manager writes only) |
| `acc_check_commit.sh ACCEPTED [COMMIT]` | laptop | the commit equals the accepted candidate up to `docs/**` |
| `NNN.json` | in the candidate tree | pre-registered acceptance file (schema below); `207.draft.json` is a draft, not in force |

## Prerequisites on the machine
- `~/alf-gpu/env/{gnu,nvhpc,intel}.sh` (the tool's env file is sourced inside `env -i`).
- Prebuilt HDF5 per compiler in `~/alf-gpu/hdf5/` (copy of `HDF5/<compiler>` from a configure run), or set
  `ALF_HDF5_DIR`.
- The Python of `~/alf-gpu/venv` (pyALF, h5py, numpy, scipy) on `PATH` for the driver and verdicts.

## Typical T2 sequence (manager)

```sh
# laptop
testsuite/gpu_reference/acceptance/acc_ship.sh nividia 207 <cand-sha> <base-sha>
# machine (one ssh session; ≤ 15 min per run; unattended ~2 h for a 207-sized ticket)
. ~/alf-gpu/env/gnu.sh
C=~/alf-gpu/work/t207/cand/<sha12>
python3 $C/src/testsuite/gpu_reference/acceptance/acc_driver.py $C/src/testsuite/gpu_reference/acceptance/207.json \
        --ticket 207 --sha12 <sha12> --who mgr
bash $C/src/testsuite/gpu_reference/acceptance/acc_gpunone.sh 207 <sha12> $C/runs/gpunone
```

Then `audit.py` against the implementer's report, the `acceptance-verifier`, and `acc_check_commit.sh` before the
commit. Raw outputs stay until the owner's gate (VERIFICATION.md §9).

## Acceptance file (schema `alf-acc/1`)
- `ticket`, `candidate_sha` (`"self"`: the commit that contains the file), `base_sha`, `created`, `approvals` (owner
  quotes, dated), `machine`, `gpu_exclusive`, `timeout_s`, `disk_min_free_GB`, `cpu_parallel`, `false_alarm_budget`.
- `builds`: name → `tool`, `configure`, optional `flags_ext`, `gpu` (false = CPU-only, may run alongside).
- `configs`: name → `base` (a name in `configs.py`/`configs_stat.py`) and `overrides`.
- `runs`: id → `build`, `config`, `seeds` (`dev` = the committed seeds of the config; `accN` = fresh, drawn by the
  driver), `env` (`VAR=value`), `allow_debug`, `refcache`.
- `gates` and `controls`: `{id, type, …}` with the types of `verdict.py` (`selftest`, `faults`, `same_seed`, `replay`,
  `paired`, `precision`, `compare2`, `log_regex`); controls carry `expect` (normally `FAIL`). Every number a gate
  checks against is in the file, with a `_reason`.
