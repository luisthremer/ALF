# CG-ALF T4: production guards of the Intel GPU production code

Read-only review, 2026-11-26. Nothing built, nothing merged. Source: fork tags `stp/production-mpi-intel-unstable-guard` (U, newest;
is a descendant-state of `...-large-prompt` (L) and `...-team-combined` (T) as far as production features go), all cut from origin master
521ddfca. Paths are relative to the tag. Prod scripts: `docs/supermuc_large_prompt/prod/` ("prod/").
CM = `Code-modularity-CPU` / `-GPU` at b0effdbd (GPU) and the CPU tip of the fork.

**Where the value is.** The Fortran side of U has few guards, and the three that matter are already on the CM lines. The real
production value is in the **shell layer** `prod/` (restart after SIGKILL, snapshots, retry, ledger, parameter checks). None of it
is Fortran, none of it touches the numerics, so porting is cheap and risk-free for the results.

## Table

Columns: ported? y/n/partly. Value for HPC production: H/M/L. Cost: S (<1 day), M, L. Rec: take/leave. Line: CPU, GPU, both.

| # | Feature | What it does | Where in U | On CM? | Value | Cost | Rec | Line |
|---|---|---|---|---|---|---|---|---|
| 1 | `ALF_UNSTABLE_THRESHOLD` | env var lowers the abort threshold of the G precision check (10 -> e.g. 1e-4); bad value exits 2 | control_mod.F90:97-121 | **y** (same code, control_mod.F90:97-113) | H | done | - | both |
| 2 | Production job sets threshold 1e-4 | healthy runs max 1e-10..3e-8 up to L=44; a corrupted G stops the run before it enters `data.h5` | prod/job_preempt.sh:39 | n (a job-script default, not code) | H | S | take (in template, see #9) | both |
| 3 | NaN abort in the precision check | NaN in G (new or wrapped) aborts with message; device variant `Control_PrecisionG_vals` | control_mod.F90 (`_vals`), gpu_stab_mod.F90:1306-1364 `isnan_bits` | **y** (bit-pattern NaN test, H1 fix; `_vals` on GPU) | H | done | - | both |
| 4 | Env allow-list: production refuses `ALF_GPU_*` debug/fault variables | a stray `ALF_GPU_FAULT*`, `_SELFTEST`, `_NS=0` in a job script cannot silently change a production run | gpu_state_mod.F90:344-420 | **y** (gpu_state_mod.F90:274, 451) | H | done | - | GPU |
| 5 | `ERROR_UNSTABLE_MATRIX` exit code and message "Unstable matrix" | distinct abort reason the job script can grep | runtime_error_mod.F90:60,99 | **y** (same enumerator) | M | done | - | both |
| 6 | `alf_exit(code)` | standard replacement for `call exit` | runtime_error_mod.F90:138 | coming (manager says so) | L | - | leave: ours maps every code != 0,2 to `error stop 1` | both |
| 7 | **Retry on "Unstable matrix"** | job loop: if exit 3 and the err log has a new "Unstable matrix" line, run recovery and restart ALF from the last bin, at most `MAX_UNSTABLE_RETRIES` (10) per job; each restart logged | prod/job_preempt.sh:70-100 | n | H on HPC (sporadic GPU aborts cost whole jobs; they had 'Unstable matrix' on a faulty node and from a queue-sync race, see #14) | S | **take** | both |
| 8 | **Restart after SIGKILL: `rc_recover`** | RUNNING file present -> every `Temp_i/data.h5` and `confout_*` read with `h5ls`; all bin datasets must have the same bin count; damaged set restored from newest of snapshot slots A/B and segment-start backup; copy to `*.rtmp`, `mv` only after all copies work; `.restoring` marker survives a kill during the restore | prod/recover.sh:341 (`rc_recover`), :435 (`rc_prepare`) | n (upstream only refuses to start while RUNNING exists, main.F90:204-225) | H for preempt queues, M otherwise | M (bash, needs `h5ls`; tested by #11) | **take** | both (PP/MPI runs) |
| 9 | Slurm job template for a preempt/requeue queue | `--requeue --open-mode=append`, mail on REQUEUE, `--exclude=<bad node>`, tool check (`h5ls`, `mpiexec`) before touching files, `MAX_RESTARTS` cap (100, exit 4), stop guard (exit 0 when every set has `TARGET_BINS`), `target_bins` override file | prod/job_preempt.sh:1-60 | n | H | S (needs SuperMUC/Helma module names per machine) | **take** as `manager_tools`/docs template, one per cluster | both |
| 10 | **Rolling snapshots** | background loop, every 30 s polls; per set takes a snapshot every 1800 s only if the bin is complete (all confouts present, newer than `data.h5` minus slack, files quiet), copies to `.snap/tmp`, re-checks sources unchanged, writes COMPLETE, replaces the older of two slots. Bounds the loss of a corrupting kill to ~30 min per set | prod/recover.sh:181 (`rc_snap_try`), :242 (`rc_snapshot_loop`) | n | M (only matters when kills corrupt files; 4 copies of `data.h5` per set on disk) | M | take, **second** after #8 | both |
| 11 | Restart test with the real binary | L=8, kill -9 after half the bins, truncated confout, truncated data.h5 during a snapshot, kill during a copy; result compared with an uninterrupted run (`cmp_h5dump.py`) | prod/restart_test.sh, recover_unit_test.sh, cmp_h5dump.py | n | H (proves #8/#10; prerequisite for using them) | M | take together with #8 | both |
| 12 | **HDF5 file locking pitfall** | a bare `h5ls`/`h5dump` on a live `data.h5` takes a read lock and the next `H5Fopen(RDWR)` of ALF fails; all reads go through `HDF5_USE_FILE_LOCKING=FALSE` | prod/README.md "HDF5 file locking", recover.sh `rc_h5ls` | n (not in CM docs) | H (knowledge: anyone monitoring a run by hand can kill it) | S | **take as a doc note now** | both |
| 13 | Duplicate-bin rule | a kill between `data.h5` append and `confout` write gives at most one duplicated bin per kill and set; analysis removes the second of two identical consecutive bins (`cmp_h5dump.py`) | prod/README.md, cmp_h5dump.py | n | M (analysis correctness; same issue exists in any ALF restart) | S | take (doc + script) | both |
| 14 | Sync diagnosis build | `ALF_OM_SYNC_ALWAYS`: `omp taskwait` + queue wait after every oneMKL call; found while hunting sporadic garbage G on Intel | diag/alf_la_mod_diag_sync_host.patch | n (oneMKL path not on CM yet) | M for step 2 gpu-intel | S | leave in the tag; **hand to gpu-intel worker** as a known symptom | GPU (Intel) |
| 15 | oneMKL error slots: `check_slots` | every dispatched call stores `info`; failures (`info<0`, or `>0` where not tolerated) abort with the routine name; test hook `om_poison_last` | alf_la_mod.F90:274 (`om_fail`), :502 | partly (linalg_mod on CM has its own `info` handling, per-vendor backends are step 2) | M | M (per backend) | take as a **design requirement** for gpu-intel/nvidia/amd notes, not as code | GPU |
| 16 | Device init guards | abort if no OpenMP offload device, no interop queue, queue not on Level Zero; `OMP_TARGET_OFFLOAD=MANDATORY` in the job script (no silent host fallback) | alf_la_mod.F90:361-382; prod/job_preempt.sh:37 | n for the device checks; `MANDATORY` is job-script level | H on GPU (a silent CPU fallback looks like a slow healthy run) | S | take `MANDATORY` into the template (#9); device checks go with step 2 | GPU |
| 17 | GPU-hour ledger | per run from `segments.log` x nodes + manual test list (sacct unreliable) | prod/ledger.sh, ledger_tests.txt | n | M (budget limits are real: 1000 GPU-h) | S | take as a tool in `manager_tools` | both |
| 18 | `segments.log` | one line per segment start/end: job id, restart count, bins per set, recovery action, exit code | recover.sh `rc_log_segment`, `rc_segment_end` | n | H (needed by #17 and for every post-mortem) | S (part of #8) | take with #8 | both |
| 19 | Parameter / seed checkers | `check_params.py` (only L, Beta, N_Phi, U, NSweep, N_delay, CPU_MAX, NBin may differ from the reference; N_FL, h0, Mz, Nwrap fixed), `check_seeds.py` (all seeds positive and distinct across run directories) | prod/check_params.py, check_seeds.py | n (project-specific reference; skill `jakob_params` exists) | M (distinct seeds across campaigns is a silent-correctness guard) | S | take `check_seeds.py` generic; `check_params.py` as a pattern only | both |
| 20 | `CPU_MAX` handling | the tag uses ALF's own `CPU_MAX` (a bit below the 7 days, clean stop at a bin boundary); nothing added | main.F90 (unchanged logic) | y (apply_CPU_MAX_policy, make_truncation) | - | - | - | both |
| 21 | MPI robustness | `gpu_state_setup` called on every rank (was rank 0 only); PP runs have independent walkers, the only traffic is `MPI_REDUCE` of bin averages and the seed distribution, so one rank dying only fails its own set | main.F90 (setup call), FINDINGS.md section 8 | partly (GPU line setup call to check, see "Not verified") | M | S | note only | GPU |
| 22 | Large-N phase | LU phase from the raw product at every CGR: unsafe for N >~ 2950 (production used L <= 44) | gpu_stab_mod.F90 | **y** (fixed on CM GPU, step 1d/1f) | - | - | **leave**, do not port | - |
| 23 | Lowrank / delayed scheme, Timers_mod, Trace_updates_mod, `Op_V_is_symmetric` fast path | performance and instrumentation, not guards | Prog/Lowrank_update_mod.F90, Timers_mod.F90, ... | Trace_updates_mod y; Timers on `cm/dev/timers`; rest is T6/T8 territory | - | - | out of T4 | - |

## Proposed features, in priority order (each a later dev/canvas feature, nothing implemented)

1. **`prod-restart-kit`** (#8, #9, #11, #12, #13, #18): `recover.sh` + `restart_test.sh` + `cmp_h5dump.py` + Slurm template under
   `Prog/..`-free place (`tools/production/` or `manager_tools/`), plus a short doc section with the HDF5 locking pitfall and the duplicate-bin rule.
   Acceptance: `restart_test.sh` reports PASS on one CPU machine (melmac or Helma) and one GPU machine. Machine-specific module lines become variables.
   No Fortran change, so no bitwise check needed. Value is high on every preempt queue and costs no performance.
2. **`prod-retry-on-unstable`** (#7, #2, #16 `MANDATORY` part): the retry loop around `mpiexec`, with `ALF_UNSTABLE_THRESHOLD` default 1e-4 and
   `OMP_TARGET_OFFLOAD=MANDATORY` in the template. Depends on #1 (needs `recover.sh`). Caution: a retry hides a systematic instability,
   so log every retry in `segments.log` and stop at `MAX_UNSTABLE_RETRIES`; a threshold of 1e-4 must be validated per model (U, L) before it is a default.
3. **`prod-snapshots`** (#10): rolling snapshots as an opt-in (`SNAP_INTERVAL`), after #1 has run once in production. Disk cost 4 x `data.h5` per set: check the quota first.
4. **`prod-ledger-seeds`** (#17, #19): `ledger.sh` and a generic `check_seeds.py` in `manager_tools/` or `tools/`.
5. **Design requirements for step 2 vendor backends** (#15, #16): per-call `info` check with routine name, device-presence checks, `MANDATORY`
   offload, no silent host fallback. Send to gpu-nvidia, gpu-amd, gpu-intel as a checklist; each backend note says how it meets them.

## Not taken, with reason

- #6 `alf_exit`: ours is not standard-conforming for exit codes (everything except 0 and 2 becomes 1). The coming CM version decides.
- #22 raw-product LU phase: replaced on the CM lines (step 1d/1f), do not port anything from the tag's stabilisation.
- #14: a diagnosis aid for one backend; only the symptom matters for gpu-intel (sporadic garbage G, "Unstable matrix", reproducible only on some nodes; one faulty node was found by excluding it).

## Bugs and risks noticed (KNOWN_ISSUES-style, none new on the CM lines)

- **C (production tooling, not code):** `recover.sh` keeps a set whose confouts belong to different bins after a kill. The walkers are
  independent chains, so no bias, but a **tempering** run (replicas exchange) would not be independent: the kit must refuse or document
  `TEMPERING` before it is used there. Class C now, becomes A if used with tempering unchanged.
- **C:** `ledger.sh` counts a job cancelled with `scancel` up to the next start or "now" (authors note ~1 h too much).
- No class A or B found on the CM lines in this review.

## Not verified

Read only; no build, no run, no check of the claims in the tag's own docs (restart test PASS, decision numbers). Whether the CM GPU line's
`gpu_state_setup` is called on every MPI rank was not checked. The `team-combined` and `large-prompt` tags were compared only by diffstat
against U (by diffstat L has the prod/ scripts in an older state, T has no prod/ directory; U adds the unstable retry, `ledger.sh`, `recover_unit_test.sh` and the recover fixes).
