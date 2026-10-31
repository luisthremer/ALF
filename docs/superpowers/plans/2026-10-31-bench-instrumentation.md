# Benchmark instrumentation for the cluster campaigns — Implementation Plan (revision 2)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Everything the Phase 8 campaigns need that can be built and tested on the Mac without cluster access: a per-sweep
time split, a BLAS ceiling microbenchmark, and a Slurm job generator with complete provenance. Revision 2 folds in the
plan review (15 findings; ledger).

**Architecture:** Exclusive (stack) timers compiled in only with `-DALF_TIMERS`, built as a separate binary like the trace
build (`builds/ALF_timers.out`); default output byte-identical. A standalone BLAS microbenchmark. A campaign generator
(python, run locally) that writes bash-only Slurm job scripts; provenance is captured on the compute node; logs are rsynced
back and collected/reported locally.

**Tech Stack:** Fortran 2008, bash (3.2-compatible on the Mac), python3 (local only), Slurm (cluster only).

**Requirements (user prompt, "Benchmarks (Phase 8)"):** time split per sweep (update, flush, wraps, stabilisation,
measurements); ZGEMM/ZGERU microbenchmark per layout; provenance in every log (machine, node, CPU, clock setting, compiler and
flags, MKL version and threading layer, fakeintel on/off, MPI, ranks × threads, model, L, β, Δτ, N_delay fixed/auto,
commit, repetitions); `MKL_NUM_THREADS = OMP_NUM_THREADS`; pinning per site documentation; whole exclusive nodes; ≥ 3
repetitions with median and spread; numbers in reports only from logs written by the scripts.

## Global Constraints
- Default build (no `-DALF_TIMERS`): **output** byte-identical (regress 13/13 IDENTICAL; the binary may differ); ctest all.
  `Timers_mod.F90` compiles with and without the flag (real module, its routines no-ops without the flag) and is linked into
  test 36. The flagged binary is built in a separate copy (`builds/timers/`), never in `Prog/` (it would poison `Prog/*.o`).
- Flagged build: trajectory IDENTICAL or ROUNDOFF (both acceptable under -ffast-math; state which); overhead ≤ 2 % s/sweep on
  the Mac for Hubbard L = 6, 12, 20 and tV L = 8 (worst case is small N); timer call cost measured and logged (T2).
- No cluster access; per-machine values are placeholders `[[...]]` in per-machine profiles; the generator exits non-zero if any
  `[[` remains in a non-local profile.
- Usual rules (configure.sh in the same command, no `set -u`; never push; do not touch USER_TESTS; commit trailer).

## Review Focus
1. Nested categories (flush inside update/hook/read point; hooks inside wraps) — exclusive timers; stack depth 0 at sweep end;
   categories + other = sweep total, and the sweep totals sum to info's CPU Time within 1 %. Pinned T1.
2. MPI aggregation collective on `Group_Comm` before the rank-0 write; per-rank other; rank 0, max, mean. Pinned T1.
3. No provenance field may silently be "unknown" where the machine can provide it — local tests use PATH shims and fixtures and
   assert the fields are filled. Pinned T3.

---

### Task T1: Timing split (`-DALF_TIMERS`)
**Files:** create `Prog/Timers_mod.F90`; modify `Prog/Makefile` (SRCS), `Prog/Wrapgr_mod.F90`, `Prog/Lowrank_update_mod.F90`,
`Prog/main.F90`, `testsuite/Prog.tests/CMakeLists.txt` (link for test 36), `testsuite/lowrank/README.md` (how to build
`builds/ALF_timers.out`, e.g. a `timers_build.sh` modelled on the trace build in compare_trajectory.sh).
**Design:**
- Exclusive stack timers: `timer_push(id)` pauses the current top and starts id; `timer_pop(id)` (asserts id is the top) stops
  it and resumes the parent; int64 `system_clock` counts (record `count_rate` in the output; if the rate is < 1e6/s, use
  `omp_get_wtime`). All calls under `#ifdef ALF_TIMERS`, also the `use Timers_mod`.
- Ids: `T_UPDATE` (ham%S0 + proposal ranf + Upgrade2/Upgrade_lowrank, excluding nested flushes), `T_FLUSH` (inside
  `Lowrank_flush` after its "nothing pending" return, so every caller is covered), `T_VWRAP` (Op_Wrapup/Op_Wrapdo + hooks),
  `T_HOP` (hopping wraps incl. cross-slice hops), `T_STAB` (WRAPUR/WRAPUL/CGR + precision control), `T_MEAS` (Obser,
  ObserT, Tau_m/Tau_p as a whole (tau modules not split), Hop_mod_Symm, reconstruction), `T_GLOBAL` (global moves,
  Global_tau, tempering exchange), `T_LANGEVIN_HMC` (whole Langevin/HMC branches), `T_SWEEP` (body of the `DO NSW` loop);
  per-bin overhead (Init_obs/Pr_obs/output) reported separately as `bin_overhead`.
- Counters: proposals, accepts, flushes by reason (capacity, slice end, ineligible vertex, read point/cross-slice), mean
  pending rank at flush.
- Output: per bin (timers reset at each bin start, so bin 1 can be dropped) one line per id
  `TIMER <bin> <id> <rank0> <max> <mean>` and `COUNTER <bin> <name> <rank0> <max> <mean>` in a file `timers.txt` (group rank
  0 writes after a collective `MPI_Allreduce`/`MPI_Reduce` on `Group_Comm` that every rank calls), plus `TIMER_SWEEPS <bin>
  <n>` (sweeps actually done; CPU_MAX may shorten bins) and `TIMER_RATE <count_rate>`; `info` gets one line
  `Timing split: see timers.txt (ALF_TIMERS)`.
- [ ] Implement. Verify: default build regress 13/13 IDENTICAL, ctest all; flagged build (separate copy) trajectory start and
  tv × delayed/submatrix IDENTICAL or ROUNDOFF; stack depth 0 at every sweep end (assert); Σ categories + other = T_SWEEP
  exactly by construction and Σ_bins (T_SWEEP + bin_overhead) within 1 % of info's CPU Time; overhead ≤ 2 % (bench.sh with
  `BENCH_BIN` flagged vs unflagged, Hubbard L = 6, 12, 20, tV L = 8, 3 reps); 2-rank MPI run prints rank0/max/mean.
  Commit `feat: exclusive timing split (-DALF_TIMERS)`.

### Task T3 (before T2): Campaign generator and provenance
**Files:** `testsuite/lowrank/campaign/`: `make_campaign.py` (campaign description as a python dict file: machine profile,
partition, nodes, layouts ranks×threads (assert ranks×threads = cores per node unless the point says otherwise), models with
keys, L list, schemes, N_delay list incl. 0/−1, repetitions, NSweep/NBin per point, wall-time estimate, memory estimate);
`profiles/{local,supermuc,helma}.sh` (module loads, account, partition, pinning, launcher, fakeintel wrapper, EAR option;
placeholders `[[...]]`); `provenance.sh`; `job_template.sh`; `collect.py`; `report.py`; `README.md`; `tests/` (PATH shims
`lscpu`, `scontrol`, `srun`, `sbatch`, `numactl`; MKL_VERBOSE fixture lines for SuperMUC AVX-512, Helma with and without
fakeintel).
**Design (bash-only on the cluster):**
- Build provenance at build time (`build_info.txt`: `$ALF_FC --version`, ALF_FLAGS_PROG, ALF_LIB, `module list`, binary
  sha256, `ldd`), referenced by every job.
- Job provenance on the compute node: `SLURM_CLUSTER_NAME`, `SLURMD_NODENAME`, `scontrol show job` (Exclusive, node list),
  `lscpu`, `numactl -H`, `free`, governor and `scaling_cur/max_freq` at start and end, boost flag, EAR option; MPI
  `mpirun --version`, `I_MPI_ROOT`; pinning (`I_MPI_DEBUG=4` map, `OMP_DISPLAY_AFFINITY=TRUE`, `taskset -cp`, `I_MPI_PIN*`,
  `KMP_*`, `OMP_PLACES/PROC_BIND`); `OMP_NUM_THREADS = MKL_NUM_THREADS`, `MKL_DYNAMIC=FALSE`; MKL version, ISA and threading
  layer from the first `MKL_VERBOSE=1` line of a tiny `blas_bench` run in the same environment (never in production runs);
  fakeintel: `LD_PRELOAD` set only in a per-rank wrapper, library sha256, check in `/proc/self/maps`; the job fails if
  fakeintel = on and the MKL ISA line is not AVX-512; model, L, β, Δτ, N_delay requested and used (parse info tolerantly:
  fixed `N_delay = n`, `auto (static n0 = …)`, `auto+tune`, `Auto N_delay: … chosen n`), commit (from info), repetition,
  setup time (start → first sweep), max RSS (`/usr/bin/time -v` or `sacct`), status file `OK`/`FAILED rc` written by a trap.
- One job = one layout × all schemes of a point (interleaved, same node); repetitions are separate jobs; `--exclusive`,
  `--no-requeue`; `submit_all.sh` idempotent (skips runs with status OK), respects a max-jobs-in-queue setting; run
  directories are never deleted.
- `collect.py`: one row per run with all required keys (fail on missing or duplicate keys; record log path + sha256,
  campaign id, script hash); excludes non-OK runs from statistics but lists them. `report.py`: medians, ranges (n = 3) and
  n_valid; refuses to mix commits, binaries or machines; node throughput (ranks / s per sweep) for layout comparisons;
  timing-split tables (dropping bin 1).
- [ ] Implement; local test with the `local` profile and shims: a mini campaign (Hubbard L = 8, standard and delayed, layouts
  1×1 and 2×2, 2 reps, the flagged binary) end to end; assert every provenance field is filled from the shims/fixtures
  (none "unknown"); a FAILED run is excluded and listed; a second `submit_all.sh` skips finished runs; the generator refuses a
  profile with `[[`. Commit `bench: campaign generator with provenance`.

### Task T2: BLAS ceiling microbenchmark
**Files:** `testsuite/lowrank/blas_bench.F90`, `testsuite/lowrank/blas_bench.sh` (built with ALF's compiler/flags/BLAS).
- [ ] Kernels: ZGEMM N×N×k (k ∈ {16, 64, 256}), ZGEMM N×N×N, ZGERU N×N, ZGEMV N×N, ZGEQP3 N×N; N up to 8000; warm-up reps,
  then median of ≥ 5; rotating buffers for the BLAS-2 kernels; GFLOP/s and GB/s; timer start/stop cost (T1's routines) in ns.
  "Per layout": `blas_bench.sh` launches one process per rank of the layout concurrently (MPI barrier or a start-time file)
  with `MKL_NUM_THREADS = OMP_NUM_THREADS = threads`, each reporting; the provenance of T3. Local run: N ∈ {400, 1600, 3600},
  1 thread and 1×4 layout. Commit `bench: BLAS ceiling microbenchmark`.

### Campaign-side notes for the cluster request (no task; used by the controller)
- Memory: N = 8000 needs GR 1 GB per flavour plus UDV stacks (tens of GB per rank): ranks per node limited by RAM (96 GB
  SuperMUC thin nodes); the generator checks the estimate against node RAM.
- Time: setup and initial stabilisation dominate at large N; pilot runs calibrate NSweep per point before a campaign.
