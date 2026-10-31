# Benchmark instrumentation for the cluster campaigns — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Everything the Phase 8 campaigns need that can be built and tested on the Mac without cluster access: a per-sweep
time split, a BLAS ceiling microbenchmark, and a Slurm job generator with complete provenance.

**Architecture:** Timers are compiled in only with `-DALF_TIMERS` (default build unchanged, byte-identical); they accumulate
`system_clock` (int64) intervals in a small module and print a "Timing split" block to `info` at the end of the run. The
microbenchmark is a standalone Fortran program linked against the same BLAS. The job generator writes Slurm scripts from a
campaign description; it is tested locally by generating scripts and running their payload under bash (no Slurm).

**Tech Stack:** Fortran 2008, bash, python3; Slurm (cluster only).

**Spec (requirements):** the user's prompt, section "Benchmarks (Phase 8)": time split per sweep (update, flush, wraps,
stabilisation, measurements); ZGEMM/ZGERU microbenchmark per layout; provenance in every log (machine, node, CPU, clock
setting, compiler and flags, MKL version and threading layer, fakeintel on/off, MPI, ranks × threads, model, L, β, Δτ,
N_delay fixed/auto, commit, repetitions); `MKL_NUM_THREADS = OMP_NUM_THREADS`; pinning; whole exclusive nodes; numbers in
reports only from logs written by the scripts.

## Global Constraints
- Default build (no `-DALF_TIMERS`) byte-identical: regress 13/13 IDENTICAL; ctest all.
- With `-DALF_TIMERS`: same Metropolis decisions (trajectory IDENTICAL) and overhead < 2 % s/sweep on the Mac (L=12, 20).
- No cluster access in this plan; nothing machine-specific hard-coded except in clearly marked per-machine config blocks
  with placeholders (`[[SuperMUC ...]]`, `[[Helma ...]]`) that the controller fills after the user's answer.
- Usual rules (build in one command after sourcing configure.sh, no `set -u`; never push; do not touch USER_TESTS; trailer).

## Review Focus
1. Timer categories must not overlap or double count (e.g. flush inside the update): pinned in T1 by a check that the
   categories add up to within 2 % of the total sweep time.
2. MPI: per-rank timers; the info block reports rank-0 values and the max over ranks of the group (load imbalance). Pinned T1.
3. Generator: every provenance field present in every log even when a value is unknown ("unknown" written, not omitted). Pinned T3.

---

### Task T1: Timing split (`-DALF_TIMERS`)
**Files:** create `Prog/Timers_mod.F90` (module `Timers_mod`: `timer_start(id)`, `timer_stop(id)`, `timers_report(file)`,
ids `T_UPDATE` (propose/accept incl. Upgrade2/Upgrade_lowrank), `T_FLUSH` (Lowrank_flush), `T_VWRAP` (Op_Wrapup/Op_Wrapdo
+ hooks), `T_HOP` (hopping wraps in WRAPGRUP/WRAPGRDO, incl. cross-slice hops), `T_STAB` (WRAPUR/WRAPUL/CGR and precision
control), `T_MEAS` (Obser, ObserT, Tau_m/Tau_p, Hop_mod_Symm, reconstruction), `T_GLOBAL` (global moves, Global_tau,
tempering exchange), `T_SWEEP` (whole sweep)); modify `Prog/Wrapgr_mod.F90`, `Prog/Lowrank_update_mod.F90`, `Prog/main.F90`,
`Prog/Makefile` (new module), every call under `#ifdef ALF_TIMERS`.
- [ ] Implement; `info` gains (only with the flag) `Timing split [s] (rank 0 / max over group): update ..., flush ..., vertex
  wraps ..., hopping ..., stabilisation ..., measurements ..., global ..., other ..., total ...` where other = total − sum.
- [ ] Verify: default build regress IDENTICAL; flagged build trajectory IDENTICAL (start, tv; delayed, submatrix); categories
  sum within 2 % of T_SWEEP (|other| / total ≤ 0.02 for Hubbard L=12 and tV L=8); overhead ≤ 2 % (bench.sh, 3 reps, with
  and without flag); MPI 2 ranks prints rank 0 and max. Commit `feat: timing split (-DALF_TIMERS)`.

### Task T2: BLAS ceiling microbenchmark
**Files:** create `testsuite/lowrank/blas_bench.F90` and `testsuite/lowrank/blas_bench.sh` (build with the same compiler
and BLAS flags as ALF from configure.sh; `blas_bench.sh N k reps` runs ZGEMM N×N×k (k = flush block, e.g. 16..256) and
N×N×N, ZGERU N×N, ZGEMV N×N; prints GFLOP/s per kernel with the provenance line of T3).
- [ ] Implement; run on the Mac for N ∈ {400, 1600, 3600}, k ∈ {16, 64, 256}, 1 thread: plausible numbers (ZGEMM ≫ ZGERU);
  commit `bench: BLAS ceiling microbenchmark`.

### Task T3: Campaign generator and provenance
**Files:** create `testsuite/lowrank/campaign/` with `make_campaign.py` (reads a campaign description `*.toml`-like or
python dict file: machine, partition, nodes, layouts ranks×threads, models with keys, L list, schemes, N_delay list incl.
0/−1, repetitions, NSweep/NBin per point, time limit), writes one Slurm script per job and a `submit_all.sh`, and a
`provenance.sh` sourced by every job that records: hostname, `lscpu`/`/proc/cpuinfo` model name, clock policy (EAR status on
SuperMUC via `[[...]]` placeholder command; `cpupower`/`--cpu-freq` on Helma), compiler and flags (from `configure.sh` env
and `ALF.out`'s `info`), MKL version (`MKL_VERBOSE=1` one-shot line) and threading layer (`MKL_THREADING_LAYER`), fakeintel
on/off (`LD_PRELOAD` contains the fakeintel path), MPI (`mpirun --version`), ranks × threads, `OMP_NUM_THREADS`,
`MKL_NUM_THREADS`, pinning variables, model, L, β, Δτ, N_delay (requested and used, from `info`), commit, repetition index;
`collect.py` gathers all job logs of a campaign into one raw table (one line per run, all provenance fields + the timing
split + s/sweep), and `report.py` produces Markdown tables with medians and spreads (numbers only from those logs).
Per-machine config blocks with placeholders for module loads, accounts, partitions, pinning.
- [ ] Implement; local test: a "local" machine profile (bash instead of sbatch, `mpirun` local) runs a 2-point mini campaign
  (Hubbard L=8, standard and delayed, 1×1 and 2×2 layouts, 2 reps) end to end: generated scripts, logs with every field
  (unknown values written as `unknown`), collect + report produce a table. Commit `bench: campaign generator with provenance`.
