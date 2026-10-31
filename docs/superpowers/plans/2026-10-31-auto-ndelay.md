# Automatic N_delay — Implementation Plan (revision 2)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `N_delay = 0` chooses N_delay by a static rule (default, deterministic), `N_delay = -1` adds warm-up timing,
`N_delay ≥ 1` stays a fixed override (project 2, Phase 6).

**Architecture:** Rule, timing state and decision live in `Prog/Lowrank_update_mod.F90`; `Wrapgr_mod.F90` brackets the
sequential part of each WRAPGRUP/WRAPGRDO with begin/end calls (no-ops unless tuning); `main.F90` calls
`Lowrank_autotune_finish` once per completed sweep on every rank and `Lowrank_report_final` once at the end. Nothing changes
for `standard`. Revision 2 folds in the plan review (12 findings; ledger).

**Tech Stack:** Fortran 2008, MPI (Open MPI/MPICH locally), CMake/CTest, bash harnesses.

**Spec:** `docs/superpowers/specs/2026-10-31-auto-ndelay-design.md` (revision 2, A1–A8)

## Global Constraints
- As in `docs/superpowers/plans/2026-10-31-rankN.md`: regress IDENTICAL for every set; Upgrade2 untouched; RNG draws
  unchanged; TOL_UNIT/TOL_ILL; build in one shell command after sourcing configure.sh (no `set -u`); ≤ 8 cores; never push;
  commit trailer `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`; do not read or touch USER_TESTS.
- **N_delay semantics:** ≥ 1 fixed (unchanged); 0 = auto, static rule (new default in `QMC_runtime_var_mod.F90`); −1 = auto
  with warm-up timing; < −1 rejected by `validate_sequential_scheme` and by `Lowrank_init`.
- **Static rule (A2), in this order:**
  1. `raw = 2**nint(log2(Ndim/20.d0))` (delayed) or `2**floor(log2(Ndim/12.d0))` (submatrix); for Ndim < 12/20 use raw = 1.
  2. `n = max(16, min(U, raw))`, U = 128 (delayed) or 256 (submatrix).
  3. `n = min(n, max_pend)`, `max_pend = Σ_{n_op = Nt_start..Nt_end, eligible} max_f Op_V(n_op,f)%N_non_zero` (the most pending
     rank one slice can produce; G is flushed at every slice end); submatrix also `n = min(n, Ndim)`.
  4. Round **down** to a power of two; `n = max(n, 1)` (covers max_pend = 0).
- **Timing (A3), only N_delay = −1:** candidates `{n0/2, n0, 2·n0}`, each passed through steps 3–4 of the rule and ≥ 1,
  duplicates dropped; < 2 distinct → no timing ("single candidate"). Sweep 1 untimed; sweeps 2 and 3: the i-th call of
  WRAPGRUP/WRAPGRDO (counted from the start of sweep 2) uses candidate `mod(i-1, ncand)+1`; per call add the wall time
  (`system_clock` with `integer(8)` counts and rate) of the sequential part (from the start of the vertex loop — WRAPGRUP:
  after its hopping wrap; WRAPGRDO: after its Global_tau moves — to after the final slice-end flush, **excluding**
  Global_tau and the hopping wraps; fix round 1 of P6-1) and the number of `Upgrade_lowrank` calls (counter incremented
  inside `Upgrade_lowrank`). Decision after sweep 3: under `#ifdef MPI`, `MPI_Allreduce` (sum) of times and counts over
  `Group_Comm`, then group rank 0 picks and `MPI_Bcast`s the chosen index over `Group_Comm`; without MPI decide locally.
  Pick = fastest seconds-per-proposal if it is < 0.95 × n0's, else n0; any candidate with zero proposals → n0.
- **Log lines (exact),** written by group rank 0 to the scheme's info file (`info`, or `Temp_<igroup>/info` in tempering
  builds — store the file name at init in a module variable):
  `Sequential scheme: <NAME>, N_delay = <n>` (fixed, unchanged) / `..., N_delay = auto (static n0 = <n0>, Ndim = <Ndim>)` /
  `..., N_delay = auto+tune (static n0 = <n0>, Ndim = <Ndim>)`; after tuning:
  `Auto N_delay: candidates <a> <b> <c>: <t_a> <t_b> <t_c> s/proposal, chosen <n>` (ES10.3 times; fewer candidates if
  deduplicated) or `Auto N_delay: tuning skipped (<reason>), N_delay = <n0>` with reason ∈ {single candidate, no proposals,
  cross-slice, tuning incomplete}; submatrix with any non-diagonal eligible vertex appends ` (K_ND = 1: non-diagonal blocks
  flush immediately)` to the auto line.

## Review Focus
1. Every rank of a `Group_Comm` must call `Lowrank_autotune_finish` the same number of times (MPI collective inside). Pinned:
   P6-3 MPI and tempering runs.
2. With the cross-slice key on (if `Lowrank_cross_slice` exists in the tree), no timing (reason "cross-slice"). Pinned: P6-1 Step 5.
3. Restart from confin: counters start at 0 again → tuning repeats; allowed and logged. Pinned: P6-3 edge runs.
4. `max_pend < 16` (L=4 sets) and `max_pend = 0` (no eligible vertex). Pinned: P6-1 unit table.
5. bench rows with N_delay ≤ 0 must not enter "best fixed"/"Default N_delay" computations. Pinned: P6-1 Step 4.

---

### Task P6-1: Static rule, timing, interface, logging

**Files:** `Prog/Lowrank_update_mod.F90`, `Prog/Wrapgr_mod.F90`, `Prog/main.F90`, `Prog/QMC_runtime_var_mod.F90`,
`Scripts_and_Parameters_files/Start/parameters` (N_delay = 0 with a comment on 0/−1/≥1), `testsuite/Prog.tests/36-lowrank-unit.F90`,
`testsuite/Prog.tests/CMakeLists.txt`, `testsuite/lowrank/bench.sh`, `testsuite/lowrank/bench_report.py`,
`testsuite/lowrank/compare_trajectory.sh`, `testsuite/lowrank/README.md`, `CLAUDE.md` (supersede the Phase-1 "N_delay = 0 → error" note).

**Interfaces (produced):**
- `integer function Lowrank_static_N_delay(scheme, Ndim, max_pend)` — public, pure; scheme ∈ {SCHEME_DELAYED, SCHEME_SUBMATRIX}.
- `Lowrank_init(scheme_name, N_delay, Nt_start, Nt_end)` — the last two optional (default 1..size(Op_V,1)); WRAPGRUP/WRAPGRDO's lazy
  call passes `get_Nt_sequential_start()`, `get_Nt_sequential_end()`. For N_delay = −1 the arrays are sized for the largest
  candidate (2·n0 after steps 3–4), for N_delay = 0 for n0 (fix round 1 of P6-1), so `Lowrank_set_N_delay` never reallocates. `Lowrank_finalize` resets the tuning state.
- `subroutine Lowrank_set_N_delay(n)` — requires n_pend = 0 (else Terminate_on_error), 1 ≤ n ≤ allocated capacity.
- `subroutine Lowrank_autotune_slice_begin()` / `subroutine Lowrank_autotune_slice_end()` — called in WRAPGRUP/WRAPGRDO around
  the sequential part (see Timing); no-ops unless tuning is active.
- `subroutine Lowrank_autotune_finish(nsweeps_done)` — `integer(8)` count of completed sweeps (1-based, across bins). main.F90
  calls it on every rank immediately before the `ENDDO` that closes the `DO NSW` sweep loop, with a main-program counter
  incremented there (saturating at 4 is enough). Before Lowrank_init (Langevin, HMC-only, Sequential=.F.) it returns at once.
  Action: nsweeps_done = 1 → timing starts; = 3 → decision (collectives), `Lowrank_set_N_delay`, log line; otherwise nothing.
- `subroutine Lowrank_report_final()` — main.F90 calls it once after the bin loop (before the final output); if tuning was
  requested but not decided: log "tuning incomplete"; for P6-2: writes the resolved N_delay to data.h5 (HDF5 builds).
- `integer function Lowrank_current_N_delay()` — the value in use (diagnostic, tests, bench).

- [ ] **Step 1 (RED, not committed):** test 36 section `auto`: static-rule table for both schemes, Ndim ∈ {16, 64, 144, 400,
  1600, 3600, 8000} × max_pend ∈ {0, 4, 12, 1000} (expected values computed by hand from the rule, written as literals with
  a comment); `Lowrank_set_N_delay` down (n0/2) and up (2·n0) after an auto init, then a short delayed/submatrix run of the
  existing driver with the changed value (errors ≤ TOL_UNIT). Pending-state check: add argument parsing to the test program
  (`get_command_argument(1)`); argument `pending` initialises, makes one accept pending and calls `Lowrank_set_N_delay` (must
  terminate); register `add_test(36-lowrank-unit-pending 36-lowrank-unit pending)` with `WILL_FAIL TRUE`. Build: compile
  error or failures expected.
- [ ] **Step 2:** Implement (Global Constraints + Interfaces). GREEN: test 36 and the pending test; full ctest; regress IDENTICAL.
- [ ] **Step 3:** `compare_trajectory.sh`: add per-run keys `A:KEY=VALUE` / `B:KEY=VALUE` (applied only to run A/B; plain
  KEY=VALUE as now for both); then `compare_trajectory.sh --rebuild start delayed delayed A:N_delay=16 B:N_delay=-1` and
  `start standard delayed N_delay=0`, `... N_delay=-1`, the same for submatrix and for tv: IDENTICAL or ROUNDOFF.
- [ ] **Step 4:** bench: `bench.sh` accepts N_delay 0/−1 (run dir `nd0`/`ndm1`); when N_delay ≤ 0 the output line gets one more
  final field `nd_used=<n>` (parsed from the info log lines); `bench_report.py` accepts it (field counts 10, 13, 15, 16) and
  excludes N_delay ≤ 0 rows from "best fixed" and "Default N_delay"; existing reports regenerate byte-identically.
- [ ] **Step 5:** Audit: grep `tiers.sh`, `stat_runs.sh`, `compare_trajectory.sh`, `sets.sh`, `bench.sh`, `parallel_tier.sh` for
  runs that relied on the old default 32 (the Start template now says 0): keep their results meaningful (pass N_delay
  explicitly where pairing or comparability needs a fixed value; document). If `Lowrank_cross_slice` exists: with
  `N_delay=-1` and the key on, the log says "tuning skipped (cross-slice)".
- [ ] **Step 6:** Commit `feat: automatic N_delay (static rule by default, opt-in warm-up timing)`.

### Task P6-2: HDF5 provenance (code only; compiled on the clusters)
- [ ] In `main.F90`'s existing `#if defined(HDF5)` parameter block (written by `Irank_g == 0`, before the first sweep), add the
  scheme name and the **requested** N_delay as attributes of the root group of `data.h5`; in `Lowrank_report_final`, under
  `#ifdef HDF5` and group rank 0, the **resolved** N_delay. Use `h5fopen_f(file, H5F_ACC_RDWR_F, file_id, ierr)`,
  `h5ltset_attribute_string_f(file_id, "/", "Sequential_scheme", trim(name), ierr)`,
  `h5ltset_attribute_int_f(file_id, "/", "N_delay_requested"/"N_delay_used", buf, int(1,size_t), ierr)` with
  `integer :: buf(1)`, `h5fclose_f`. Follow how ALF already opens data.h5 there. The noMPI build without HDF5 must be
  unchanged (regress IDENTICAL). Record in `CLAUDE.md` under "Open questions"/cluster items: "HDF5 provenance: compile and
  check `h5dump -A data.h5` on the cluster". Commit `feat: scheme and N_delay in data.h5 (HDF5 builds)`.

### Task P6-3: Validation (controller-run)
- [ ] Spec Sec. 4 items 2–5: trajectories (P6-1 Step 3 plus kondo, jz, pv_adiab); edge cases (NSweep=1 NBin=1 with −1;
  Langevin; Nt_sequential_end < start; all vertices ineligible; a projector set; CPU_MAX truncation during tuning); local MPI
  2 ranks and tempering 2 groups with −1 (one decision per group, identical on all ranks); Mac timing (Mz Hubbard L = 12, 16,
  20; tV L = 8, 12; fixed {8, 16, 32, 64, 128} vs 0 and −1; 3 repetitions; pass if within max(5 %, 2 × spread) of the best
  fixed; how often −1's choice differs between repetitions). Results and commands in `CLAUDE.md` ("Phase 6 results"); commit.
