# Automatic N_delay — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `N_delay = 0` chooses N_delay by a static rule (default, deterministic), `N_delay = -1` adds warm-up timing,
`N_delay ≥ 1` stays a fixed override (project 2, Phase 6).

**Architecture:** Rule, timing state and decision live in `Prog/Lowrank_update_mod.F90`; `Wrapgr_mod.F90` brackets each
WRAPGRUP/WRAPGRDO with begin/end calls (no-ops unless tuning); `main.F90` calls `Lowrank_autotune_finish` after every sweep
on every rank. Nothing changes for `standard`.

**Tech Stack:** Fortran 2008, MPI (Open MPI locally), CMake/CTest, bash harnesses.

**Spec:** `docs/superpowers/specs/2026-10-31-auto-ndelay-design.md` (revision 2, decisions A1–A8)

## Global Constraints
- Same as `docs/superpowers/plans/2026-10-31-rankN.md` (regress IDENTICAL for every set; Upgrade2 untouched; RNG draws
  unchanged; tolerances; build commands; ≤ 8 cores; never push; Co-Authored-By line).
- `N_delay` semantics: ≥ 1 fixed; 0 static automatic; −1 automatic with timing; < −1 rejected in `validate_sequential_scheme`.
- Static rule and clamps exactly as spec A2; timing exactly as A3 (sweeps 2–3, round-robin, int64 `system_clock`,
  Allreduce over `Group_Comm`, > 5 % rule, zero proposals → n0).

## Review Focus
1. MPI: a rank that skips the collective (e.g. Sequential = .F. on one group in PARALLEL_PARAMS, CPU_MAX) — every rank of a
   `Group_Comm` must call `Lowrank_autotune_finish` with the same global sweep count. Pinned: Task P6-3 MPI/tempering runs.
2. `N_delay = -1` with cross-slice on (if kept): no timing, static n0 (A7). Pinned: Task P6-1 unit/run check.
3. Restart from `confin`: the global sweep counter starts at 0 again → tuning repeats (allowed); logged. Pinned: P6-3 edge runs.
4. Tiny models (`max_pend < 16`): value = max_pend rounded down, timing skipped if < 2 distinct candidates. Pinned: P6-1 unit table.
5. `bench_report.py` with `auto` rows: label, no numeric sort crash. Pinned: P6-1 step.

---

### Task P6-1: Static rule, timing, interface

**Files:** Modify `Prog/Lowrank_update_mod.F90`, `Prog/Wrapgr_mod.F90`, `Prog/main.F90`, `Prog/QMC_runtime_var_mod.F90`
(`validate_sequential_scheme`, default `N_delay = 0`), `Scripts_and_Parameters_files/Start/parameters` (comment),
`testsuite/Prog.tests/36-lowrank-unit.F90`, `testsuite/lowrank/bench.sh`, `testsuite/lowrank/bench_report.py`,
`testsuite/lowrank/README.md`, `CLAUDE.md` (supersede the Phase-1 "N_delay = 0 → error" note).

**Interfaces (produced):**
- `integer function Lowrank_static_N_delay(scheme, Ndim, max_pend)` (public, pure, A2) — unit-testable.
- `subroutine Lowrank_set_N_delay(n)` (flushed state required; terminates otherwise).
- `subroutine Lowrank_autotune_slice_begin()`, `subroutine Lowrank_autotune_slice_end()` (WRAPGRUP/WRAPGRDO entry/exit,
  no-ops unless tuning), `subroutine Lowrank_autotune_finish(isweep_global)` (every rank, every sweep).
- `integer function Lowrank_current_N_delay()` (diagnostic; the value in use).
- info line format (exact): `Sequential scheme: <NAME>, N_delay = auto (static n0 = <n0>, Ndim = <Ndim>)` for 0;
  `..., N_delay = auto+tune (static n0 = <n0>, Ndim = <Ndim>)` for −1, followed after tuning by
  `Auto N_delay: candidates <a> <b> <c>: <t_a> <t_b> <t_c> s/proposal, chosen <n>` (or the reason it was skipped).

- [ ] **Step 1:** Unit tests first (test 36 section `auto`): static rule table (both schemes; Ndim 64, 144, 400, 1600, 3600,
  8000; max_pend 4, 12, 1000), `Lowrank_set_N_delay` down/up within the allocation (ok) and with pending state (separate
  WILL_FAIL ctest or a guarded subprocess check), RED.
- [ ] **Step 2:** Implement; GREEN; full ctest; regress IDENTICAL; `compare_trajectory.sh start standard delayed N_delay=0`
  and `... N_delay=-1` IDENTICAL; same for submatrix and for tv.
- [ ] **Step 3:** `bench.sh`/`bench_report.py`: N_delay ≤ 0 accepted, run dir `nd0`/`ndm1`, output adds the chosen value
  (from info) as field `nd_used` — only when N_delay ≤ 0 (existing lines unchanged); reports regenerate byte-identically.
- [ ] **Step 4:** Commit `feat: automatic N_delay (static rule by default, opt-in warm-up timing)`.

### Task P6-2: HDF5 provenance (code only; compiled on the clusters)
- [ ] Behind `#ifdef HDF5`, write the scheme name and the N_delay in use (and whether auto) as attributes of `data.h5`
  at the point where ALF writes its parameters (find it in `main.F90` / the ham's `write_parameters_hdf5`). Must compile
  in the noMPI build without HDF5 (preprocessed out). Record "HDF5 provenance: compile + check on cluster" as a Phase-7
  item in `CLAUDE.md`. Commit `feat: scheme and N_delay in data.h5 (HDF5 builds)`.

### Task P6-3: Validation (controller-run)
- [ ] Spec Sec. 4 items 2–5 (trajectory, edge cases, local MPI/tempering, Mac timing with `bench.sh`), results with commands
  in `CLAUDE.md` ("Phase 6 results"); commit.
