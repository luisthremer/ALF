# Delayed and Submatrix Local Updates — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add runtime-selectable `delayed` and `submatrix` variants of ALF's sequential local update. They must be exact up to roundoff and faster, and `standard` must stay bit-identical.

**Architecture:** A new module `Prog/Lowrank_update_mod.F90` holds the pending low-rank state, the ratio/accept kernel, the wrap hook and the BLAS-3 flush. `WRAPGRUP`/`WRAPGRDO` keep their proposal loop. With a non-standard scheme they route eligible vertices to the new kernel, call a hook after each diagonal vertex wrap, and flush before any ineligible vertex, before Global_tau moves and before returning. A validation harness is built first (Phase 2) and run standard-vs-standard; each scheme is then added under it.

**Tech Stack:** Fortran 2008 (gfortran 16.2, `-cpp -O3 -ffast-math -fopenmp`), BLAS/LAPACK (Accelerate), CMake/CTest, bash, python3 + numpy (run with `python3 -I`), ALF's `Analysis/ana.out`.

**Spec:** `docs/superpowers/specs/2026-10-30-delayed-submatrix-design.md` (approved 2026-10-30, all four open questions accepted as proposed). Math: `docs/algorithm-notes.md`.

## Global Constraints

- Priorities: Markov-chain correctness > default behaviour unchanged > speed > elegance.
- `Sequential_scheme = "standard"` (default) must give **bit-identical** observable and `confout_0` files to baseline commit `f81b5cbc`.
- Parameters in `&VAR_QMC`: `Sequential_scheme` (`"standard"|"delayed"|"submatrix"`, case-insensitive), `N_delay` (integer ≥ 1, default 32 until Phase 5 retunes it). Both are optional.
- Eligible vertex: for every flavour `nf`, `Op_V(n,nf)%N_non_zero <= 1` and `Op_V(n,nf)%diag`. Ineligible vertices → flush + standard `Upgrade2`. One message is printed once if any vertex is ineligible.
- Ratio, Metropolis weight, ranf draw and phase update use the same expressions, in the same order, as `Upgrade2` FINAL mode (notes §1). The RNG consumption order is identical to standard.
- No pending state may survive a return from `WRAPGRUP`/`WRAPGRDO` (`Lowrank_assert_flushed`, always on).
- Only these existing files change: `Prog/Wrapgr_mod.F90`, `Prog/QMC_runtime_var_mod.F90`, `Prog/Makefile` (SRCS), `Prog/upgrade_mod.F90` (`#ifdef ALF_TRACE_UPDATES` only), `Scripts_and_Parameters_files/Start/parameters`, `testsuite/Prog.tests/CMakeLists.txt`. **Not** `main.F90`, not the Hamiltonians.
- Unit tolerance `TOL_UNIT = 1.0d-10` (relative; justified in Task 3). Trajectory tolerance for observable bins: `1e-8` relative.
- Runs: minutes; longer jobs in the background; ask the user before anything expected to take more than 30 min. Never push. Don't install software. Use `/bin/rm` (`rm` is aliased to `trash`).
- Commit after each task whose checks pass, ending the message with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Build: `cd ALF && source configure.sh GNU noMPI NO-INTERACTIVE >/dev/null 2>&1 && make program -j8` (env vars do not persist between shell calls). Extra flags via `ALF_FLAGS_EXT="-D..."` exported before sourcing.

## Review Focus

1. **Down-sweep with submatrix**: the un-wrap with the new field touches the just-added site → `Γ̃⁻¹` must be rescaled. Expect exact agreement in both sweep directions (Task 11 unit test runs both directions).
2. **Same site twice within one block** (more than one vertex per site per slice): submatrix must flush first; delayed needs nothing. Expect agreement (Tasks 7 and 11 unit tests include a vertex list with repeated sites).
3. **N_delay extremes**: `N_delay=1` (flush after every accept) and `N_delay` ≥ the number of vertices per slice (flush only at slice end). Expect identical results (Tasks 7 and 11 unit tests loop over N_delay ∈ {1, 7, 32, 1000}).
4. **Models with ineligible vertices** (e.g. tV bond operators, rank 2): scheme `delayed`/`submatrix` must fall back for each such vertex, print one message, and, if no vertex is eligible, give output bit-identical to standard (Task 8 integration check).
5. **Continuous (type 3) fields and `Propose_S0`**: these change how many random numbers are drawn before the update. Expect trajectories identical to standard (Task 8/12 trajectory runs include a continuous-HS set).

---

## Phase 2 — Validation harness (no algorithm code)

### Task 1: Baseline build and bit-identity regression script

**Files:**
- Create: `testsuite/lowrank/regress_vs_baseline.sh`, `testsuite/lowrank/make_params.py`, `testsuite/lowrank/README.md`
- Modify: `CLAUDE.md` (Commands section)

**Interfaces:**
- Produces: `make_params.py OUT_DIR KEY=VALUE ...`: copies `Scripts_and_Parameters_files/Start/{parameters,seeds}` to OUT_DIR and sets each `KEY = VALUE` in whichever namelist contains KEY (case-insensitive; a key not present, such as `Sequential_scheme`, is appended to `&VAR_QMC`). `SEED=<int>` replaces the first seed in `seeds`.
- Produces: `regress_vs_baseline.sh [extra KEY=VALUE ...]`: ensures a worktree `../ALF_baseline` at `f81b5cbc` exists and is built, then for each named parameter set below runs baseline `ALF.out` and branch `ALF.out` in `runs/regress/<set>/{base,branch}` and `cmp`s every file except `info*`. Exit 0 iff all are identical; prints one line per set (`<set> IDENTICAL|DIFFERENT`). Each set is defined as common keys plus optional branch-only keys (keys the baseline cannot parse). The env var `BRANCH_EXTRA="K=V ..."` adds branch-only keys to every set.
- Parameter sets (short: `NSweep=5 NBin=2`): `start` (Start defaults), `doped` (`ham_chem=-0.5 L1=4 L2=4 Beta=4`), `density` (`Mz=.F.`), `continuous` (`Continuous=.T.`), `tv` (`ham_name="tV"` with `Model="tV"`, the tV defaults).

- [ ] **Step 1:** Write `make_params.py` and `regress_vs_baseline.sh` as specified (bash `set -euo pipefail`; the baseline is built in its own dir with its own `source configure.sh`).
- [ ] **Step 2:** Run `testsuite/lowrank/regress_vs_baseline.sh`. Expected: 5 lines `IDENTICAL`, exit 0 (the branch has no code changes yet).
- [ ] **Step 3:** Sanity-check the harness: run once with an extra `NSweep=6` applied to the branch side only (`BRANCH_EXTRA=NSweep=6`); expected: `DIFFERENT` for every set, exit ≠ 0.
- [ ] **Step 4:** Add the regression command to `CLAUDE.md`, then commit `test: add bit-identity regression harness vs baseline`.

### Task 2: Runtime parameters `Sequential_scheme`, `N_delay`

**Files:**
- Modify: `Prog/QMC_runtime_var_mod.F90` (variables, namelist, defaults, getters/setters, MPI bcast), `Scripts_and_Parameters_files/Start/parameters` (document both keys with defaults, in `&VAR_QMC`)

**Interfaces:**
- Produces: `character(len=64) function get_Sequential_scheme()`, `subroutine set_Sequential_scheme(val)`, `integer function get_N_delay()`, `subroutine set_N_delay(val)`. Defaults `"standard"`, `32`. MPI: `MPI_BCAST(Sequential_scheme,64,MPI_CHARACTER,...)`, `MPI_BCAST(N_delay,1,MPI_INTEGER,...)` inside `broadcast_QMC_runtime_var`.

- [ ] **Step 1:** Add a regression set `explicit` = `start` with branch-only keys `Sequential_scheme="standard" N_delay=32`. Run it: expected FAIL (the branch aborts on an unknown namelist key).
- [ ] **Step 2:** Implement the parameters following the existing `Amplitude`/`ham_name` patterns.
- [ ] **Step 3:** Rebuild; run `regress_vs_baseline.sh`. Expected: all 6 sets `IDENTICAL`.
- [ ] **Step 4:** Rebuild the CTest suite and run `ctest`. Expected: 39/39 pass. Commit `feat: add Sequential_scheme and N_delay runtime parameters`.

### Task 3: Unit-level harness (reference = real `Upgrade2`), standard vs standard

**Files:**
- Create: `testsuite/Prog.tests/36-lowrank-unit.F90`
- Modify: `testsuite/Prog.tests/CMakeLists.txt` (new target linking all `Prog/*.o` except `main.o`, plus `Prog/Hamiltonians/*.o`, via `file(GLOB)`; `add_test(36-lowrank-unit ...)`)
- Modify: `CLAUDE.md` (Tolerances section)

**Interfaces:**
- Consumes: `Upgrade2` (INTERMEDIATE mode returns the full ratio in `Prev_Ratiotot` and always accepts), `Op_Wrapup`, `Op_Wrapdo`, `Predefined_Int_U_MZ`, `Predefined_Int_U_SUN`, `Predefined_Int_U_MZ_continuous_HS`, the `Hamiltonian_main` globals (`Op_V, nsigma, Ndim, N_FL, N_FL_eff, Calc_Fl_map, N_SUN, reconstruction_needed=.false.`), `Ranset`.
- Produces (inside the test program, reused by Tasks 7 and 11): a driver `run_case(case_id, scheme, N_delay, direction, max_ratio_err, max_G_err)` over these cases:
  `mz` (N_FL=2, type 2, Ndim=64), `density` (N_FL=1, N_SUN=2, complex g from `Predefined_Int_U_SUN` with U>0, Ndim=64), `continuous` (type 3, Ndim=48), `repeat` (Mz vertices with each site listed twice per slice, Ndim=32).
  Each case has 3 slices. Directions `up` (per vertex: `Op_Wrapup` type 1 old φ → hook → propose → decide → `Op_Wrapup` type 2) and `down` (`Op_Wrapdo` type 2 → propose → decide → `Op_Wrapdo` type 1 new φ → hook). Between slices: flush, then the same random dense similarity applied to both G's. Proposals come from `nsigma%flip` drawn once and shared; decisions come from a fixed accept mask (about 60 % accept, seeded). Reference ratio: `Upgrade2(...,"Intermediate")` on G_ref, or on a copy for rejected proposals, with `nsigma%f(n,nt)` restored. Initial G: `(1+A)^{-1}` with `A = Q diag(e^{x_i}) Q^{-1}`, x_i uniform in [−4,4], random complex Q with condition number ≲ 10.
  Errors: `max_ratio_err` = max over proposals of |R_cand − R_ref|/|R_ref|; `max_G_err` = max over slice ends of max|G_cand − G_ref| / max|G_ref|.
- Scheme `standard` candidate = the same `Upgrade2` path on a second copy (proves the harness is deterministic: errors must be exactly 0).

- [ ] **Step 1:** Write the test with the `standard` candidate only. It prints a table (case, direction, N_delay, max_ratio_err, max_G_err) and `stop 1` if any error exceeds `TOL_UNIT = 1.0d-10`.
- [ ] **Step 2:** Make the harness fail on purpose to prove it can fail: build once with `-DLOWRANK_UNIT_SELFTEST`, which perturbs one candidate G entry by 1e-8. Expected: exit 1 with that case flagged. (Keep the `#ifdef`.)
- [ ] **Step 3:** Run the normal build. Expected: all errors exactly 0, exit 0, `ctest -R 36-lowrank-unit` passes.
- [ ] **Step 4:** Tolerance measurement: add a report of the **standard path's own drift** against an exact recomputation `(1+A')^{-1}` (A' tracked explicitly, as in `docs/scripts/check_formulas.py`) for every case. Record the measured drift in `CLAUDE.md` with the justification: both candidate and reference are within roughly that drift of exact, so `TOL_UNIT` must be ≥ 2× drift with margin. If the drift exceeds 1e-12, stop and report to the user before going on.
- [ ] **Step 5:** Commit `test: add unit-level low-rank update harness (standard vs standard)`.

### Task 4: Trace instrumentation and trajectory harness

**Files:**
- Create: `Prog/Trace_updates_mod.F90` (added to `SRCS` in `Prog/Makefile`), `testsuite/lowrank/compare_trajectory.sh`, `testsuite/lowrank/compare_runs.py`, `testsuite/lowrank/trace_diff.py`
- Modify: `Prog/upgrade_mod.F90` (only inside `#ifdef ALF_TRACE_UPDATES`)

**Interfaces:**
- Produces: `subroutine Trace_update_write(nt, n_op, Ratiotot, Weight, ranf_val, toggle)` in module `Trace_updates_mod`. It lazily opens `trace_updates.txt` (append) and writes one line `nt n_op Re(R) Im(R) Weight ranf toggle` (`ES24.16`). In `Upgrade2`, under `#ifdef ALF_TRACE_UPDATES` only: the ranf draw is stored in a local, then traced in FINAL mode. The `#else` branch is the original line, unchanged.
- Produces: `compare_trajectory.sh SET SCHEME_A SCHEME_B [KEY=VALUE ...]`: builds (or reuses) a trace binary `builds/ALF_trace.out` (`ALF_FLAGS_EXT=-DALF_TRACE_UPDATES`; the normal binary is rebuilt afterwards) and runs both schemes with identical seeds and `NSweep=3 NBin=2` in `runs/traj/<set>/<scheme>`. It then calls `compare_runs.py` (field-by-field `confout_0` equality; every observable bin file compared with relative tolerance `1e-8`; for a non-standard scheme, the run's `info` must contain the line `Sequential scheme: <NAME>`, upper case, so a run that silently fell back to standard fails) and `trace_diff.py` (first differing toggle; for it prints `|Weight−ranf|/Weight` for both runs and classifies it as `ROUNDOFF` if ≤ 1e-10, else `BUG`). Exit 0 iff confout is identical, the bins are within tolerance and there are no `BUG` lines.

- [ ] **Step 1:** Run the regression `regress_vs_baseline.sh` after adding the `#ifdef` code and the new module. Expected: all `IDENTICAL` (the trace is off by default).
- [ ] **Step 2:** `compare_trajectory.sh start standard standard`. Expected: identical traces (`diff` clean), exit 0.
- [ ] **Step 3:** Prove the analyser fires: in a scratch copy of one trace, flip one toggle with a large margin. Expected: `trace_diff.py` reports `BUG` and exits ≠ 0. Do the same with a margin of 1e-13; expected: `ROUNDOFF`.
- [ ] **Step 4:** Commit `test: add update trace and trajectory comparison harness`.

### Task 5: Statistical and stability harness

**Files:**
- Create: `testsuite/lowrank/stat_runs.sh`, `testsuite/lowrank/stat_compare.py`

**Interfaces:**
- Produces: `stat_runs.sh SET SCHEME SEED` runs one statistical job in `runs/stat/<set>/<scheme>_s<SEED>` (to be started with `run_in_background`), then `ana.out` on `Ener_scal Pot_scal Kin_scal Part_scal SpinZ_eq Green_eq` (plus `*_scal` files as produced).
  Sets (all Mz unless noted, `NSweep=50 NBin=40 Ltau=0`):
  `half` = 6×6, U=4, β=5, μ=0; `doped` = 6×6, U=4, β=4, μ chosen in Step 1; `density` = `Mz=.F.`, 6×6, U=4, β=4, μ=0.
- Produces: `stat_compare.py DIR_A DIR_B` reads the analysed `*J` / `*JK` / `*JR` outputs and computes `z = (a−b)/sqrt(σa²+σb²)` for every scalar (energy, potential energy → double occupancy, kinetic energy, average sign/phase) and every k-point entry of `SpinZ_eq` and `Green_eq`. Pass: all scalar |z| < 3; for arrays, the fraction with |z| > 2 is ≤ 10 % and max |z| < 4.5. It also extracts `Precision Green Mean/Max` and `Precision Phase` from both `info` files: report them, and flag (without failing) if B's max is more than 10× A's.
  The two runs being compared **use different seeds** (independent chains).

- [ ] **Step 1:** Pilot for `doped`: short runs (`NSweep=20 NBin=5`) at μ ∈ {−0.3, −0.5, −0.8}; pick the μ with the average sign in [0.3, 0.8]. Record it in `CLAUDE.md`.
- [ ] **Step 2:** Run `standard` with seeds 1 and 2 for all three sets in the background (expected a few minutes each). `stat_compare.py`: expected PASS for all three sets. This validates the z-test on two independent standard chains.
- [ ] **Step 3:** Prove the comparator fires: compare `half` standard against `half` with U=4.5. Expected: FAIL on energy.
- [ ] **Step 4:** Commit `test: add statistical comparison harness`. Phase 2 is done: run superpowers:verification-before-completion and show the outputs of Tasks 1–5.

## Phase 3 — Delayed update

### Task 6: `Lowrank_update_mod`: state, eligibility, init

**Files:**
- Create: `Prog/Lowrank_update_mod.F90` (add to `SRCS` after `upgrade_mod.F90`)
- Test: `testsuite/Prog.tests/36-lowrank-unit.F90` (eligibility section)

**Interfaces:**
- Produces (public): `Integer, parameter :: SCHEME_STANDARD=0, SCHEME_DELAYED=1, SCHEME_SUBMATRIX=2`;
  `subroutine Lowrank_init(scheme_name, N_delay)` (`character(len=*)`, `integer`; it upper-cases, validates (error-terminates on an unknown name or `N_delay < 1`), computes `eligible(:)` from `Op_V`, allocates per-flavour state; **only for a non-standard scheme** it prints `"Sequential scheme: <NAME>, N_delay = <n>"` to `output_unit` and appends it to `info` (for standard nothing is printed and nothing allocated), and if any vertex is ineligible prints once `"Sequential scheme <NAME>: <k> of <M> vertices are not diagonal rank-1; they use the standard update (first: n=<i>)"`);
  `subroutine Lowrank_finalize()`; `logical function Lowrank_active()`; `logical function Lowrank_eligible(n_op)`; `integer function Lowrank_pending()`; `subroutine Lowrank_assert_flushed(where)` (error-terminates with `where` if pending > 0).
- State per computed flavour: delayed `X(Ndim, N_delay)`, `Y(Ndim, N_delay)`; submatrix `S(N_delay)`, `slot(Ndim)`, `Dvec(N_delay)`, `Ginv(N_delay,N_delay)`; `n_pend`. Column counts are in vectors (rank-N_P-ready).

- [ ] **Step 1:** Add an eligibility test to `36-lowrank-unit`: Mz vertices → eligible; a `Predefined_Int_V_SUN` bond vertex (rank 2) → ineligible; a non-diagonal rank-1 operator built by hand (N=2, O = (1 1;1 1)/2) → ineligible. Run: expected compile failure (module missing).
- [ ] **Step 2:** Implement the module skeleton and these routines. Build and run the test: expected PASS. Run `regress_vs_baseline.sh`: all `IDENTICAL`.
- [ ] **Step 3:** Commit `feat: add Lowrank_update_mod skeleton with eligibility`.

### Task 7: Delayed kernel, wrap hook, flush

**Files:**
- Modify: `Prog/Lowrank_update_mod.F90`; `testsuite/Prog.tests/36-lowrank-unit.F90` (add the `delayed` candidate)

**Interfaces:**
- Produces:
  `subroutine Lowrank_propose(GR, n_op, nt, Hs_new, Ratiotot)`: ratio only, no RNG, no state change except caching per-flavour `Δ_f`, `d_f`. Ratio scalars exactly as `Upgrade2` (incl. `weight_reconstruction` when `reconstruction_needed`).
  `subroutine Lowrank_accept(GR, n_op, nt, Hs_new)`: appends `(Δ/d)c` to X and `e_p − r` to Y per flavour (ZGEMV), sets `nsigma%f(n_op,nt) = Hs_new`, flushes if `n_pend == N_delay`.
  `subroutine Upgrade_lowrank(GR, n_op, nt, Phase, Hs_new, S0_ratio, T0_proposal_ratio, toggle)`: propose → Weight/ranf/Phase exactly as Upgrade2 FINAL (+ trace call under `#ifdef ALF_TRACE_UPDATES`) → accept if toggled → `Control_upgrade(toggle)`, `Control_upgrade_eff(toggle)`.
  `subroutine Lowrank_wrap_hook(n_op, nt, HS_field, dir)`: `dir=+1` after `Op_Wrapup` type 1, `dir=−1` after `Op_Wrapdo` type 1. For each flavour and each `i ≤ Op%N`: `t = exp(dir·g_loc·φ(HS_field)·E(i))` (`g_loc = g_t(nt)` if allocated), `X(P(i),1:j) *= t`, `Y(P(i),1:j) /= t`.
  `subroutine Lowrank_flush(GR)`: `ZGEMM('N','T',Ndim,Ndim,j,−1,X,Ndim,Y,Ndim,1,GR(:,:,nf),Ndim)`; reset counters.
  `subroutine Lowrank_current_G(GR, nf, Gcur)`: materialise `G0 − XYᵀ` without flushing (for tests/debug).
  GR is `Complex(Kind=Kind(0.d0)), intent(inout) :: GR(Ndim,Ndim,N_FL)`.

- [ ] **Step 1:** Add the `delayed` candidate to `run_case` for all cases × {up, down} × N_delay ∈ {1, 7, 32, 1000}. It uses `Lowrank_propose`/`Lowrank_accept` with the shared forced decisions, the hook after each type-1 wrap, and `Lowrank_flush` at slice end. Run: expected FAIL (routines missing).
- [ ] **Step 2:** Implement the routines above.
- [ ] **Step 3:** Run the unit test. Expected: every `max_ratio_err`, `max_G_err` ≤ `TOL_UNIT`; print the values.
- [ ] **Step 4:** Add a kernel-equivalence check to the test: from a flushed state with the same RNG seed, `Upgrade2(FINAL)` and `Upgrade_lowrank` on 200 proposals must give identical `toggle`, `|ΔPhase| ≤ 1e-12`, and the same next `ranf_wrap()` value afterwards (same RNG consumption). Run: expected PASS.
- [ ] **Step 5:** Commit `feat: delayed update kernel, wrap hook and BLAS-3 flush`.

### Task 8: Wire the schemes into `WRAPGRUP`/`WRAPGRDO`

**Files:**
- Modify: `Prog/Wrapgr_mod.F90`

**Interfaces:**
- Consumes: Task 6/7 public routines; `get_Sequential_scheme()`, `get_N_delay()` (Wrapgr_mod gains `use QMC_runtime_var, only: ...` and `use Lowrank_update_mod`).
- Changes (exactly; the standard path executes no new arithmetic):
  `Wrapgr_alloc`: `call Lowrank_init(get_Sequential_scheme(), get_N_delay())`. `Wrapgr_dealloc`: `call Lowrank_finalize()`.
  `WRAPGRUP` loop, per `n`: if `Lowrank_active() .and. .not. Lowrank_eligible(n)`: `Lowrank_flush(GR)` before the type-1 wrap. After the type-1 wraps of an eligible n: `Lowrank_wrap_hook(n, ntau1, HS_Field, +1)`. The `Upgrade2` call becomes `if (Lowrank_active() .and. Lowrank_eligible(n)) then Upgrade_lowrank(...) else Upgrade2(...)`. After the loop and before the `N_Global_tau` block: `if (Lowrank_active()) Lowrank_flush(GR)`. At the end: `Lowrank_assert_flushed("WRAPGRUP")`.
  `WRAPGRDO`: the mirror image: flush for an ineligible n before its type-2 wrap; hook with `dir=−1` after the type-1 un-wrap with the new field; flush after the loop and before the hopping; assert at the end.

- [ ] **Step 1:** Run `compare_trajectory.sh start standard delayed`. Expected: FAIL, because `Lowrank_init` is never called, so the `Sequential scheme: DELAYED` line is missing from `info`. Record the output.
- [ ] **Step 2:** Implement the changes above.
- [ ] **Step 3:** `regress_vs_baseline.sh`. Expected: all sets `IDENTICAL` (standard bit-identical).
- [ ] **Step 4:** `compare_trajectory.sh` for sets `start`, `doped`, `density`, `continuous` with `standard delayed`, and `start` with `Propose_S0=.T.`. Expected: exit 0 for all; any toggle difference is classified `ROUNDOFF`. Paste the summaries.
- [ ] **Step 5:** Review Focus 4: run the `tv` set with `Sequential_scheme="delayed"` and compare against standard with `compare_runs.py` at tolerance 0. Expected: identical files and the fallback message printed exactly once in stdout.
- [ ] **Step 6:** `ctest` all green. Commit `feat: route sequential updates through delayed scheme`.

### Task 9: Debug mode (`-DALF_LOWRANK_DEBUG`)

**Files:**
- Modify: `Prog/Lowrank_update_mod.F90`

**Interfaces:**
- Produces (debug builds only): a shadow `Gsh(Ndim,Ndim,N_FL)`. It is set to GR whenever `n_pend` goes 0 → 1, gets the same factors as the hook, and gets a standard rank-1 update on each accept. After every hook and every accept it checks `max|Lowrank_current_G − Gsh| / max|Gsh| ≤ TOL_UNIT`; on violation `Terminate_on_error` with `n_op, nt, nf` and the error value.

- [ ] **Step 1:** Build with `ALF_FLAGS_EXT=-DALF_LOWRANK_DEBUG` and add a temporary `-DLOWRANK_DEBUG_SELFTEST` that skips one X-row scaling in the hook. Run `start` with `delayed`: expected termination with the error message. Remove the self-test define (keep the `#ifdef` code path).
- [ ] **Step 2:** Debug build, run `start` and `doped` with `delayed` (`NSweep=3 NBin=2`). Expected: normal completion.
- [ ] **Step 3:** Rebuild normally; `regress_vs_baseline.sh` all `IDENTICAL`. Commit `feat: shadow-G debug mode for low-rank updates`.

### Task 10: Delayed — statistical and stability validation, review (⏸ Checkpoint 3)

- [ ] **Step 1:** In the background: `stat_runs.sh {half,doped,density} delayed 3`. Compare each with the standard seed-1 runs from Task 5 using `stat_compare.py`. Expected: PASS ×3; precision diagnostics comparable (paste the table).
- [ ] **Step 2:** Use superpowers:requesting-code-review on the Phase 3 diff (`f9bb1707..HEAD`); address findings with superpowers:receiving-code-review; commit the fixes.
- [ ] **Step 3:** superpowers:verification-before-completion: rerun CTest, regression, the trajectory set, and show the outputs. Update `CLAUDE.md` (decisions, tolerances). Commit. **Stop for Checkpoint 3.**

## Phase 4 — Submatrix update

### Task 11: Submatrix kernel

**Files:**
- Modify: `Prog/Lowrank_update_mod.F90`; `testsuite/Prog.tests/36-lowrank-unit.F90` (add the `submatrix` candidate)

**Interfaces (same public routines; scheme-dependent bodies, notes §4):**
- `Lowrank_propose`: if `slot(p) ≠ 0` for any flavour, `Lowrank_flush(GR)` first (repeated site). Then `a = Γ̃⁻¹(D·G0(S,p))` (ZGEMV k×k), `G_pp = G0(p,p) + G0(p,S)·a`, cached for accept.
- `Lowrank_accept`: border `Ginv` with `s = d`, `Γ̃⁻¹w = −a`, `qΓ̃⁻¹ = −Δ·(G0(p,S) Γ̃⁻¹)` (one ZGEMV), rank-1 ZGERU on the k×k block; append to `S`, set `slot(p)`, `Dvec`.
- `Lowrank_wrap_hook`: for each `i` with `slot(P(i)) = m > 0`: `Ginv(m,1:k) *= t`, `Ginv(1:k,m) /= t`.
- `Lowrank_flush`: copy `C = G0(:,S)` (N×k); `M = Γ̃⁻¹·D·(E_Sᵀ − G0(S,:))` (ZGEMM); `G0 −= C·M` (ZGEMM); clear `slot`.
- `Lowrank_current_G`: the same algebra into `Gcur` without modifying GR.

- [ ] **Step 1:** Add the `submatrix` candidate to `run_case` (all cases, both directions, all N_delay; the `repeat` case is Review Focus 2). Run: expected FAIL.
- [ ] **Step 2:** Implement.
- [ ] **Step 3:** Run the unit test. Expected: all errors ≤ `TOL_UNIT` for `delayed` and `submatrix`. Add the kernel-equivalence check (Task 7 Step 4) for `submatrix`; expected PASS.
- [ ] **Step 4:** Rebuild normally; `regress_vs_baseline.sh` all `IDENTICAL`. Commit `feat: submatrix update kernel`.

### Task 12: Submatrix — integration, trajectory, debug

- [ ] **Step 1:** Extend the Task 9 shadow-G check to submatrix (it uses `Lowrank_current_G`, so only the self-test needs a submatrix variant: skip one `Ginv` rescale). Expected: the self-test fires; the normal debug runs `start`, `doped` with `submatrix` complete.
- [ ] **Step 2:** `compare_trajectory.sh` for `start`, `doped`, `density`, `continuous` with `standard submatrix` and with `delayed submatrix`; `start` with `Propose_S0=.T.`. Expected: exit 0, only `ROUNDOFF` classifications.
- [ ] **Step 3:** `tv` with `submatrix` vs standard at tolerance 0: identical + one message. `ctest` green. Commit `test: submatrix trajectory validation`.

### Task 13: Submatrix — statistical, stability, review (⏸ Checkpoint 4)

- [ ] **Step 1:** In the background: `stat_runs.sh {half,doped,density} submatrix 4`; `stat_compare.py` against standard seed 1. Expected: PASS ×3; precision table.
- [ ] **Step 2:** superpowers:requesting-code-review on the Phase 4 diff; fix findings; commit.
- [ ] **Step 3:** superpowers:verification-before-completion; update `CLAUDE.md`. **Stop for Checkpoint 4.**

## Phase 5 — Benchmarks and tuning

### Task 14: Benchmark harness and size scan

**Files:**
- Create: `testsuite/lowrank/bench.sh`, `testsuite/lowrank/bench_report.py`

**Interfaces:**
- Produces: `bench.sh SCHEME L N_delay [THREADS]` runs a Mz Hubbard job (U=4, β=2, Dtau=0.1, `Ltau=0`, `LOBS_ST=LOBS_EN=1` to minimise measurement cost, `NSweep` chosen so the run takes ≈ 30–90 s, `NBin=2`) with `OMP_NUM_THREADS=VECLIB_MAXIMUM_THREADS=THREADS` (default 1). It prints `scheme L N N_delay threads sec_per_sweep acceptance`, where sec/sweep comes from the `CPU Time` in `info` divided by the number of sweeps. `bench_report.py` collects the lines into `docs/benchmarks/results.csv` and a markdown table.

- [ ] **Step 1:** Scan L ∈ {6, 8, 10, 12, 14, 16, 20} (N = 36…400) × the three schemes at N_delay=32, single thread, in the background in chunks of ≤ 30 min (ask the user first if the estimated total exceeds 30 min).
- [ ] **Step 2:** If delayed/submatrix are not faster than standard at L ≥ 14, apply superpowers:systematic-debugging (profile first) before reporting.
- [ ] **Step 3:** Commit `bench: add benchmark harness and size scan results`.

### Task 15: N_delay scan, profile, default

- [ ] **Step 1:** At L=12 and L=20: N_delay ∈ {4, 8, 16, 32, 64, 128, 256} for delayed and submatrix. Pick the default (one value or a simple rule in N, e.g. Sun–Xu-style `2^λ`), justified by the data.
- [ ] **Step 2:** Profile the standard and best runs at L=20 with `sample <pid> 20` (macOS). Report the time fractions (flush ZGEMM, ratio/border, hopping wraps, stabilisation, measurement). Record compiler, flags, BLAS (Accelerate) and thread settings.
- [ ] **Step 3:** One multithreaded comparison (THREADS=8) at L=20.
- [ ] **Step 4:** Set the new default `N_delay` in `QMC_runtime_var_mod.F90` and in Start's `parameters`. Rerun `regress_vs_baseline.sh` (standard must stay `IDENTICAL`), the unit test and one trajectory comparison per scheme. Commit `perf: set default N_delay from scan`.

### Task 16: Final report (⏸ Checkpoint 5) and finish

**Files:**
- Create: `docs/benchmarks/report.md`
- Modify: `CLAUDE.md`

- [ ] **Step 1:** Write `report.md`: method, the validation summary (unit/trajectory/statistical/regression/stability, with numbers), time per sweep vs N (table), the N_delay scan, the profile, the recommendation, and caveats (MPI/tempering untested, projector untested, `N_FL_eff < N_FL` untested).
- [ ] **Step 2:** superpowers:requesting-code-review on the whole branch; fix; superpowers:verification-before-completion with the full suite. Commit. **Stop for Checkpoint 5.**
- [ ] **Step 3:** After approval: superpowers:finishing-a-development-branch (keep local; no push, no MR).
