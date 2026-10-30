# Design: delayed and submatrix local updates in ALF

Date: 2026-10-30 · Branch: `feature/delayed-submatrix` · Baseline: `f81b5cbc`
Math: `docs/algorithm-notes.md` (cited below as "notes §N").

## 1. Goal and success criteria (restated)

Add two faster variants of the sequential single-field update to ALF, for production runs of the built-in
Hubbard model, selected at runtime:

* `delayed`: accepted rank-1 corrections are accumulated as `G = G0 − X Yᵀ` and applied with ZGEMM.
* `submatrix`: the DQMC submatrix update of Sun & Xu [SX-S], storing only `Γ̃⁻¹` (k×k) between flushes.
* `standard` (the default) is bit-identical to the unmodified code.

Done means: unit, trajectory, statistical, regression and stability checks pass (user spec,
"Validation"), plus a benchmark and tuning report. Priorities: chain correctness > default unchanged > speed > elegance.

What you specified: the three-way runtime switch; flushing before anything else touches G; ratios computed exactly as
Upgrade2 computes them, with the same RNG order; rank-1 only, with an automatic fallback to standard and a one-time message;
data layout ready for a later Woodbury extension; BLAS-3 flushes; minimal footprint; global/Langevin/MPI out of scope.
My additions are marked **[proposal]**.

## 2. Approaches considered

**A. Hooks in the existing per-slice loop + new module (recommended).** `WRAPGRUP`/`WRAPGRDO` keep their loop:
`flip` → `S0` → `T0_proposal` ranf → update. With a non-standard scheme, the `Upgrade2` call is routed to a
low-rank kernel. After each per-vertex wrap a cheap hook transforms the pending factors, and a flush runs before
any non-eligible vertex, before Global_tau moves and before returning. Pros: the proposal and RNG logic is shared
by construction, which is the strongest guarantee for the trajectory test; the diff is small. Cons: a few `if`
branches inside the standard loop, so bit-identity must be confirmed by test (it does not change arithmetic).

**B. Separate copies of WRAPGRUP/WRAPGRDO, with dispatch at the top.** The standard loop is untouched textually,
but the proposal logic is duplicated, which risks silent RNG-order drift whenever upstream edits the loop.

**C. A "G object" abstraction used everywhere.** This is the most elegant, but it is invasive: Upgrade2, Wrapgr,
Global and Langevin would all change. It violates the minimal-footprint rule.

→ **A.**

## 3. Architecture

### 3.1 Files

| File | Change |
|---|---|
| `Prog/Lowrank_update_mod.F90` (new) | State, eligibility, ratio/accept kernel, wrap hooks, flush, debug checks |
| `Prog/Wrapgr_mod.F90` | Hooks in `WRAPGRUP`, `WRAPGRDO`; init/finalize in `Wrapgr_alloc/dealloc` |
| `Prog/QMC_runtime_var_mod.F90` | New `&VAR_QMC` entries, getters, MPI broadcast line |
| `Prog/Makefile` (+ deps) | New object |
| `Prog/upgrade_mod.F90` | Only an `#ifdef ALF_TRACE_UPDATES` trace write (off in default builds) **[proposal]** |
| `Scripts_and_Parameters_files/Start/parameters` | Document the new keys (commented default) |
| `testsuite/Prog.tests/…`, `testsuite/lowrank/…` (new) | Unit test, trajectory/statistical scripts |

`main.F90` and the Hamiltonian files are not touched.

### 3.2 Parameters **[proposal: names]**

```
&VAR_QMC
Sequential_scheme = "standard"   ! "standard" | "delayed" | "submatrix"  (case-insensitive)
N_delay           = 32           ! max pending accepted vertices per flavour before a flush (tuned in Phase 5)
/
```
Both are optional, so old parameter files keep working. Invalid values cause an error termination. The scheme
and N_delay are written to `info` only when the scheme is not standard.

### 3.3 Eligibility and fallback

Vertex `n` is eligible if, for every flavour, `Op_V(n,nf)%N_non_zero ≤ 1` **and** `Op_V(n,nf)%diag`
**[proposal — see open question 1]**. Eligibility is computed once at init; if any vertex is ineligible, one
message is printed (count and first index). Ineligible vertices are flushed and then updated with the standard
`Upgrade2`. The Hubbard model in every variant (Mz, SU(N) density channel, continuous HS, bilayer) consists of
eligible vertices only.

### 3.4 State (per computed flavour, allocated at init with `N_delay` columns)

* Common: `n_pend` (pending accepted vertices), scheme, `N_delay`, `eligible(:)`.
* Delayed: `X(Ndim, N_delay·N_P)`, `Y(Ndim, N_delay·N_P)` with `N_P = 1` now. Columns are counted in vectors,
  so a rank-N_P Woodbury extension only changes the column bookkeeping.
* Submatrix: `S(N_delay)` site list, `slot(Ndim)` (site → position, 0 if absent), `Dvec(N_delay)`,
  `Ginv(N_delay, N_delay)` (= `Γ̃⁻¹`), work vectors.

### 3.5 Control flow (up-sweep; the down-sweep mirrors notes §0)

```
WRAPGRUP: hopping (unchanged)
  do n = start, end
    if lowrank and .not. eligible(n): flush
    Op_Wrapup(type 1, old φ)  ;  if lowrank and eligible(n): hook_wrap(n, +old φ)   ! X,Y rows / Γ̃⁻¹ rescale
    flip, S0, T0 ranf  (unchanged)
    if accepted-for-proposal:
       if lowrank and eligible(n): Upgrade_lowrank(...)   ! same scalars + RNG draw as Upgrade2
       else:                         Upgrade2(...)        ! unchanged
    Op_Wrapup(type 2)  (no-op for diagonal ops)
  end do
  if lowrank: flush            ! before Global_tau moves and before return
  assert_flushed               ! debug builds: hard check
```
Submatrix extra rule: flush before proposing at a site already in `S` (notes §4, restriction 2).

### 3.6 Low-rank kernel (`Upgrade_lowrank`)

1. For each computed flavour, form `G_pp` from G0 plus the pending state (notes §3/§4).
2. Evaluate the scalar ratio, the Metropolis weight, the ranf draw and the phase update with the **same expressions in
   the same order** as `Upgrade2` (FINAL mode), including `weight_reconstruction` and the `Control_upgrade*` calls.
3. On accept: append the column/row (delayed) or border `Γ̃⁻¹` (submatrix), set `nsigma%f`, and flush if `n_pend = N_delay`.

The ratio code is duplicated from Upgrade2, not refactored, so the standard path keeps its arithmetic. A unit
test guards the equivalence.

### 3.7 Flush (BLAS-3)

* Delayed: `ZGEMM('N','T', Ndim, Ndim, j, −1, X, Y, 1, G)`.
* Submatrix: `M = Γ̃⁻¹·D·(E_Sᵀ − G0(S,:))` (ZGEMM k×k×N), then `G0 −= G0(:,S)·M` (ZGEMM N×k×N, with `G0(:,S)` copied first).

### 3.8 Stale-G protection

* By construction, pending state exists only inside one `WRAPGRUP`/`WRAPGRDO` call. Every reader outside
  (stabilization, CGR, measurements, time-displaced, global moves, end of sweep) runs after the return, which
  flushes. Inside, the only readers are `Op_Wrap*` (correct by notes §3/§4) and standard `Upgrade2` /
  `Wrapgr_Random_update` (preceded by a flush).
* `lowrank_assert_flushed()` runs at both routines' exits (always on; an integer compare) and terminates with an
  error if anything is pending.
* Debug build (`-DALF_LOWRANK_DEBUG`): a shadow G with standard rank-1 updates runs alongside. After every
  vertex, `max|G0 − XYᵀ − G_shadow|/max|G_shadow|` (or the submatrix equivalent) is checked against the
  tolerance; on violation the run aborts with n, τ, flavour.

## 4. Validation (Phase 2 builds this before any algorithm code)

1. **Unit** (`testsuite/Prog.tests/36-lowrank-unit.F90`, CTest). A synthetic state: Ndim ≈ 64–256,
   Hubbard-like eligible vertices (Mz pair with N_FL=2 type 2; SU(2) density with complex g, N_FL=1, N_SUN=2;
   type 3 continuous), `G = (1+A)⁻¹` with realistic conditioning, several slices up and down with the
   real `Op_Wrapup/do`. The reference is the real `Upgrade2`. Decisions are forced from a fixed accept mask, so they do
   not depend on roundoff; reference ratios come from `Upgrade2` in INTERMEDIATE mode on copies. Compared: every ratio
   (relative) and G after each flush and at slice end (`max|ΔG|/max|G|`). The tolerance is set in Phase 2 from the
   measured error growth with N_delay and Ndim (target 1e-10, justified there).
2. **Trajectory** (`testsuite/lowrank/compare_trajectory.sh` + python). Same seeds and parameters, a few sweeps,
   with standard vs delayed vs submatrix; compare `confout_0` field-by-field (exact) and all observable bins (relative
   ≤ 1e-8 **[proposal]**). A trace build (`-DALF_TRACE_UPDATES`) logs `(τ, n, R, ranf, accepted)` per proposal for
   both paths. The analyser reports the first differing decision and `|Weight − ranf|/Weight`: within ~1e-10 counts
   as roundoff (acceptable), otherwise it is a bug.
3. **Statistical** (background runs of minutes each, analysed with `ana.out`):
   (a) Mz, half filling, 6×6, U=4, β=5;
   (b) Mz, doped (μ chosen in Phase 2 so the average sign is ~0.3–0.8), 4×4 or 6×6;
   (c) density channel (`Mz=.F.`, complex g, N_FL=1, N_SUN=2), half filling.
   Compare energy, potential energy (double occupancy), kinetic, SpinZ_eq (k-resolved), Green_eq, and the average sign/phase
   with z-scores. Pass if all |z| < 3 and the z distribution is consistent with N(0,1).
4. **Regression.** CTest is all green. Bit-identity: build the baseline `f81b5cbc` in a separate worktree
   (`../ALF_baseline`), then run the Start example plus (b) and (c), short, with the same seeds; `cmp` all observable
   and `confout` files against the branch with default parameters.
5. **Stability.** For the statistical runs, compare `Precision Green Mean/Max` and `Precision Phase` from `info`.
   Acceptance: no systematic increase (the maxima within a factor of a few, means comparable).

## 5. Benchmarks (Phase 5)

L = 6…20 (N = 36…400) for the square-lattice Mz Hubbard model, short β (time per sweep ∝ Ltrot); OMP and
vecLib threads pinned to 1, plus a multithreaded run. Report the time per sweep for the three schemes, an N_delay scan,
profiles (`sample`/xctrace), compiler, flags and BLAS. Optionally compare OpenBLAS with Accelerate (already
installed; only link flags change). Default N_delay is chosen from the scan.

## 6. Out of scope

Global updates, Langevin/HMC, MPI/tempering (only a broadcast line is added, untested), projector-specific
validation, rank-N_P > 1, and non-diagonal operators **[pending open question 1]**.

## 7. Open questions for Checkpoint 1

1. **Scope of "rank-1".** Proposal: only *diagonal* rank-1 vertices (`Op%diag`) are eligible for both schemes. That
   covers all Hubbard variants; non-diagonal rank-1 vertices fall back to standard. Delayed could support non-diagonal
   rank-1 by applying the operator's `U` rotations to X and Y (notes §3), at extra code and test cost. Submatrix
   cannot do this cheaply.
2. **Parameter names** `Sequential_scheme` / `N_delay`. OK?
3. **Trace instrumentation in `upgrade_mod.F90`** under `#ifdef ALF_TRACE_UPDATES`. This is the only way to log the
   standard path's ratio and ranf without changing default builds. OK?
4. **Trajectory tolerance** for observable bins: 1e-8 relative is my proposal; Phase 2 will measure it.
