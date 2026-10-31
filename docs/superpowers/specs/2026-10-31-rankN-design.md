# Design: rank-N vertices, general submatrix blocks, projector, cross-slice delay (project 2, Phases 2–5)

Date: 2026-10-31. Formulas: `docs/algorithm-notes.md` Part II (Secs. 9–13). Approval: replaced by independent review
(user instruction for Phases 0–6); review outcome recorded in `CLAUDE.md`.

## 1. Goal and success criteria
Both low-rank schemes accelerate **every** vertex type ALF ships (diagonal and non-diagonal, any `N_non_zero`, field
types 1–4, `g_t`), the projector code, and keep the Markov chain of the standard update. Done means: the five tiers
of project 1 (unit ≤ 1e-10 / 1e-8, identical trajectories, pre-registered statistics, byte-identical default,
stability) pass for tV, Kondo-J, Z2 gauge (type-1 Ising, Δ ≠ 0) and Spin-Peierls (type-4 bonds), Jz (test model),
a projector set and a sign-problem/complex-hopping set — for delayed and submatrix.

## 2. Decisions
| # | Decision | Reason |
|---|---|---|
| D1 | Eligible = every vertex whose support `max_f Op_V(n,f)%N ≤ N_P_LIM = 8`. Larger supports: flush + `Upgrade2` (fallback, reported once in `info`). | All shipped vertices have N ≤ 4; a vertex with N ~ Ndim would make one accept a dense update anyway. |
| D2 | Delayed: Sec. 11 (r columns per accept, hooks after every wrap of an eligible vertex, `X(P,:) ← L X(P,:)`, `Y(P,:) ← R^T Y(P,:)`). | Exact; O(N_P² j) per hook. |
| D3 | Submatrix: M-form (notes Sec. 12): state S, dense M = Γ̃^{-1}Q; proposal 12.2; one accept formula 12.3 (pad active sites, `M ← M̃ + a Δ Mat^{-1} b^T`); pre-wrap padding for a non-diagonal wrap that partially overlaps S; `M ← L_S M R_S` in the hooks. **No flush on repeated sites or overlapping bonds.** (Revised after the independent derivation: the first draft's Γ̃^{-1}+Q form with two accept paths lost accuracy, up to 1.1e-7, with many repeated-site updates; notes 12.6.) | Flush-on-repeat would flush at almost every bond of tV/Kondo (consecutive bonds share sites); disjoint families depend on the Hamiltonian's vertex order, which we must not change (chain identity). The M-form costs O(k² r) per proposal and accept, like the rank-1 bordering, and needs no G0(S,S) gather. |
| D4 | Project 1's rank-1 submatrix state (Ginv, Dvec, flush before a repeated site) is replaced by the M-form for rank-1 as well. | One code path; Hubbard never revisits a site within a slice. Trajectory tier re-run for all sets. |
| D5 | `N_delay` bounds the **pending rank per flavour**: delayed = columns of X; submatrix = k = |S| (padded sites included). Capacity `cap = max(N_delay, N_P_LIM)` (a vertex always fits), delayed also ≤ M·N_P_LIM, submatrix ≤ Ndim. A vertex that would exceed cap flushes first. Rank-1: unchanged semantics. | Matches the papers' `k·n_d = 2^λ` (total pending rank is what sets the BLAS-3 block size). |
| D6 | Ratio: build `Mat` from the current-G block with Upgrade2's loop order and determinant branches (op_dim 1, 2 with s1/s2, `Det` for > 2), duplicated in `Lowrank_update_mod`. `Upgrade2` itself is **not** touched. | Default path must stay byte-identical; duplication of ~20 lines is the price. |
| D7 | `Mat^{-1}` / `Σ^{-1}` (r ≤ 8): LAPACK `ZGETRF`/`ZGETRS` on the r×r matrix. | Pivoted, robust; cost negligible. |
| D8 | Hooks compute L, R from `Op%U`, `Op%E`, `g`/`g_t(τ)` and the wrap's field exactly as `Op_Wrapup`/`Op_Wrapdo` do (table in Sec. 9 of the notes). New call sites in `Wrapgr_mod`: hook after the type-2 wrap of non-diagonal eligible vertices; a pre-wrap call before both wraps of non-diagonal eligible vertices (submatrix padding; no-op for delayed). | ALF's wrap must stay as is; the hooks mirror it on the pending state. |
| D9 | Projector (T = 0): the schemes run on ALF's projector G unchanged (same `WRAPGRUP`/`WRAPGRDO`); validated statistically in Phase 5 (pv_adiab rank-1, tV projector rank-2). **Submatrix-T (arXiv:2602.03656) is not implemented.** | It replaces ALF's G by the N×N_p projector representation and its own propagation: a different code structure, out of scope; [SX-D] reports only 1.3× for delay at T = 0 unless N_p ~ N. |
| D10 | Cross-slice delay: research prototype for the **delayed scheme only**, opt-in key `Lowrank_cross_slice` (logical, default `.F.`) in `&VAR_QMC`. Up/down hopping wraps applied to X and Y (Sec. 12.9). Submatrix with the key on converts to delayed form at the slice end and continues delayed until the next flush. `main.F90` flushes before every read of G (measurement slice, stabilisation, global/Global_tau moves, Langevin/HMC, end of a sweep direction). Keep/opt-in/remove decided in Phase 5 from a Mac measurement. | Gains only when accepts per slice < N_delay and hopping is cheap (checkerboard); at large N the in-slice capacity flushes dominate (audit §3). |
| D11 | Validation model for diagonal rank-2: a test Hamiltonian `LR_Jz_Test` (Mz Hubbard on the square lattice + nearest-neighbour `Predefined_Int_Jz` bonds), added through `Prog/Hamiltonians.list.d/`. Marked test-only. | No shipped Hamiltonian uses Jz or SaSb; SaSb-zz needs the spin-orbital (N_FL = 1, Ndim = 2N) layout and adds nothing that Jz does not test (diagonal multi-site, N_FL = 2 with opposite couplings). Rank-4 diagonal is covered at the unit tier (SaSb-zz operator from `Predefined_Int_SaSb`). |
| D12 | Debug mode (`-DALF_LOWRANK_DEBUG`): the shadow G receives the full rank-r Woodbury update and the general L·R wraps; checks as in project 1. | Catches stale reads and wrong hooks in real runs. |

## 3. Supported cases and fallbacks
| Vertex (model) | N | diag | r | Delayed | Submatrix |
|---|---|---|---|---|---|
| U_SUN / U_MZ (Hubbard, PV, Kondo U, tV U=0, Spin-Peierls U) | 1 | yes | 1 (0 if g = 0) | yes | yes |
| LRC | 1 | yes | 1 | yes (no effect with Global_tau) | yes |
| Jz (test model) | 2 | yes | 2 | yes | yes |
| SaSb-zz (unit tier only) | 4 | yes | 4 | yes | yes |
| V_SUN / VJ_SUN (Kondo J, tV) | 2 | no | 2 | yes | yes (padding) |
| Ising_SUN (Z2 gauge, type 1) | 2 | no | 2 | yes | yes |
| Spin-Peierls type-4 bonds | 2 | no | 2 | yes | yes |
| Z2 site-matter (g = 0, N_non_zero = 1, Δ = 0) | 1 | yes | 1 | yes (Δ = 0 → ratio 1) | yes (Q entry 0) |
| any vertex with N > 8 | > 8 | – | – | fallback | fallback |
Unsupported by design: nothing else. If a formula cannot be validated at the unit tier within three fix attempts,
that vertex class becomes ineligible (fallback) and the case goes to the open questions.

## 4. Code changes (files and interfaces)
* `Prog/Lowrank_update_mod.F90`: `N_P_LIM = 8`; eligibility D1; per-flavour state as now, submatrix `S`, `slot`,
  `M(cap,cap,N_FL_eff)` replacing `Ginv`/`Dvec`; new public `Lowrank_prewrap_hook(n_op, nt, HS_field, kind, dir)`; `Lowrank_wrap_hook(n_op, nt,
  HS_field, kind, dir)` gains `kind` (1 or 2); `Lowrank_propose`/`Lowrank_accept` rank-r; flush with dense M;
  cross-slice entry points `Lowrank_hop_up(nf_eff ...)`/`Lowrank_hop_down` (D10). Existing public names keep their meaning.
* `Prog/Wrapgr_mod.F90`: call the pre-wrap hook before, and the hook after, each wrap of an eligible vertex (type 2
  calls only for non-diagonal ones); slice-end flush skipped when cross-slice is on (then `Lowrank_hop_*`).
* `Prog/main.F90`: with cross-slice on, `Lowrank_flush` before every read of G (D10). Nothing changes when off.
* `Prog/QMC_runtime_var_mod.F90`: `Lowrank_cross_slice` key, getter, broadcast, default `.F.`.
* `testsuite/Prog.tests/36-lowrank-unit.F90`: rank-N cases (Sec. 5); `testsuite/lowrank/sets.sh`: new sets.
* New test Hamiltonian (D11).

## 5. Validation (tiers of `CLAUDE.md`; tolerances unchanged)
1. Unit (ctest 36): delayed and submatrix vs real `Upgrade2` and exact recomputation for: Jz, SaSb-zz, tV bond
   (V_SUN), Kondo VJ_SUN, Ising_SUN type 1, Spin-Peierls-like type-4 bond, a mixed set with repeated sites and partial
   overlaps, a vertex with a zero eigenvalue, `g_t`; up and down; N_delay {1, 2, 7, 32, 1000}; ill-conditioned variant.
2. Trajectory: tv, tv projector, tv with flux (complex hopping), kondo (J ≠ 0), z2 (sequential gauge Ising), spin-peierls,
   jz-test, plus all project-1 sets; standard vs delayed, standard vs submatrix, delayed vs submatrix.
3. Statistical: tv, kondo, jz-test, a projector set (tv projector), doped (sign problem) — pre-registered rules.
4. Regression: `regress_vs_baseline.sh` 9/9 IDENTICAL; ctest all pass.
5. Stability: paired precision diagnostics ≤ standard × (project-1 observed factor ≤ 2).
Debug build: shadow checks over the trajectory sets.

## 6. Risks
* Submatrix forward error grows with cond(Γ̃) (project 1: up to 2.6e-10 in mz_illcond with the explicit Γ̃^{-1});
  the M-form's behaviour there is measured at the unit tier. Mitigation: unit tier at the ill-conditioned case; if > TOL_ILL, cap k for non-diagonal
  vertices or route them to delayed form (decision recorded then).
* Performance of submatrix for tV-like models: padding inflates k (≈ 2–5 sites per accept). Measured in Phase 5;
  if submatrix is slower than delayed there, document it (correctness first).
* Cross-slice: main.F90 read points must all flush; the debug build asserts flushed G at every read.
