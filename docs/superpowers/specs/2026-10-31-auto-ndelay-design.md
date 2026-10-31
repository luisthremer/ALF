# Design: automatic N_delay (project 2, Phase 6)

Date: 2026-10-31. Approval replaced by independent review (user instruction for Phases 0–6); outcome in `CLAUDE.md`.

## 1. Goal and success criteria
`N_delay` (the pending rank per flavour before a flush, spec D5) is chosen automatically for the delayed and submatrix
schemes, with a fixed value as override. Done means: on the Mac the automatic value is within the timing noise of the best
fixed value of a scan (Mz Hubbard L = 12, 16, 20; tV L = 8, 12) for both schemes; the chain with automatic N_delay makes
the same Metropolis decisions as a fixed-N_delay run except at roundoff near-ties (trajectory tier); the chosen value is
logged in `info`; default (standard) path byte-identical.

## 2. Facts this rests on
* N_delay changes only the summation order of the exact update (notes Secs. 3, 11, 12), never the ratio formula; so any
  N_delay — also one that changes during the run — gives the same chain up to roundoff near-ties.
* Papers' rules (notes Sec. 8): delay [SX-D] `n_d = 2^λ, λ = min(6, ⌊log2(N/20)⌋)`; submatrix [SX-S] `k·n_d = 2^⌊log2(N/12)⌋`.
  Mac scan (project 1, Mz Hubbard, 1 thread): best 16 (N = 144) and 32 (N = 400) for both schemes; N_delay = 4 is
  +48/+60 %, 256 +33/+44 % slower than 32 at N = 400 — the rules give 8/16 (delay) and 8/32 (submatrix) there, i.e. up to a
  factor 2 off: a static rule alone is not enough, a cheap measurement must correct it.
* G is flushed at every slice end, so a pending rank beyond the pending columns one slice can produce (`M·r̄·acceptance`)
  is never reached; submatrix bond blocks are capped at K_ND = 8 anyway (D15).

## 3. Decisions
| # | Decision | Reason |
|---|---|---|
| A1 | `N_delay = 0` means automatic and becomes the default; `N_delay ≥ 1` is a fixed override (unchanged semantics); negative values are rejected at the parameter read. | The prompt asks for automatic choice with a fixed override; 0 is otherwise meaningless. `standard` ignores N_delay (default path unchanged). |
| A2 | Static start value `n0 = 2^round(log2(Ndim/20))` (delayed) and `2^floor(log2(Ndim/12))` (submatrix), clamped to `[16, 256]` and to `max_pend = size(Op_V,1)·r_max` (the most columns one slice can produce). `Ndim` (matrix dimension, which is what the BLAS-3 block acts on) replaces the papers' N. | Papers' rules as the prior; lower clamp 16 because our scan never favours < 16 and the lowest candidate n0/2 must not fall below 8 (= N_P_LIM, one rank-8 vertex); upper clamp 256 = the papers' flat region; 2^k keeps BLAS block sizes cache-friendly. |
| A3 | Warm-up refinement: during the first `N_tune = 2` sweeps (up and down), the slices are assigned round-robin to the candidates `{n0/2, n0, 2·n0}` (each clamped as A2): before each WRAPGRUP/WRAPGRDO the candidate of that slice is set (G is flushed there), the wall time of the call (`system_clock`, as ALF's CPU time) and its number of proposals are accumulated per candidate. Then each rank computes seconds per proposal per candidate; the group's rank 0 (`Group_Comm`) chooses the fastest, **switching away from n0 only if it is > 5 % faster than n0**, and broadcasts the choice in `Group_Comm` at the end of sweep `N_tune` (a point every rank of the group reaches). Single runs (noMPI) decide locally. | Interleaving cancels drift; per-proposal time removes the dependence on how many vertices a slice proposes; the 5 % rule (from the prompt) avoids flip-flopping on noise (our repeat noise: 2.3 %); per group because PARALLEL_PARAMS groups may simulate different sizes. |
| A4 | No second round, no re-tuning later. | Timing is stationary after the first sweeps; a single decision keeps the log simple. |
| A5 | Logging: `info` gets `Auto N_delay: static n0 = <n0> (Ndim = <Ndim>), candidates <a> <b> <c>: <t_a> <t_b> <t_c> s/proposal (rank 0 of group), chosen <n>` (one line, group rank 0). With a fixed N_delay nothing new is written. `data.h5`: deferred to the cluster phase (no local HDF5 build allowed); listed in Phase 9. | Provenance for the benchmark campaigns. |
| A6 | Interface: `Lowrank_set_N_delay(n)` (requires a flushed state, re-derives capacities, reallocates if the capacity grows), `Lowrank_autotune_slice_begin()` / `Lowrank_autotune_slice_end(n_prop)` called by WRAPGRUP/WRAPGRDO only while tuning, `Lowrank_autotune_finish()` called by main.F90 after sweep `N_tune` (contains the broadcast). Allocation at init for `2·n0` (max candidate) so tuning never reallocates. | Keeps the tuning logic in the low-rank module; main.F90 gets one guarded call. |
| A7 | Cross-slice key (D10), if kept: tuning runs with it as configured (it changes the timing, not the chain). | Consistency. |

## 4. Validation
1. Unit: `Lowrank_set_N_delay` with pending state → terminates; with flushed state → capacities as for a fresh init
   (ctest 36 extension, a few rows at N_delay {8, 16, 64}).
2. Trajectory: `compare_trajectory.sh <set> standard <scheme> N_delay=0` for start, tv, kondo, jz, pv_adiab, both schemes:
   IDENTICAL or DIVERGED-ROUNDOFF with margin (the tuning changes N_delay between slices).
3. Mac timing: for Mz Hubbard L = 12, 16, 20 and tV L = 8, 12 (1 thread), scan fixed N_delay ∈ {8, 16, 32, 64, 128} and
   run N_delay = 0, 3 repetitions each (`bench.sh`; Task 6a); pass if the automatic run's s/sweep is within
   max(5 %, 2 × spread) of the best fixed value. Report the chosen values.
4. MPI: one 2-rank run with N_delay = 0: both ranks use the same value (info line once, from group rank 0).
5. Regress 13/13 IDENTICAL; ctest all.

## 5. Risks
* Timing noise on a loaded machine picks a worse value: bounded by the 5 % rule and the factor-2 candidate spacing; the
  fixed override remains.
* Warm-up sweeps are part of the production chain (no extra sweeps): correct because the chain is N_delay-independent;
  the measurement is per proposal, so the slice mix does not bias it.
