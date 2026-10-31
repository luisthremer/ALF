# Design: automatic N_delay (project 2, Phase 6)

Date: 2026-10-31. Approval replaced by independent review (user instruction for Phases 0–6); outcome in `CLAUDE.md`.
Revision 2 (same day): review folded in (I1–I6, M1–M5, validation gaps); the default is now deterministic.

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
  factor 2 off before clamping; with the clamps of A2 the static value is within 5 % of the best on the Mac, so the
  timing is insurance for other machines and caches.
* G is flushed at every slice end, so a pending rank beyond the pending columns one slice can produce (`M·r̄·acceptance`)
  is never reached; submatrix bond blocks are capped at K_ND = 8 anyway (D15).

## 3. Decisions
| # | Decision | Reason |
|---|---|---|
| A1 | `N_delay = 0` (new default): automatic, **static rule only** (deterministic). `N_delay = -1`: automatic with warm-up timing (opt-in). `N_delay ≥ 1`: fixed (unchanged). `N_delay < -1`: rejected at the parameter read. `info` says `N_delay = auto (static n0 = …)` / `auto+tune (…)`, never "0". | Review I1: timing-dependent choices break bit-reproducibility ("same seed ⇒ same trajectory", restarts, user regression checks); the user asked for timing as *optional*. |
| A2 | Static value: delayed `2^round(log2(Ndim/20))`, submatrix `2^floor(log2(Ndim/12))` (papers' rules with Ndim, the dimension the BLAS-3 block acts on), then clamped to `[16, U]` with U = 128 (delayed; the paper caps n_d at 64, our L=20 data show 256 at +33 %) or 256 (submatrix), then to `max_pend` = Σ over the sequential vertex range of `r_max` (one slice's maximal pending rank) and, for submatrix, to `Ndim`; finally rounded **down** to a power of two. If `max_pend < 16`, the value is `max_pend` rounded down (≥ 1). | Prior from the papers; clamps from our Mac scan (never < 16 best; 4 and 256 clearly worse) and from the slice-end flush. Mac check: n0 = 16/16 (L=12) and 16/32 (L=20, delayed/submatrix), within 5 % of the best (worst +4.9 %, delayed L=20). |
| A3 | Timing (only `N_delay = -1`): candidates `{n0/2, n0, 2·n0}` (each clamped as A2, duplicates dropped; < 2 distinct → no timing). Sweep 1 is not timed (allocation, first-touch, equilibration); sweeps 2 and 3 assign slices round-robin to the candidates; per call of WRAPGRUP/WRAPGRDO the wall time (`system_clock` with int64 counts) and the number of `Upgrade_lowrank` proposals are added to the slice's candidate. After sweep 3 the sums are `MPI_Allreduce`d over `Group_Comm`; the fastest seconds-per-proposal wins only if it beats n0 by > 5 %; zero proposals for any candidate → n0. Every rank computes the same decision from the reduced sums (no separate broadcast needed). | Review I2–I4, M3: robust, no 0/0, no first-sweep bias, ms clock too coarse at L=12. |
| A4 | One round, no re-tuning (cannot reach an optimum 4× away from n0; logged). | Simplicity; A2 already lands within 5 % on the Mac. |
| A5 | Logging: one line in `info` (`Temp_<g>/info` in tempering builds) by group rank 0: static n0, Ndim, candidates, seconds per proposal, chosen value (or "tuning incomplete"/"no proposals"/"cross-slice: no timing"). `data.h5`: the chosen N_delay and the scheme as attributes, written behind `#ifdef HDF5` now; compiled and validated on the clusters (no local HDF5). | Provenance; review I6. |
| A6 | Interface: `Lowrank_set_N_delay(n)` (flushed state required; capacities re-derived; arrays allocated at init for the largest candidate, so no reallocation), `Lowrank_autotune_slice_begin()` / `_end()` in WRAPGRUP/WRAPGRDO (no-op unless tuning), `Lowrank_autotune_finish(isweep_global)` called by main.F90 after every sweep on every rank, unconditionally; it acts only once, when the global sweep counter (across bins) reaches 3, and decides only from parameters and reduced sums. If the run ends earlier: n0 kept, "tuning incomplete" logged. Sequential = .F. or Langevin (no Lowrank_init): finish is a no-op on all ranks. | Review I2: no rank can skip the collective; decisions independent of local state. |
| A7 | With the cross-slice key on (if kept), timing is skipped (static n0). | Review I5: the slice-end flush that the timing relies on is skipped there. |
| A8 | Submatrix bond models: with K_ND = 1 (spec D15) all candidates behave alike for non-diagonal blocks; logged, harmless under the 5 % rule. | Review M4. |
Compatibility: `validate_sequential_scheme` and `Lowrank_init` accept 0 and −1; `Start/parameters` gets `N_delay = 0` with a comment; `bench.sh`/`bench_report.py` label N_delay ≤ 0 as `auto`/`tune` and record the chosen value; `tiers.sh`, `stat_runs.sh` and ctest keep explicit fixed values where pairing needs them. CLAUDE.md's Phase-1 note ("N_delay = 0 → error") is superseded.

## 4. Validation
1. Unit: `Lowrank_set_N_delay` with pending state → terminates; set down and back up within the allocation; static rule
   values for a table of (scheme, Ndim, max_pend) incl. max_pend < 16.
2. Trajectory: `compare_trajectory.sh <set> standard <scheme> N_delay=0` and `... N_delay=-1` for start, tv, kondo, jz,
   pv_adiab, both schemes; and `<scheme> N_delay=<n0>` vs `<scheme> N_delay=-1` (same scheme): IDENTICAL or ROUNDOFF.
3. Edge cases: NSweep=1 NBin=1 with N_delay=-1 ("tuning incomplete"), Langevin, Nt_sequential_end < start, all vertices
   ineligible (support > 8), a projector set, CPU_MAX truncation during tuning: no crash, sensible log line.
4. MPI (local builds): 2 ranks and tempering with 2 groups, N_delay = -1: identical decision on all ranks of a group,
   one info line per group.
5. Mac timing: Mz Hubbard L = 12, 16, 20 and tV L = 8, 12 (1 thread): fixed N_delay ∈ {8, 16, 32, 64, 128} vs N_delay = 0
   and −1, 3 repetitions each; pass if within max(5 %, 2 × spread) of the best fixed; report how often −1's choice differs
   between repetitions.
6. Regress 13/13 IDENTICAL; ctest all.

## 5. Risks
* Timing noise (only with N_delay = −1) picks a worse value: bounded by the 5 % rule and the factor-2 spacing; the default
  is deterministic and the fixed override remains.
* Warm-up sweeps are part of the production chain (no extra sweeps): correct because the chain is N_delay-independent;
  the measurement is per proposal, so the slice mix does not bias it.
