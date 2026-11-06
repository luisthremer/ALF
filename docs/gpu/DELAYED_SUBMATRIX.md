# Delayed and submatrix updates (test branch `gpu-port-delayed-submatrix`)

Test branch, numerical acceptance not pursued (owner, 2026-11-06). Source: `fork/feature/delayed-submatrix-clean`
(source files only: `Prog/Lowrank_update_mod.F90`, `Operator_mod`, `QMC_runtime_var_mod`, `upgrade_mod`, `Wrapgr_mod`,
`Timers_mod`, `Trace_updates_mod`, `main.F90` hooks, parameters). Keys: `Sequential_scheme = "standard"|"delayed"|"submatrix"`,
`N_delay` (>= 1 fixed, 0 automatic static rule).

## CPU
Merged with the GPU hooks in `main.F90` / `Wrapgr_mod.F90` (GPU early returns come first, the low-rank code runs on the host
path). GNU, 6x6 Hubbard, 40x10 sweeps: same trajectory as `standard` to ~1e-13, identical acceptance.

## GPU (`Prog/gpu_lowrank_mod.F90`, called from `gpu_update_mod`)
Scope = what the GPU path supports: rank-1 vertices, diagonal operators. Other vertices flush and use the standard update.
- **delayed**: `G = G0 - X Y^T` on the device. Proposal: one single-thread region (2 reals D->H). Accept: ONE kernel per flavour
  (`c`, `r` from `G0`, `X`, `Y`; O(N j)). Wrap: `G0` kernel plus rows of `X`, `Y`. Flush: one `la_gemm` (ZGEMM).
- **submatrix**: M-form with `M` (k x k) on the host. Proposal: gather `G0(p,p)`, `G0(S,p)`, `G0(p,S)` (2k+1 elements), host O(k^2).
  Accept: host only, no device work. Flush: gather + two `la_gemm`.
- Flushes: capacity, before an N > 1 vertex, end of every WRAPGRUP / WRAPGRDO (hopping, stabilization, measurements read a flushed `G_d`).

## Measurements (nividia, A30, nvfortran 26.9, N_FL = 2 Hubbard U=4 beta=5, Ltau=1, one run each)
| lattice, sweeps | scheme | wall s | update launches | accepted / proposals |
|---|---|---:|---:|---|
| 6x6, 400 | standard | 204 | 2 470 712 | 617 678 / 1 440 000 |
| 6x6, 400 | delayed (N_delay 16) | 234 | 1 344 192 | same |
| 6x6, 400 | submatrix (16) | 259 | n/a | same |
| 12x12, 20 | standard / delayed / submatrix | 42.8 / 51.5 / 48.1 | | |
| 20x20, 10 | standard | 79.2 | 690 340 | 172 585 / 400 000 |
| 20x20, 10 | delayed (16) | 87.6 | 367 686 | same |
| 20x20, 10 | submatrix (32) | 92.8 | 35 466 | same |

Acceptance is identical in all schemes; `Pot_scal` of the first bin agrees to 1e-13 (relative); Precision Green max unchanged (~6e-11).

## Reading
The accept/update launches fall by 1.9x (delayed) and 19x (submatrix), but wall time does not improve: the per-site cost that
remains is the wrap launch (2 per site, 800 000 here) plus the ratio round trip (one D->H wait per proposal and flavour). Both
are unchanged by delayed/submatrix and dominate at N <= 400. Next levers: wraps deferred/fused for the pending sites (the wrap of
`G0` is only needed at the flush for the submatrix scheme), and the device-side accept decision (ticket 406).

---

# Speed-up round (branch `gpu-port-delayed-submatrix-speedup`, 2026-11-06)

Goal: remove the bottlenecks of the GPU delayed/submatrix version. Method: `nsys` (launch/copy counts, kernel time) and wall-clock
timers (`GPU time in sequential updates / TAU_M / tau_obs` lines in `info`, plus the `ALF_TIMERS` build). The first profile showed the
GPU almost idle (kernels 1.4 s of 24 s under `nsys`; 432 000 stream synchronisations, 46 000 small D->H copies): the run was bound by
the host-device round trips, not by kernels.

## Bottlenecks found and what was done

| # | Bottleneck | Measured | Fix |
|---|---|---|---|
| 1 | **Per-site device round trips** in the submatrix update: per vertex and flavour a wrap launch, an H->D copy of S, a gather launch, a D->H copy (about 8 synchronous operations per site) | 117 us per proposal, 5.0 s of 7.9 s | **Block mode** (`gpu_lowrank_mod`: `blk_*`): a block of N_delay consecutive vertices is handled on the host. Between flushes G0 changes only by the diagonal wrap scalings, so the host copies the principal block G0(B,B) once, scales it with every wrap, and decides all proposals/accepts from it. At the block end the device applies the scalings to G0 and flushes M (one upload of S, M, factors). About 12 device operations per block instead of ~2 per site: 0.53 s (10x) |
| 2 | Same for **delayed**: wrap, G(p,p), accept were separate synchronous regions per flavour | about 38 us per site | **Fused step** (`k_dly_up/down`): wrap + G(p,p) of all flavours in one launch with one synchronisation; accept of both flavours in one launch; field / wrap work deferred into the next launch's arguments (no H->D copy per accept): update time at 12x12: 1.9 s (fused step, first version) -> 1.65 s |
| 3 | **Host symmetrisation of the measurements** (`Symm = .T.`, the default parameter file): `gpu_obs_setup` refused the device measurement with Symm, so G was copied to the host every slice and the 4 tau matrices every tau slice | 0.8 s of 3.3 s | **`gpu_symm`**: Hop_mod_Symm on the device; observables work on the symmetrised matrices. Symm(M) = A M B with fixed A, B built once from the hop elements: 2 GEMMs instead of ~14 launches per matrix |
| 4 | **Launch-bound hop products** (`mmthr`, `mmthl_m1`, ... one launch per checkerboard group, ~8 per product, 4 per matrix per slice) | tau_slice 0.17 s, hops in every slice | **Dense hop matrices** (`hopd_*`): each sequence is a fixed linear map, built once by applying the kernels to the identity; one GEMM + copy (Ndim <= 800) |
| 5 | tau observables kernels: 144 threads each with a serial sum over 144 terms | small | nested `parallel do reduction` over the sum (no measurable gain at N = 144; helps larger N) |

nvfortran 26.9 note: `blk_end` is mis-compiled at -O2/-O3 (wrong G after the first stabilisation, found by bisecting the optimisation
level per routine: -O1 and -O0 correct). It is pinned to -O1 with `!pgi$r opt=1`. The cause is not understood; the routine is a few
host loops around target regions and `target update`s.

## Results (nividia, A30; Hubbard U = 4, beta = 5, N_FL = 2, Ltau = 1, Symm = .T., default Nwrap; one run each)

Same inputs, same seeds; acceptance is identical in every run and `Pot_scal` agrees to 12 digits.

| 20x20, 6 sweeps (NSweep 3, NBin 2) | before (branch gpu-port-delayed-submatrix) | after (-speedup) | speed-up |
|---|---:|---:|---:|
| standard | 42.2 s | 25.5 s | 1.7x |
| delayed (N_delay auto) | 46.5 s | 15.5 s | 3.0x |
| submatrix (N_delay auto, 32) | 49.7 s | **8.7 s** | **5.7x** |

12x12, 3 sweeps: submatrix 7.9 s -> 1.9 s (4.1x), standard 6.5 s -> 5.0 s, delayed now 3.2 s. Sequential-update time alone at 12x12 (3 sweeps):
submatrix 5.0 s -> 0.39 s (13x), standard 3.7 s -> 3.3 s, delayed now 1.65 s.

Overall the submatrix run is 4.9x faster than the *old standard* run and 2.9x faster than the *new standard* run.

## Not done (next candidates)
- Delayed per-site latency (about 38 us/site: one launch + one D->H wait + one accept launch). A block mode like the submatrix one would
  need either the M-form flush (identical to submatrix) or a replay of the accepted columns on the device.
- Stabilisation: cuSOLVER QR of small matrices is ~1 ms each (302a/306), CGR2_2 in the tau pass (0.17 s of 2.0 s at 12x12).
- Fixed costs: ~0.6 s start-up (mirror, first stabilisation) and the per-bin output (0.1 s at 12x12, 0.7 s at 20x20 per bin; amortised
  over more sweeps per bin in real runs).
- The N > 1 / non-diagonal vertices still take the standard path (flush first), as before.
