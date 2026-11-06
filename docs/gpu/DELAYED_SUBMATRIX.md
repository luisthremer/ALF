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
