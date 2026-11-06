# Not implemented yet, by group

State 2026-11-06, branch `gpu-port-delayed-submatrix`. Numbers in brackets are the global priority numbers of `PARKING.md`
(there: column `x`). "Parked" = owner decision of 2026-11-06, see the section "Parked" in `PARKING.md`.

## A Correctness
Nothing open.

## B Model and feature coverage
- Rank > 1 and non-diagonal vertices (#15): would open tV, Kondo, Spin_Peierls
- Flavour-symmetry reconstruction `Calc_Fl` (#16): would open Plain_Vanilla, Spin_Peierls
- Device observables for other models, entanglement, custom observables (#19)
- MPI together with GPU, HDF5 with GPU untested (#20, parked)
- HS field types 3 and 4 (#29)
- Global-in-τ and global space-time moves (#30)
- `Propose_S0` (#31)
- Time-dependent couplings and hoppings (#33)
- STAB1/2/3 and QRREF stabilization (#34)
- Cross-slice delay on the GPU (#39)
- Parallel tempering (#41)
- Projector, T = 0 (#42)
- Langevin and HMC (#44)

## C Performance at small and medium N
- Several chains per GPU (#9, parked)
- Stabilization latency (#12)
- Block mode for the delayed scheme (#17)
- Per-bin output `Pr_obs` (#22)
- Remaining launch overhead (#28)
- Device-side accept decision (#32)
- Start-up (#36)

## D Performance at large N
- Sparse `Symm` (#10)
- Memory footprint, about 1650·N² bytes (#25)
- τ-pass `CGR2_2` block structure (#40)
- Cholesky-QR variants (#43)

## E Intel and portability
- Helma and SuperMUC builds (#14): no access, NVHPC missing on julia2
- Unified accuracy gates, the pre-registered criteria predate the new schemes (#21)

## F Housekeeping
- Validation suite 209, final benchmarks 210, owner gate for the new code (#24)

## Open but with code in place (no `x`)
#11 `N_delay` at large N, #13 Ponte Vecchio, #18 56×56 run, #23 remove instrumentation, #26 parameter choices (owner's decision),
#27 B580 FP64, #35 dense-hop threshold, #37 documentation, #38 ifx pragma.
