# N_delay scan: delayed and submatrix

Mz Hubbard model on an L x L square lattice, U=4, Beta=2, Dtau=0.1 (Ltrot=20), Ltau=0, LOBS_ST=LOBS_EN=1, NBin=2, default Nwrap=10, N_global=1, same seed for all runs, 1 thread. sec/sweep = `CPU Time` / (NSweep*NBin), wall time of the sweep phase only (setup excluded). Speedup = standard / scheme. Precision = `Precision Green Max` / `Mean` and `Precision Phase Max` from `info` (the Green precision compares G propagated through the sweep with G recomputed at each stabilisation; trajectories are identical for all runs, so differences are pure round-off). The capacity is min(N_delay, vertices per slice[, Ndim for submatrix]) and G is flushed every slice: the effective delay is min(N_delay, accepts per slice). Produced by `testsuite/lowrank/bench.sh` and `testsuite/lowrank/bench_report.py --ndelay`; raw data `raw_ndelay_scan.txt`.

## L=12 (N=144), NSweep=140, acceptance 0.432044, about 62 accepts per slice

Standard (reference): 0.1842 s/sweep; Precision Green mean/max 1.9009612428683341E-011 / 5.9076009323311141E-007; Precision Phase max 8.7707618945387367E-015.

| N_delay | delayed [s/sweep] | speedup | Green max | max / standard | Green mean | Phase max | submatrix [s/sweep] | speedup | Green max | max / standard | Green mean | Phase max |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | 0.1404 | 1.31 | 5.93e-07 | 1.00 | 1.89e-11 | 1.04e-14 | 0.1458 | 1.26 | 5.92e-07 | 1.00 | 1.9e-11 | 1.05e-14 |
| 8 | 0.1267 | 1.45 | 5.93e-07 | 1.00 | 1.89e-11 | 1.01e-14 | 0.1363 | 1.35 | 5.93e-07 | 1.00 | 1.9e-11 | 8.88e-15 |
| 16 | 0.1223 | 1.51 | 5.92e-07 | 1.00 | 1.9e-11 | 9.1e-15 | 0.1169 | 1.58 | 5.92e-07 | 1.00 | 1.91e-11 | 9.88e-15 |
| 32 | 0.1243 | 1.48 | 5.93e-07 | 1.00 | 1.9e-11 | 8.88e-15 | 0.1178 | 1.56 | 5.92e-07 | 1.00 | 1.92e-11 | 9.55e-15 |
| 64 | 0.1323 | 1.39 | 5.94e-07 | 1.00 | 1.91e-11 | 1.03e-14 | 0.1307 | 1.41 | 5.93e-07 | 1.00 | 1.95e-11 | 9.21e-15 |
| 128 | 0.1309 | 1.41 | 5.93e-07 | 1.00 | 1.91e-11 | 8.77e-15 | 0.1345 | 1.37 | 5.93e-07 | 1.00 | 1.95e-11 | 1.07e-14 |
| 256 | 0.1309 | 1.41 | 5.93e-07 | 1.00 | 1.91e-11 | 8.77e-15 | 0.1315 | 1.40 | 5.93e-07 | 1.00 | 1.95e-11 | 1.07e-14 |

## L=20 (N=400), NSweep=12, acceptance 0.430935, about 172 accepts per slice

Standard (reference): 2.592 s/sweep; Precision Green mean/max 8.8184956334657815E-012 / 5.6328444356692738E-007; Precision Phase max 1.4099832412739488E-014.

| N_delay | delayed [s/sweep] | speedup | Green max | max / standard | Green mean | Phase max | submatrix [s/sweep] | speedup | Green max | max / standard | Green mean | Phase max |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | 1.731 | 1.50 | 5.63e-07 | 1.00 | 8.75e-12 | 1.18e-14 | 1.776 | 1.46 | 5.63e-07 | 1.00 | 8.75e-12 | 1.23e-14 |
| 8 | 1.381 | 1.88 | 5.63e-07 | 1.00 | 8.74e-12 | 1.61e-14 | 1.379 | 1.88 | 5.63e-07 | 1.00 | 8.76e-12 | 1.28e-14 |
| 16 | 1.23 | 2.11 | 5.63e-07 | 1.00 | 8.75e-12 | 1.22e-14 | 1.209 | 2.14 | 5.63e-07 | 1.00 | 8.77e-12 | 1.18e-14 |
| 32 | 1.172 | 2.21 | 5.63e-07 | 1.00 | 8.75e-12 | 1.27e-14 | 1.112 | 2.33 | 5.63e-07 | 1.00 | 8.82e-12 | 1.18e-14 |
| 64 | 1.232 | 2.10 | 5.63e-07 | 1.00 | 8.79e-12 | 1.45e-14 | 1.124 | 2.31 | 5.63e-07 | 1.00 | 8.9e-12 | 1.24e-14 |
| 128 | 1.393 | 1.86 | 5.63e-07 | 1.00 | 8.82e-12 | 1.42e-14 | 1.297 | 2.00 | 5.63e-07 | 1.00 | 8.98e-12 | 1.08e-14 |
| 256 | 1.56 | 1.66 | 5.63e-07 | 1.00 | 8.87e-12 | 1.44e-14 | 1.606 | 1.61 | 5.63e-07 | 1.00 | 9.12e-12 | 9.99e-15 |

## Multithreaded comparison (L=20)

| scheme | N_delay | threads | NSweep | s/sweep | speedup vs standard, same threads |
|---|---:|---:|---:|---:|---:|
| standard | 32 | 8 | 12 | 2.596 | 1.00 |
| delayed | 32 | 8 | 12 | 1.185 | 2.19 |
| submatrix | 32 | 8 | 12 | 1.121 | 2.32 |

## Provenance

- Git commit (of the binary, from `info`): 26b354c7
- Compiler: GNU Fortran (Homebrew GCC 16.2.0) 16.2.0
- Flags: `-std=f2008 -cpp -O3 -ffree-line-length-none -ffast-math -fopenmp`
- BLAS/LAPACK: Apple Accelerate (vecLib), linked as `-llapack -lblas`
- Threads: OMP_NUM_THREADS = VECLIB_MAXIMUM_THREADS = 1, 8
- Machine: Apple M1 Pro (8 performance + 2 efficiency cores)
- Date (START TIME in `info`): 2026-10-31

## Default N_delay

**Chosen default: N_delay = 32** (a single value; the optimum does not move with N). At L=12 (62 accepts per slice) the
best points are N_delay=16 and 32 for both schemes, within 1-2 % of each other, which is inside the run-to-run noise (about
2 % between the three standard L=12 runs); at L=20 (172 accepts per slice) 32 is the best for both schemes (delayed 1.172,
submatrix 1.112 s/sweep, N_delay=64 within 1 % for submatrix). Beyond 64 the time rises again (L=20: N_delay=256 is 33 %
slower than 32 for delayed and 44 % slower for submatrix), so a large value is not free in time even though the per-slice
cap limits the memory: the pending-vertex work per proposal grows with the number of pending vertices, so the plateau
expected from the cap is not seen; the cap only prevents N_delay > vertices per slice from mattering. Below 16 the flush
count is too high (N_delay=4 is 1.3-1.5x slower than 32 at L=20).

Precision against N_delay: `Precision Green Max` is 5.9e-7 (L=12) and 5.6e-7 (L=20) for every run of both schemes and for standard
(ratio 1.00 to three digits); it is set by one event of the trajectory (identical for all runs) and carries no information
about N_delay. The mean is more sensitive: the submatrix Green mean grows slightly with N_delay (L=20: 8.75e-12 at 4 to
9.12e-12 at 256, +4 %; L=12: 1.90e-11 to 1.95e-11, +2.6 %), the delayed mean less (8.75e-12 to 8.87e-12, +1.4 %); at the chosen
N_delay=32 both are within 0.1-0.4 % of standard (8.82e-12 vs 8.82e-12 at L=20 for submatrix, 8.75e-12 for delayed).
`Precision Phase Max` scatters between 0.9e-14 and 1.6e-14 without trend (standard: 0.9e-14 at L=12, 1.4e-14 at L=20).
The ratio of 1.9 between submatrix and standard Green precision max seen on the density set at Checkpoint 4 does not
appear for this Mz Hubbard model at any N_delay in {4,...,256}.
