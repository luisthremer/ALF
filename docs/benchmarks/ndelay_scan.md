# N_delay scan: delayed and submatrix

Mz Hubbard model on an L x L square lattice, U=4, Beta=2, Dtau=0.1 (Ltrot=20), Ltau=0, LOBS_ST=LOBS_EN=1, NBin=2, default Nwrap=10, N_global=1, same seed for all runs, 1 thread (except the multithreaded table). sec/sweep = `CPU Time` / (NSweep*NBin), wall time of the sweep phase only (setup excluded). Speedup = standard / scheme. Precision = `Precision Green Max` / `Mean` and `Precision Phase Max` from `info` (the Green precision compares G propagated through the sweep with G recomputed at each stabilisation; trajectories are identical for all runs, so differences are pure round-off). The capacity is min(N_delay, vertices per slice[, Ndim for submatrix]) and G is flushed every slice: the effective delay is min(N_delay, accepts per slice). Produced by `testsuite/lowrank/bench.sh` and `testsuite/lowrank/bench_report.py --ndelay`; raw data `raw_ndelay_scan.txt`, `raw_threads.txt`. Regeneration: see the docstring of `bench_report.py`.

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

## Default N_delay

**Default N_delay = 32** (scanned at N = 144 and 400). The lowest time is at N_delay = (L=12 delayed 16; L=12 submatrix 16; L=20 delayed 32; L=20 submatrix 32); the points within the run-to-run noise (2.3 % spread of 3 runs of standard at L=12) of the best are: L=12 delayed 16 and 32; L=12 submatrix 16 and 32; L=20 delayed 32; L=20 submatrix 32 and 64. So 32 is best or within noise of the best for both schemes at both sizes. It was chosen as the single default for N = 144-400 and is untested outside that range.

Precision check: the submatrix Green mean rises mildly with N_delay (+2.6 % at L=12, +4.2 % at L=20 from N_delay=4 to 256), so precision mildly favours the lower end; at N_delay=32 the Green mean differs from standard by +1.2 % at L=12, 0.0 % at L=20 (delayed: +0.1 % at L=12, -0.7 % at L=20), so there is no visible degradation at the chosen value. `Precision Green Max` of every run is within 0.5 % of standard, presumably set by one event of the (identical) trajectory, so it does not discriminate between N_delay values.

Time relative to N_delay=32 (s/sweep at N_delay / s/sweep at 32 - 1; the numbers behind the choice):

| L | scheme | 4 | 8 | 16 | 32 | 64 | 128 | 256 |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 12 | delayed | +12.9 % | +1.9 % | -1.6 % | 0.0 % | +6.5 % | +5.3 % | +5.3 % |
| 12 | submatrix | +23.8 % | +15.7 % | -0.8 % | 0.0 % | +10.9 % | +14.2 % | +11.6 % |
| 20 | delayed | +47.6 % | +17.8 % | +4.9 % | 0.0 % | +5.1 % | +18.8 % | +33.0 % |
| 20 | submatrix | +59.7 % | +24.0 % | +8.7 % | 0.0 % | +1.1 % | +16.6 % | +44.4 % |

Interpretation (not measured; no profile against N_delay was taken): beyond the optimum the time rises, which we attribute to costs that grow with the number k of pending vertices outgrowing the gain in flush efficiency. For delayed, forming the new column G(:,p) = G0(:,p) - X Y(p,:)^T is O(N k) per accept (BLAS-2); summed over a slice this grows with k, while the BLAS-3 flush total does not. For submatrix, the per-proposal solve with the Gamma-tilde^-1 matrix is O(k^2) and the border update per accept is also O(k^2).

At L=12 (about 62 accepts per slice) N_delay = 64, 128, 256 are all above the cap; N_delay=128 and 256 give 0.1309 and 0.1309 s/sweep (delayed) and 0.1345 and 0.1315 (submatrix), flat to within 2.3 % (noise 2.3 %): this is the plateau, the effective delay is capped by the accepts per slice. At L=20 (about 172 accepts per slice) only N_delay=256 is above the cap, so the grid has no plateau pair there and the rise from 32 to the next values happens below the cap.

## Provenance

- Git commit (of the binary, from `info`): 26b354c7
- Compiler: GNU Fortran (Homebrew GCC 16.2.0) 16.2.0
- Flags: `-std=f2008 -cpp -O3 -ffree-line-length-none -ffast-math -fopenmp`
- BLAS/LAPACK: Apple Accelerate (vecLib), linked as `-llapack -lblas`
- Threads: OMP_NUM_THREADS = VECLIB_MAXIMUM_THREADS = 1, 8
- Machine: Apple M1 Pro (8 performance + 2 efficiency cores)
- Date (START TIME in `info`): 2026-10-31
