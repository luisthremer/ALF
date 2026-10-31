# Size scan: standard vs delayed vs submatrix

Mz Hubbard model on an L x L square lattice, U=4, Beta=2, Dtau=0.1 (Ltrot=20), Ltau=0, LOBS_ST=LOBS_EN=1, NBin=2, default Nwrap=10, N_Global=1, same seed for all schemes. sec/sweep = `CPU Time` / (NSweep*NBin); NSweep differs per L (same for the three schemes). Speedup = standard / scheme. Timing = wall time of the sweep phase only (`CPU Time` in `info`; the timer starts after the initial Green function computation, setup is excluded). Produced by `testsuite/lowrank/bench.sh` and `testsuite/lowrank/bench_report.py`.

| L | N | N_delay | threads | NSweep | standard [s/sweep] | delayed [s/sweep] | submatrix [s/sweep] | delayed speedup | submatrix speedup | acceptance |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 | 36 | 32 | 1 | 3600 | 0.007239 | 0.00661 | 0.007006 | 1.10 | 1.03 | 0.431789 |
| 8 | 64 | 32 | 1 | 1000 | 0.02362 | 0.02003 | 0.02098 | 1.18 | 1.13 | 0.431750 |
| 10 | 100 | 32 | 1 | 350 | 0.07246 | 0.05442 | 0.0546 | 1.33 | 1.33 | 0.432338 |
| 12 | 144 | 32 | 1 | 140 | 0.1815 | 0.1227 | 0.1185 | 1.48 | 1.53 | 0.432044 |
| 14 | 196 | 32 | 1 | 66 | 0.368 | 0.219 | 0.2103 | 1.68 | 1.75 | 0.432286 |
| 16 | 256 | 32 | 1 | 35 | 1.098 | 0.7277 | 0.6997 | 1.51 | 1.57 | 0.431994 |
| 18 | 324 | 32 | 1 | 14 | 1.423 | 0.7139 | 0.6767 | 1.99 | 2.10 | 0.431256 |
| 20 | 400 | 32 | 1 | 12 | 2.56 | 1.178 | 1.121 | 2.17 | 2.28 | 0.430935 |

## Provenance

- Git commit (of the binary, from `info`): 26b354c7
- Compiler: GNU Fortran (Homebrew GCC 16.2.0) 16.2.0
- Flags: `-std=f2008 -cpp -O3 -ffree-line-length-none -ffast-math -fopenmp`
- BLAS/LAPACK: Apple Accelerate (vecLib), linked as `-llapack -lblas`
- Threads: OMP_NUM_THREADS = VECLIB_MAXIMUM_THREADS = 1
- Machine: Apple M1 Pro (8 performance + 2 efficiency cores)
- Date (START TIME in `info`): 2026-10-31
- The commit/date of standard L=14 and L=16 come from the `info` of the re-run of those two points (same binary, same day), not from the original runs (rows without precision fields are from the original scan).

## Notes

- Acceptance is identical between the three schemes at every L (same seed, identical Metropolis trajectories).
- N_delay=32 is the maximum, not always the effective delay: the capacity is capped at min(N_delay, vertices per slice[, Ndim for submatrix]) and G is flushed at the end of every slice, so the effective delay is min(N_delay, accepts per slice); at L=6 (36 vertices per slice, about 16 accepts) N_delay=32 is never reached.
- From L=14 to L=16 (N=196 to 256) the time per sweep grows by x3.0 (standard), x3.3 (delayed) and x3.3 (submatrix); the N^2.46 trend fitted to standard at L=8-14 would give x1.9. In relative terms the two low-rank schemes jump more than standard, which is why their speedup dips (delayed 1.68 -> 1.51, submatrix 1.75 -> 1.57 from L=14 to L=16).
- L=18 (N=324) tests the guess that the jump at N=256 is a power-of-two / cache effect (leading dimension 256). Local exponents d ln t / d ln N, standard / delayed / submatrix: L=14-16: 4.09 / 4.50 / 4.50; L=16-18: 1.10 / -0.08 / -0.14; L=18-20: 2.79 / 2.38 / 2.40. Speedups at L=14, 16, 18, 20: delayed 1.68, 1.51, 1.99, 2.17, submatrix 1.75, 1.57, 2.10, 2.28. Observation: the speedup of both low-rank schemes recovers at L=18, consistent with a localized slowdown at N=256 (an observation, not proof; one run per point).
- Repeats (same binary and parameters, run later in the session; raw_repeat.txt, whose second row is the standard L=12 row of the N_delay scan, raw_ndelay_scan.txt, copied there): standard at L=12: 3 runs, 0.1815, 0.18, 0.1842 s/sweep (spread max/min - 1 = 2.3 %). This is the only noise estimate (one point); differences of a few percent between delayed and submatrix are not significant.
- L=20 has only 24 sweeps in total (NSweep=12, NBin=2).
