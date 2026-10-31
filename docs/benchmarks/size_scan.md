# Size scan: standard vs delayed vs submatrix

Mz Hubbard model on an L x L square lattice, U=4, Beta=2, Dtau=0.1 (Ltrot=20), Ltau=0, LOBS_ST=LOBS_EN=1, NBin=2, default Nwrap=10, N_Global=1, same seed for all schemes. sec/sweep = wall time of the run (`CPU Time` in `info`) / (NSweep*NBin); NSweep differs per L (same for the three schemes). Speedup = standard / scheme. Produced by `testsuite/lowrank/bench.sh` and `testsuite/lowrank/bench_report.py`.

| L | N | N_delay | threads | NSweep | standard [s/sweep] | delayed [s/sweep] | submatrix [s/sweep] | delayed speedup | submatrix speedup | acceptance |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 | 36 | 32 | 1 | 3600 | 0.007239 | 0.00661 | 0.007006 | 1.10 | 1.03 | 0.431789 |
| 8 | 64 | 32 | 1 | 1000 | 0.02362 | 0.02003 | 0.02098 | 1.18 | 1.13 | 0.431750 |
| 10 | 100 | 32 | 1 | 350 | 0.07246 | 0.05442 | 0.0546 | 1.33 | 1.33 | 0.432338 |
| 12 | 144 | 32 | 1 | 140 | 0.1815 | 0.1227 | 0.1185 | 1.48 | 1.53 | 0.432044 |
| 14 | 196 | 32 | 1 | 66 | 0.368 | 0.219 | 0.2103 | 1.68 | 1.75 | 0.432286 |
| 16 | 256 | 32 | 1 | 35 | 1.098 | 0.7277 | 0.6997 | 1.51 | 1.57 | 0.431994 |
| 20 | 400 | 32 | 1 | 12 | 2.56 | 1.178 | 1.121 | 2.17 | 2.28 | 0.430935 |

## Provenance

- Git commit: 26b354c7
- Compiler: GNU Fortran (Homebrew GCC 16.2.0) 16.2.0
- Flags: `-std=f2008 -cpp -O3 -ffree-line-length-none -ffast-math -fopenmp`
- BLAS/LAPACK: Apple Accelerate (vecLib), linked as `-llapack -lblas`
- Threads: OMP_NUM_THREADS = VECLIB_MAXIMUM_THREADS = 1
- Machine: Apple M1 Pro (8 performance + 2 efficiency cores)
- Date: 2026-10-31

## Notes

- Acceptance is identical between the three schemes at every L (same seed, identical Metropolis trajectories).
- The standard time at L=16 (N=256) lies above the N^2.4 trend of its neighbours (1.098 s/sweep; a re-run gave 1.093),
  so it is reproducible and not a timing outlier of one run; it lowers the speedups at L=16 relative to the trend.
- Single runs, no repetitions, so timings carry a few-percent noise. Total scan wall time: 1017 s.
