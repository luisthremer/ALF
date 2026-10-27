# Benchmarks

## Phase 0: CPU baseline (nividia, upstream code, no GPU)

Machine: 16 Zen 3 cores of an EPYC 7763 VM (MACHINES.md). This is **not** a full HPC node, so
GPU-vs-CPU comparisons later must say so. Hubbard square lattice, U=4, β=4, Δτ=0.1 (L=40),
Nwrap=10, Mz, Ltau=1. GNU 13.3 + OpenBLAS, default stabilization, 1 OpenMP/BLAS thread per rank.

| N (lattice) | hopping | s/sweep per rank, 16 ranks (node full) | s/sweep, 1 rank alone* | node throughput, chain-sweeps/s |
|---|---|---:|---:|---:|
| 36 (6×6) | cb | 0.029 | – | 550 |
| 36 (6×6) | dense | 0.034 | – | 470 |
| 144 (12×12) | cb | 0.61 | 0.60 | 26 |
| 144 (12×12) | dense | 0.86 | 0.80 | 19 |
| 400 (20×20) | cb | 9.1 | 8.0 | 1.8 |
| 400 (20×20) | dense | 14.8 | 13.7 | 1.1 |
| 1024 (32×32) | cb | – | 193 | (≈0.08) |
| 1024 (32×32) | dense | – | 294 | (≈0.05) |

\* From the profiling runs (PROFILE.md): includes setup and gdb sampling overhead, and other jobs were
running on the node at the same time, so these are indicative only. Node throughput is 16 / (s/sweep
per rank); the 32×32 values in brackets assume 16 × the single-rank rate.

- ifx 2025.3 + MKL: 0.78 s/sweep at 12×12 cb, i.e. 1.3× slower than GNU on this AMD host.
- STABLOG costs the same as the default (0.616 vs 0.608 s/sweep).
- Per-kernel breakdown: PROFILE.md.

**What a GPU must beat:** one GPU chain needs about 16 × (single-rank s/sweep) less time per sweep to
match this node, unless the GPU runs several chains at once (HARD.md, "The bar a GPU has to clear").
