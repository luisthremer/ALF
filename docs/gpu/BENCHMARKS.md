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

## Phase 1: kernel benchmark (ticket 101, nividia A30, 2026-10-27)

Standalone `gpu_bench/` (nvfortran 26.9, OpenMP target + cuBLAS/cuSOLVER) vs one host core (NVHPC BLAS/LAPACK).
Medians in µs; full tables in `gpu_bench/results/nividia_A30_2026-10-27.txt` (reproduced by the manager, within
noise). Host load ≈10 during the runs (other tickets), so host numbers are slightly pessimistic.

| Kernel (per call) | N=144 device / host | N=400 device / host | N=1024 device / host |
|---|---|---|---|
| Checkerboard B_T (G ← B G B⁻¹, 8 family kernels) | 96 / 191 | 139 / 1485 | 392 / 9983 |
| Dense B_T, 2 ZGEMM | 53 / 1095 | 172 / 22941 | 2757 / 388561 |
| Dense B_T, DGEMM on real view + transposes | 59 / – | 168 / – | 1890 / – |
| Op_V wrap, one launch per site × N | 1621 / 199 | 4211 / 1087 | 10924 / 10002 |
| Op_V wrap, all sites fused | 13 / 36 | 17 / 276 | 40 / 1855 |
| Stabilization QR (pre-pivot + geqrf + ungqr) vs host ZGEQP3 + ZUNGQR | 1058 / 1834 | 2822 / 25721 | 9956 / 373241 |
| Rank-1 update (cuBLAS zgeru) | 4.7 / 6.8 | 6.6 / 71 | 31.6 / 353 |
| Rank-1 + host round trip of one scalar | 28 / – | 29 / – | 53 / – |
| Rank-k, per update, k=128 | 0.21 / – | 0.37 / – | 1.25 / 189 |

What this means:
- **Stabilization, dense B_T and fused element-wise kernels are fast on the device** (5–40× one core at N ≥ 400).
- **Per-site work is launch- and latency-bound:** one OpenMP target region ≈ 12–15 µs (7–9 µs with `nowait`),
  cuBLAS call ≈ 5 µs, scalar D→H→D round trip 16–23 µs. A rank-1 sweep with host decisions costs
  ≈ 2·L·N × 28 µs: 0.32 s/sweep at N=144 (L=40), versus 0.038 s per chain-sweep for the 16-core node,
  i.e. ~8× slower. At N=1024, ≈ 4.4 s vs ≈ 12 s per chain-sweep for the node.
- **Rank-k (k=128) is 25–40× cheaper per update than rank-1 on the device:** this quantifies the value of
  delayed or submatrix updates.
- Strided-batched ZGEMM over 2 flavours helps only at N ≤ 256.
- cuBLAS needs 16-byte-aligned host alpha/beta (cuDoubleComplex); Fortran complex(8) constants are 8-byte
  aligned, which caused crashes (see `gpu_bench/README.md`).

## Phase 1: launch overhead, OpenMP target vs OpenACC (ticket 102, nividia A30, nvfortran 26.9)

Medians in µs; full tables in `gpu_bench/results/nividia_A30_launch_2026-10-27.txt` (re-run by the manager, within noise).
Per site = per proposed flip, ≈ 50 % accepted; async rows and device decision: one wait per 256 sites.

| Variant | N=144 OpenMP / OpenACC | N=400 | N=1024 |
|---|---|---|---|
| Empty kernel, sync | 12.4 / 8.4 | 12.2 / 8.5 | 12.1 / 8.4 |
| Empty kernel, async (`nowait` / `async`) | 11.0 / 2.7 | 9.9 / 2.6 | 9.9 / 2.6 |
| Rank-1 kernel, async | 6.1 / 3.6 | 8.6 / 6.9 | 34 / 31 |
| cuBLAS zgeru | 3.9 / 3.9 | 6.5 / 6.7 | 31 / 31 |
| Host decision, per site (G(i,i) D→H + conditional launch) | 15.2 / 15.4 | 18.0 / 17.9 | 34 / 33 |
| Device decision, 2 async kernels, per site | 16.9 / 6.5 | 17.8 / 8.8 | 38 / 28 |

- The model difference is the asynchronous launch: OpenMP `nowait` costs ≈ 10 µs (task creation), OpenACC `async` 2.6 µs.
- Projected single-site update part only (2·L·N sites, L = 40) at N = 144: host decision 0.17 s/sweep, OpenACC device
  decision 0.075 s/sweep; the 16-core node does the **whole** sweep in 0.038 s per chain. At N = 1024 everything is
  bandwidth-bound (≈ 2.3–3.1 s/sweep vs 12 s per chain-sweep for the node).
- Decision (ADR-001): OpenMP only; the per-site lever is delayed/submatrix updates.
