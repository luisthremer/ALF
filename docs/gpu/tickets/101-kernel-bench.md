# 101 – Kernel benchmark harness (HOST + CUDA backends)

- **Phase:** 1
- **Owner agent:** fortran-gpu-dev
- **Depends on:** –
- **Status:** open

## Goal
Write a standalone Fortran benchmark, outside ALF's build, that times the operations the GPU port
depends on. It runs on device-resident complex(8) N×N matrices with nvfortran OpenMP target offload
plus cuBLAS/cuSOLVER, and is checked against a HOST backend. Its numbers feed ADR-001/002
(programming model, `alf_la` design) and BENCHMARKS.md.

## Context
- Read `docs/gpu/HARD.md` (items 1, 2, 3, 5) and `docs/gpu/PLAN.md` §1 first. They explain why each
  kernel matters.
- Smoke tests already done on nividia (`docs/gpu/MACHINES.md`):
  - nvfortran 26.9 `-mp=gpu -gpu=cc80` works;
  - `use cublas` + `cublasZgemm` on OpenMP-mapped arrays via `!$omp target data use_device_addr`
    works;
  - calling the **host** `zgemm` name on device addresses fails, so device library calls must go
    through explicit backend routines.
- Environment: `. ~/alf-gpu/env/nvhpc.sh` (NVHPC 26.9, CUDA, cuBLAS, cuSOLVER; nvfortran modules
  `cublas`, `cusolverDn`; host BLAS/LAPACK from NVHPC `-lblas -llapack`). GPU: one A30 (sm_80, 24 GB).
- The design should preview the planned `alf_la` layer, because ticket 103 adds an ifx/oneMKL
  backend to this harness:
  - a module `bench_la` with backend-neutral routines (`la_gemm`, `la_gemm_batched`, `la_geqrf`,
    `la_ungqr`, `la_ormqr`, `la_ger`, …) taking device-resident arrays;
  - implementations selected by `-DBENCH_LA_CUDA` / `-DBENCH_LA_HOST`;
  - kernels written as OpenMP `target teams distribute parallel do` on plain arrays.
  The HOST build must also work with `-mp` (no GPU) or plain host code.
- Problem shapes (Hubbard, square lattice L×L, N = L², F = 2 flavours, all complex(8)):
  - **B_T, checkerboard:** 4 bond families (x-even, x-odd, y-even, y-odd). Bonds within a family
    are disjoint. Each bond (i,j) applies a real 2×2 [[c, s],[s, c]] (c = cosh(dτ t), s = sinh(dτ t),
    dτ t = 0.1) to rows i, j of G (left multiply), and the inverse [[c, −s],[−s, c]] to columns i, j
    (right multiply). Families are applied in a fixed order; the order must match the HOST reference.
  - **B_T, dense:** real symmetric E = exp(−dτ T) for the square-lattice hopping T (compute E on the
    host once, e.g. by diagonalising T with LAPACK `dsyev`). Operation: G ← E G E⁻¹. Variants:
    1. ZGEMM with a complex copy of E;
    2. DGEMM on the interleaved real view, which works directly for right-multiplication (complex
       N×N as real 2N×N); for left-multiplication find the cheapest correct layout trick or explain
       why ZGEMM is needed.
  - **Op_V diagonal wrap:** for one site i: row i of G × e^{α}, column i × e^{−α}. Variants:
    (a) one kernel launch per site, sequentially for all N sites (as the MVP would do, interleaved
    with updates); (b) all N sites of a slice in one fused kernel (G_kl ← e^{α_k} G_kl e^{−α_l}).
  - **Stabilization:**
    - Input: a graded matrix A = Q₁·diag(d)·R₁ with Q₁ unitary, R₁ random upper-triangular
      well-conditioned, and d log-spaced from 1 to 1e-12 (also 1e+12…1e-12).
    - Device path: column norms → stable sort descending → permute columns → `geqrf` → `ungqr` (form
      Q) → also time `ormqr` (apply Qᴴ to an N×N matrix).
    - Host path A: single-thread LAPACK `ZGEQP3` + `ZUNGQR` (what ALF does today).
    - Host path B: the same pre-pivot + `ZGEQRF` + `ZUNGQR` as the device.
  - **Updates:** G ← G + x yᵀ.
    - Rank-1 loop: n = 256 sequential rank-1 updates, one kernel or cuBLAS `zgeru` per update.
      Variant "with decision": before each update, copy one complex scalar from device to host and
      back (e.g. G(i,i)) to emulate the host accept/reject round trip.
    - Rank-k: accumulate k updates as X (N×k), Y (k×N) and apply one ZGEMM, for
      k = 1, 2, 4, …, 128 (time per update = total / k).
  - **Flavour batching:** two independent N×N ZGEMMs as two calls vs one
    `cublasZgemmStridedBatched` (or `cublasZgemmBatched`).
- Sizes N = 64, 144, 256, 400, 576, 1024 (L = 8, 12, 16, 20, 24, 32).
- Timing:
  - `omp_get_wtime` around `≥ 5` repetitions after 1 warm-up, with device synchronisation
    (`!$omp taskwait` / `cudaDeviceSynchronize` / a target-region barrier) before reading the clock;
  - report the median;
  - HOST timings single-threaded (`OMP_NUM_THREADS=1`, NVHPC BLAS sequential or
    `OPENBLAS_NUM_THREADS=1`); 16 cores = 16 independent ALF chains is the CPU bar (HARD.md).

## Files allowed
- `gpu_bench/**` (new directory at the repo root: sources, `Makefile`, `README.md`, `run.sh`,
  `results/` with plain-text result tables)

Nothing in `Prog/`, `Libraries/`, `configure.sh` or the testsuite. Anything else: stop and report.

## nividia rules for this ticket
- Work dir `~/alf-gpu/work/t101/` (rsync `gpu_bench/` there from the laptop repo). Don't touch
  other `~/alf-gpu/work/*` dirs: other tickets run concurrently.
- **Budget: the GPU + 2 CPU cores.** Other agents use the other 14 cores, so host timings run one
  at a time, single-threaded.
- Each run ≤ 15 min. No processes left behind; free the GPU when done (`nvidia-smi` shows no
  process of yours). Install nothing.

## Acceptance
- `make BACKEND=cuda` and `make BACKEND=host` build cleanly with nvfortran 26.9 on nividia.
- Correctness (printed by the program, PASS/FAIL per kernel and N): device vs HOST results with
  relative Frobenius difference
  - ≤ 1e-12 for propagation, wrap and updates;
  - for QR: (i) ‖A·P − Q·R‖/‖A‖ ≤ 1e-12 on the device; (ii) |R_ii| from device vs host path B agree
    to ≤ 1e-10 relative; (iii) report (don't fail on) the ratio of host ZGEQP3 |R_ii| to pre-pivoted
    |R_ii| (min/max over i), because that is the numerical difference ticket 104 studies.
- Results: one table per kernel, N × variant, median time in µs, plus GFLOP/s or GB/s where
  meaningful, written to `gpu_bench/results/nividia_A30_<date>.txt`.
- `gpu_bench/README.md`: how to build and run, what each kernel emulates in ALF (with ALF file:line
  from HARD.md).
- The report states: total runtime of the full run, the per-update latency of rank-1 with and
  without the host round trip, the rank-k k at which the time per update stops improving, and the
  crossover N where device beats single-core host for each kernel.

## Out of scope
- OpenACC (ticket 102 decides later), the ifx/oneMKL backend (103), integration into ALF, MAGMA,
  mixed precision.
