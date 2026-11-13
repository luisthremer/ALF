# Apple silicon (M1 Pro) — `GPU=apple` prototype

Status: **experimental prototype, not recommended for production** (no measurable sweep speedup, Green-function
precision ~30–40× worse than the CPU build). Full study: `plan/APPLE_GPU_FEASIBILITY.md` next to the clone
(benchmarks in `plan/apple_gpu_bench/`).

## What it is

Apple GPUs have no float64 (`'double' is not supported in Metal`). `GPU=apple` is the plain CPU build
(`ALF_GPU` stays unset, `Prog/` compiles exactly as for `GPU=none`, BLAS/LAPACK = Accelerate) plus one Objective-C
object, `Libraries/Modules/alf_la_metal_shim.m`, that defines the Fortran symbol `zgemm_`:

* ZGEMMs with all of M, N, K ≥ `ALF_METAL_MIN` (default 512) are split by columns of C: the GPU computes the first
  `ALF_METAL_FRAC` (default 0.35) of the columns in float-float arithmetic (~48 bits, runtime-compiled Metal kernel,
  `MTLMathModeSafe`), Accelerate computes the rest concurrently. Host memory is wrapped without copies
  (`newBufferWithBytesNoCopy` on the enclosing pages).
* All other ZGEMMs go to Accelerate (`cblas_zgemm`) unchanged. ZGEMMs called inside LAPACK are not affected.
* Inputs whose max-norm lies outside [2^-40, 2^40] fall back to the CPU (checked on the GPU).
* `ALF_METAL_VERBOSE=1` prints call counts and times at exit; `ALF_METAL_FASTMATH=1` is a planted-defect control.

## Build and run

```sh
source configure.sh GNU noMPI NO-INTERACTIVE GPU=apple && make lib program
ALF_METAL_VERBOSE=1 ./ALF.out
```

## Results (2026-11-13, M1 Pro 16-core GPU, Hubbard defect parameters `docs/gpu/scripts/defect/params/`, 2 sweeps)

| case | CPU build (s) | GPU=apple (s) | notes |
|---|---|---|---|
| L=16 (N=256) | 12.48, 12.58 | 12.72, 12.59 | no ZGEMM reaches the threshold (0 split calls) |
| L=24 (N=576) | 116.19, 115.19, 110.58 | 115.97, 114.97, 113.98 | 718 split calls: 2.91 s → 2.28 s; difference within noise |

Fixed-seed comparison L=24 (`compare.py`, tol 1e-8): PASS, worst relative deviation 2.1e-12, identical field
configurations and acceptance. Controls: `ALF_METAL_FRAC=0` gives 0.0 deviation (bit-identical); the planted defect
`ALF_METAL_FASTMATH=1` FAILs (3.3e-7). Precision Green mean/max: CPU 8.4e-15 / 1.6e-12, GPU=apple 2.4e-13 / 6.7e-11
(phase 1.9e-12 → 8.2e-11): within ~1e-10 but **not** at CPU level.

Why no speedup: the split ZGEMM itself is 1.2–1.45× faster than Accelerate for square N ≥ 576, but large square
ZGEMMs are only ~2.5 % of an L=24 checkerboard sweep; the time goes to checkerboard hopping, QR with pivoting
(ZGEQP3) and the submatrix updates (rank-32 flushes, where the GPU is not faster). The float-float GPU kernel alone
reaches 290 GF/s vs 460–520 GF/s for Accelerate's AMX ZGEMM.
