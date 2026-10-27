# gpu_bench: kernel benchmark harness (ticket 101)

Standalone Fortran micro-benchmarks for the operations the ALF GPU port depends on.
Outside ALF's build; nvfortran 26.9 only (ticket 103 adds an ifx/oneMKL backend to the same layer).
All data are complex(8) N x N, sizes N = 64, 144, 256, 400, 576, 1024 (L = 8..32).

## Build and run (nividia)

```sh
ssh nividia; . ~/alf-gpu/env/nvhpc.sh; cd ~/alf-gpu/work/t101   # rsync gpu_bench/ here first
make BACKEND=cuda         # OpenMP target offload + cuBLAS/cuSOLVER   (build/cuda/bench)
make BACKEND=host         # same kernels with nvfortran -mp on the host + host BLAS/LAPACK
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
OMP_TARGET_OFFLOAD=MANDATORY ./build/cuda/bench [64,144,...]   # fails loudly if it would run on the host
./build/host/bench [64,144,...]
BENCH_SECTIONS=qr,upd ./build/cuda/bench 256                    # subset: cb,dense,opv,qr,upd,batch
./run.sh                                                         # both backends -> results/nividia_A30_<date>.txt
```
Exit status 1 if any check FAILs (all checks are printed as `PASS/FAIL name N value tol`).
Host timings are single-threaded (the CPU bar in HARD.md is 16 independent single-core chains).

## Layout

| File | Role |
|---|---|
| `src/bench_la.F90` | **Backend layer** (preview of `alf_la`): `la_gemm`, `la_gemm_cr` (complex x real via DGEMM on the interleaved real view), `la_gemm_batched`, `la_geru`, `la_geqrf/ungqr/ormqr` (+ `_lwork`), `la_sync`. `-DBENCH_LA_CUDA`: cuBLAS/cuSOLVER through `bind(C)` interfaces with device pointers from `!$omp target data use_device_addr`. `-DBENCH_LA_HOST`: host BLAS/LAPACK. Takes device-resident arrays. Nothing else knows the backend. |
| `src/bench_kern.F90` | Device kernels, `!$omp target teams distribute parallel do` on plain explicit-shape arrays (checkerboard family, Op_V wrap, rank-1, transpose, column norms, rank-count sort, column permutation). No library calls. |
| `src/bench_ref.F90` | HOST reference (plain Fortran + host LAPACK, never offloaded): checkerboard, Op_V, `dsyev` hopping exponential, ZGEQP3 / pre-pivot+ZGEQRF paths, graded test matrix. |
| `src/bench_sec.F90` | The benchmark sections (setup, device-vs-host check, timing). |
| `src/bench_util.F90`, `bench_main.F90` | Timing (`omp_get_wtime`, 1 warm-up + 21 reps, 11 if >50 ms, 5 if >0.5 s, median, `cudaDeviceSynchronize` before each clock read), checks, tables. |

Rules the kernels follow: arrays are mapped once (`target enter data`), never per kernel call; the only
host<->device transfers inside timed regions are the deliberate `+roundtrip` variants; all other
transfers (setup, verification) are marked in the source.
To add a backend (103): implement the `la_*` routines under a new `#ifdef` in `bench_la.F90`; the
kernels and sections stay untouched.

## What each kernel emulates in ALF

| Section | Emulates | ALF location |
|---|---|---|
| `cb` checkerboard B_T | G -> B G B^-1 with 4 bond families (8 launches: 4 left, 4 right inverse), one kernel per family, all N/2 disjoint bonds in parallel | `Prog/OpTTypes_mod.F90:145-199` calling `ZDSLSYMM` (`Libraries/Modules/Mat_subroutines_mod.F90:967`); HARD.md item 5 |
| `dense` B_T | G -> E G E^-1 with real symmetric E = exp(-dtau T): (1) 2 ZGEMM with complex copies; (2) DGEMM on the real view for the right product, left product as transpose + DGEMM + transpose (E symmetric: E H = (H^T E)^T); plus single-op variants | `ZLARCM` / `ZLACRM` calls in `Mat_subroutines_mod.F90:1236, 1401`; HARD.md item 5 |
| `opv` | `Op_Wrapup` for a diagonal Op_V on every site: (a) one launch per site (what interleaving with updates would do), (b) one fused kernel | `Prog/Operator_mod.F90:743`, loop in `Prog/Wrapgr_mod.F90:111-120`; HARD.md item 3 |
| `qr` | Stabilization on a graded matrix A = Q diag(d) R (d from 1 to 1e-12, or 1e12 to 1e-12; mat3 = columns of mat1 scaled randomly; mat4 = mild, d 1 to 1e-3). Device: column norms -> rank-count stable descending sort -> column permutation -> `geqrf` -> `ungqr`; `ormqr` applies Q^H. Host A: ZGEQP3 + ZUNGQR (what ALF does); host B: same pre-pivot + ZGEQRF + ZUNGQR. Timed operations include an N x N device copy (in-place factorization needs a fresh input). | `Prog/QDRP_decompose_mod.F90:78-82`, `Prog/udv_state_mod.F90:503-542`, `Prog/cgr1_mod.F90:305-312`, `Prog/cgr2_2_mod.F90:346,405,420`; HARD.md item 2 |
| `upd` | Rank-1 G += x y^T, 256 sequential updates: kernel, kernel with `nowait`+`depend` (asynchronous launches, no decision in between), cuBLAS ZGERU, and `+roundtrip` (one complex scalar device->host->device per update = the accept/reject decision on the host). Rank-k: one ZGEMM with k = 1..128, time per update = total/k. | `Prog/upgrade_mod.F90:222` (accept test), `:269` (ZGERU); HARD.md item 1 |
| `batch` | Two flavours: two ZGEMM calls vs `cublasZgemmStridedBatched`, and one ZGEMM | flavour loops, e.g. `Prog/Wrapgr_mod.F90:110-117`; HARD.md item 3 |

Rates: GFLOP/s uses complex-equivalent flops (8 N^3 per N x N complex product, so the DGEMM variants are
comparable); GB/s counts the nominal traffic (read+write of G once per launch) and is "effective".

## Correctness (device vs HOST reference)

Propagation, wrap and updates: relative Frobenius difference <= 1e-12. QR: (i) |AP - QR|/|A| <= 1e-12;
(ii) |R_ii| device vs host path B <= 1e-10 relative; (iii) ratio of ZGEQP3 to pre-pivoted |R_ii| is only reported.
Extra checks: `ormqr` (Q^H A P = R), permutation equality, `(ii-b)` the diagonal difference normalised by the
column norm (the quantity Householder QR actually controls). See `results/` and the ticket report for the
(ii) discussion: it cannot hold for strongly graded input, even between two host LAPACK variants.

## Gotchas found while writing this

- **cuBLAS reads host alpha/beta as 16-byte-aligned cuDoubleComplex.** A Fortran `complex(8)` PARAMETER or
  local is only 8-byte aligned: random SIGSEGV inside libcublas/libcublasLt (layout dependent, looked like heap
  corruption). `bench_la` copies alpha/beta into a malloc'ed array first.
- A device->host->device `!$omp target update` loop whose result is never used on the host is optimised away by
  nvfortran (measured 0.15 us per round trip instead of ~9 us). The `+roundtrip` variants consume the value.
- OpenMP target launches cost ~12-15 us each (synchronous); cuBLAS ZGERU ~5.5 us; see `upd` table.
