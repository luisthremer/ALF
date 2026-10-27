# gpu_bench: kernel benchmark harness (ticket 101)

Standalone Fortran micro-benchmarks for the operations the ALF GPU port depends on.
Outside ALF's build; nvfortran 26.9 (`cuda`, `host` backends) and ifx 2025.3 + oneMKL (`onemkl` backend, ticket 103).
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

## ifx / oneMKL backend (ticket 103)

```sh
ssh nividia; . ~/alf-gpu/env/intel.sh; cd ~/alf-gpu/work/t103/gpu_bench       # rsync gpu_bench/ here first
make BACKEND=onemkl                 # BLAS + QR family through "!$omp dispatch"      -> build/onemkl_dispatch/bench
make BACKEND=onemkl LAPACK=shim     # BLAS via dispatch, geqrf/ungqr/unmqr via the C++ SYCL shim -> build/onemkl_shim/bench
make BACKEND=onemkl shimtest        # shim alone vs host LAPACK incl. getrf/getri   -> build/shim_test
./run_onemkl.sh build/onemkl_dispatch/bench [N-list]    # host fallback (no Intel GPU): sets the runtime environment
./run_onemkl.sh build/shim_test
```

Exact flags (all verified to build and run on nividia):

| Step | Command |
|---|---|
| Fortran compile | `ifx -O2 -g -traceback -fiopenmp -fopenmp-targets=spir64 -qmkl -I$MKLROOT/include -DBENCH_LA_ONEMKL [-DBENCH_LA_SHIM] -module <dir> -c src/X.F90` (JIT `spir64`; no `ocloc` on nividia, so no `spir64_gen` AOT) |
| shim compile | `icpx -fsycl -O2 -g -std=c++17 -I$MKLROOT/include -c src/bench_sycl_shim.cpp` (no device kernels, so the object is a plain host object and ifx can link it) |
| link, dispatch build | `ifx <same flags> *.o -o bench` (`-qmkl` together with `-fopenmp-targets` adds `-lmkl_sycl_*`, `libsycl`, `libomptarget` by itself, checked with `ldd`; no `-lOpenCL`/`-lstdc++` needed) |
| link, shim build | `ifx <same flags> *.o bench_sycl_shim.o -o bench -lmkl_sycl_lapack -lsycl -lstdc++` |

`include "mkl_blas_omp_offload.f90"` / `"mkl_lapack_omp_offload.f90"` at file scope in `bench_la.F90`; the interface modules
`onemkl_blas_omp_offload_lp64` / `onemkl_lapack_omp_offload_lp64` are `use`d inside each `la_*` routine
(`only:` the one routine), because ifx rejects a module that re-exports `bind(C)` names which another module declares `external`.

Source changes needed for ifx (also built with nvfortran, unchanged behaviour): `bench_sec.F90` no longer re-declares
`external :: zgemm, zgeru` (they arrive from `bench_ref` by use association; ifx error #6406); `rand_c` uses allocatable
instead of automatic arrays (2 N^2 reals overflow ifx's 8 MB stack: SIGSEGV at N = 576).

### Runtime environment on nividia (`run_onemkl.sh`)

* **ifx 2025.3 code must run on the 2025.3 OpenMP runtime.** The env file puts `compiler/2026.2/lib` first
  (oneAPI `latest`: MKL 2026.2, icpx 2026.2). With it, `!$omp interop` fails at run time with
  `undefined symbol: __tgt_create_interop` (the 2026.2 `libomptarget` exports `__tgt_interop_init` instead).
  `run_onemkl.sh` prepends `compiler/2025.3/lib`; `libsycl.so.9` (needed by `libmkl_sycl_*`) exists only in 2026.2 and is
  still found. On SuperMUC with one consistent oneAPI release this is not needed. The dispatch build without
  interop runs with either runtime.
* **Do not set `ONEAPI_DEVICE_SELECTOR=opencl:cpu` for the whole process.** The OpenMP runtime honours it, treats the
  OpenCL CPU device as an offload device and aborts (`getExtensionsInfoForDevice failed (Extension DeviceAttributeQuery
  is required ...)`). The shim therefore selects the SYCL CPU device itself (`sycl::cpu_selector_v`, OpenCL CPU runtime of
  oneAPI works on nividia); `ALFSH_FALLBACK=default` switches to the default selector.
* nividia has no Intel GPU: `omp_get_num_devices() = 0`, so every `!$omp target` region and every `!$omp dispatch` runs
  on the host (the variant is not selected, the base MKL routine is called). `OMP_TARGET_OFFLOAD=DISABLED` behaves the
  same; `MANDATORY` aborts (as it should). **Timings in host fallback carry no meaning for PVC**; they only show that the
  code path runs (the device-resident arrays are the host arrays).

### Where each routine runs on PVC (what the backend calls, and evidence)

| Routine (bench use) | OpenMP dispatch variant exists? | Class on PVC | Evidence / caveat |
|---|---|---|---|
| ZGEMM, DGEMM, ZGEMM_BATCH_STRIDED, ZGERU (`la_gemm`, `la_gemm_cr`, `la_gemm_batched`, `la_geru`) | yes | **dispatch-device** | `$MKLROOT/include/mkl_blas_omp_offload_lp64.f90`: `!$omp declare variant(... mkl_blas_<x>_omp_offload_lp64) match(construct={dispatch}, device={arch(gen)}) ... adjust_args(need_device_addr:a,b,c)`; the symbols are referenced by the built binary (`nm -D`). That they are fully on the device is the documented BLAS behaviour (R3), not testable here. |
| ZGEQRF (`la_geqrf`) | yes (`mkl_lapack_zgeqrf_omp_offload_lp64` in `libmkl_sycl_lapack.so`) | **dispatch-hybrid** (R3) | Header + exported symbol confirm the variant exists. **No installed header or doc says which part runs on the host**; "hybrid, synchronous" is R3's claim and stays unverified until a PVC run. |
| ZUNMQR (`la_ormqr`) | yes | **dispatch-hybrid** (R3) | same |
| ZGETRF (not in the bench sections; shim_test only) | yes | **dispatch-hybrid** (R3) | same |
| ZUNGQR (`la_ungqr`) | yes (interface present) | **host** (R3: "interface only, computed on the CPU") | Variant symbol exists, so device pointers are accepted, but R3 says the work is done on the CPU (with transfers). Unverified. |
| ZGETRI | yes | unknown (R3 does not list it) | `mkl_lapack_zgetri_omp_offload_lp64` exists; behaviour unknown |
| ZGEQP3 | **no** (no `declare variant`, no exported symbol, no `geqp3` in `oneapi/mkl/lapack/lapack.hpp`) | **host** only | and no SYCL API either: no pivoted QR on the device on either path, hence the pre-pivot design (HARD.md item 2) |
| geqrf/ungqr/unmqr/getrf/getri via **shim** (`oneapi::mkl::lapack::*` USM API, `bench_sycl_shim.cpp`) | n/a | **shim** (device; USM API, R3/HARD.md: "USM device versions", not hybrid-labelled) | headers `oneapi/mkl/lapack/lapack.hpp` have the USM overloads and `*_scratchpad_size`; verified on the SYCL **CPU** device only (`shim_test`: all five routines, N = 64, 144, 400). On PVC unverified. |

### How the SYCL queue is obtained

`la_init` (shim build): `!$omp interop init(prefer_type("sycl"), targetsync: iop)`; if `iop /= omp_interop_none`, check
`omp_get_interop_int(iop, omp_ipr_fr_id) == omp_ifr_sycl (4)`, then `qh = omp_get_interop_ptr(iop, omp_ipr_targetsync)`, which
the shim casts to `sycl::queue*` (`omp_ipr_targetsync = -9`, `omp_ipr_fr_id = -1`, values from `compiler/2026.2/.../omp.h`).
ifx 2025.3's `omp_lib` has the interop kinds and `omp_ifr_*` constants but **not** `omp_get_interop_ptr/int`, so they are
declared locally with `bind(C)`. Device pointers for the shim come from `!$omp target data use_device_addr(...)`.
What was verified on nividia: the directive compiles with ifx 2025.3, links (`__tgt_create_interop`, `omp_get_interop_ptr`),
and with no offload device returns `omp_interop_none`; the shim then uses its private CPU-device fallback queue (host
pointers). **What could not be verified without a GPU:** that the interop object on PVC really carries a `sycl::queue*`
(and not a Level Zero handle) for `prefer_type("sycl")`, that it is the in-order queue (`omp_ipr_is_inorder = 11` can be queried), and that
USM pointers from `use_device_addr` are valid in that queue's context. Fallback if it fails: `prefer_type("level_zero")` +
`sycl::make_queue<backend::ext_oneapi_level_zero>` from the returned L0 handles (not tried).

### Things only SuperMUC (PVC) can verify

1. Real device execution: `OMP_TARGET_OFFLOAD=MANDATORY` run of the dispatch and shim builds; `ZE_AFFINITY_MASK`/stack choice.
2. The interop-to-`sycl::queue*` path above (`alfsh_describe` prints `queue from OpenMP interop, device ... (gpu)`).
3. Which of zgeqrf/zunmqr/zungqr/zgetrf really run on the host under dispatch (time `dev geqrf/ungqr/ormqr` for both builds and watch
   PCIe traffic or `unitrace`), and the shim's speed relative to them.
4. Whether `lwork = -1` queries (done on the host here) give a workspace sufficient for the dispatch variant on the device.
5. Accuracy checks with real device arithmetic (all correctness tolerances unchanged; host fallback only proves the plumbing).
6. All performance numbers (launch latency of `!$omp target` / `dispatch` / shim on PVC, rank-1 round trip, GEMM rates).
7. AOT (`-fopenmp-targets=spir64_gen -Xopenmp-target-backend "-device pvc"`) build and JIT compile time; `-qmkl=sequential` vs parallel
   choice; link flags with the SuperMUC module set.

## Layout

| File | Role |
|---|---|
| `src/bench_la.F90` | **Backend layer** (preview of `alf_la`): `la_gemm`, `la_gemm_cr` (complex x real via DGEMM on the interleaved real view), `la_gemm_batched`, `la_geru`, `la_geqrf/ungqr/ormqr` (+ `_lwork`), `la_sync`. `-DBENCH_LA_CUDA`: cuBLAS/cuSOLVER through `bind(C)` interfaces with device pointers from `!$omp target data use_device_addr`. `-DBENCH_LA_HOST`: host BLAS/LAPACK. `-DBENCH_LA_ONEMKL` (ticket 103): oneMKL `!$omp dispatch` interfaces, with `-DBENCH_LA_SHIM` the QR family goes through `bench_sycl_shim.cpp`. Takes device-resident arrays. Nothing else knows the backend. |
| `src/bench_sycl_shim.cpp` | C++ SYCL/oneMKL shim (`extern "C"`, queue as `void*`, USM device pointers, scratchpad owned by the caller, errors returned as status): `alfsh_z{geqrf,ungqr,unmqr,getrf,getri}` + `_ws` scratchpad-size queries. |
| `src/shim_test.F90` | Standalone check of the shim against host LAPACK. |
| `run_onemkl.sh` | Runtime environment for the ifx binaries on nividia (see above). |
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

## Launch overhead: OpenMP target vs OpenACC (ticket 102)

`src/bench_launch.F90` (standalone; reuses `bench_util`/`bench_ref`/`bench_la` unchanged). One source, model chosen by macro:

```sh
. ~/alf-gpu/env/nvhpc.sh; cd ~/alf-gpu/work/t102/gpu_bench
make launch MODEL=omp   # nvfortran -mp=gpu  -DMODEL_OMP   -> build/launch_omp/bench_launch
make launch MODEL=acc   # nvfortran -acc=gpu -DMODEL_ACC   -> build/launch_acc/bench_launch
OMP_TARGET_OFFLOAD=MANDATORY ./build/launch_omp/bench_launch [N-list]   # default 64,144,256,400,1024
./build/launch_acc/bench_launch
```
Flags (checked with `nvfortran -help -acc`/`-gpu`): `-O2 -g -traceback -gpu=cc80 -cudalib=cublas,cusolver`; the support
modules are always built with `-mp=gpu`, the final link adds `-acc=gpu` for MODEL=acc. Data are mapped once with the model's
own `enter data`; the cuBLAS reference gets device pointers from `host_data use_device` / `use_device_addr`.
Variants: empty kernel (sync; async = 1000 launches, one wait: `nowait depend(inout)`+`taskwait` / `async(1)`+`wait(1)`), rank-1
(sync, async, cuBLAS zgeru; 256 updates, one wait), single-site path (host decision = D2H of G(i,i) + launch if accepted;
device decision with 2 async kernels: 1-thread decide + conditional update; device decision with ONE gang/team: decide,
barrier, conditional update), Op_V per-site launches (sync, async). Decision: r = |1+G(i,i)|, accept iff u < r/(1+r), u
pre-generated (accept rate ~ 50%). All checks: G vs host `zgeru` sequence <= 1e-12 and the accept pattern of host- and
device-decision variants identical to the host reference (0 mismatches). Results: `results/nividia_A30_launch_<date>.txt`.
Race note: a one-kernel multi-gang decision+update is not possible (the updating gangs overwrite G(i,i) that the others
read), hence the 2-kernel and single-gang forms.
