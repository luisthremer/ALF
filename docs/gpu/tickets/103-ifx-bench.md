# 103 – ifx / oneMKL backend for the kernel benchmark (host fallback on nividia)

- **Phase:** 1
- **Owner agent:** fortran-gpu-dev
- **Depends on:** 101 (done, reviewed)
- **Status:** done (2026-10-27; host-fallback only, manager re-ran dispatch, shim and shim_test: 0 failures)

## Goal
Make the 101 benchmark compile with ifx for the Intel GPU target (`spir64` JIT) with a oneMKL
backend, and run it correctly in host-fallback mode on nividia. Prototype the C++ SYCL shim for the
LAPACK routines that oneMKL's Fortran `!$omp dispatch` would run on the host.

## Context
- Read `docs/gpu/HARD.md` items 2 and 9, `docs/gpu/MACHINES.md` (ifx smoke tests), `gpu_bench/README.md`.
- Toolchain on nividia: `. ~/alf-gpu/env/intel.sh` (ifx 2025.3.3 at `compiler/2025.3/bin`, MKL
  2026.2 incl. SYCL, icpx 2026.2). `ocloc` is missing, so **no** `spir64_gen` AOT. JIT `spir64` only.
- oneMKL OpenMP offload from Fortran: `include "mkl_omp_offload.f90"`, `use onemkl_blas_omp_offload_lp64`,
  data via `target data` / `enter data`, `!$omp dispatch` directly before **one** MKL call. Link:
  `-fiopenmp -fopenmp-targets=spir64 -qmkl -lmkl_sycl -lsycl -lOpenCL -lstdc++`
  (verify against the MKL 2026 docs installed under `/opt/intel/oneapi/mkl/latest/`; record what
  actually works).
- Per R3 §3 [verify against installed MKL docs/headers]:
  - `zgeqrf`, `zunmqr` and `zgetrf` are "hybrid" (part on the host) under OpenMP dispatch;
  - `zungqr` is computed on the host;
  - `zgeqp3` is not offloadable.
  The shim: small C++ (`icpx -fsycl`) functions with `extern "C"` / `bind(C)` interfaces calling
  `oneapi::mkl::lapack::geqrf/ungqr/unmqr/getrf/getri` (USM device pointers, scratchpad from
  `*_scratchpad_size`), taking the device pointer from `use_device_addr`. Obtaining the SYCL queue
  that matches the OpenMP device is the open question: try `omp_get_interop` /
  `omp_get_interop_ptr(…, omp_ipr_targetsync)` and document what works on ifx 2025.3.
- Without an Intel GPU, `OMP_TARGET_OFFLOAD=DISABLED` (or default fallback) runs everything on the
  host. The SYCL shim then needs a host/CPU SYCL device (`ONEAPI_DEVICE_SELECTOR=opencl:cpu` or
  similar). Record which selector works.

## Additional task from the 101 review: representative QR test matrices
101's criterion (ii) was ill-posed for its row-mixed graded matrices (see 101 "Resolution"). Add to the `qr`
section, for all backends:
- **mat5:** column-graded A = M·diag(d), M a random well-conditioned complex matrix (cond ≲ 10), d log-spaced
  1e+12…1e-12 in random column order. This is the structure of ALF's `decompose` (B·U·D).
- **mat6:** two-sided graded A = diag(d₁)·M·diag(d₂) (the structure of CGR's U_R†U_L + D_R(V_R V_L)D_L).
- Checks: (i) and (ii-b) as in 101. (ii) relative |R_ii| device vs host path B ≤ 1e-10 on mat5. On mat6, report
  the value and the same quantity for host path A (ZGEQP3) vs host path B, but do not gate on it: that is the
  pre-pivot question, which ticket 104 answers at the ALF level.
Also apply this to `BACKEND=cuda` and `BACKEND=host` (rerun on nividia, GPU allowed for this re-run only, ≤ 5 min).

## Files allowed
- `gpu_bench/**` (add `-DBENCH_LA_ONEMKL` backend, `bench_sycl_shim.cpp`, Makefile target
  `BACKEND=onemkl`, README section)

## nividia rules
Work dir `~/alf-gpu/work/t103/`. Budget 4 CPU cores, no GPU needed. Runs ≤ 15 min, no leftovers.

## Acceptance
- `make BACKEND=onemkl` builds with ifx 2025.3 for `spir64`.
- With host fallback on nividia, every correctness check from 101 passes (same tolerances).
- `make BACKEND=cuda` and `BACKEND=host` are unchanged and still pass.
- README documents:
  - the exact compile/link flags;
  - per routine, whether it would run on the device on PVC: dispatch-device / dispatch-hybrid /
    host / shim;
  - how the SYCL queue is obtained;
  - what can only be verified on SuperMUC (performance, real device execution) → a short list for
    the owner's later run.
- Timings in host fallback are reported but carry no meaning for PVC; say so.

## Out of scope
PVC performance, AOT builds, integration into ALF.
