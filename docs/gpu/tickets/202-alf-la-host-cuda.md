# 202 – `alf_la` module: HOST and CUDA backends with unit tests

- **Phase:** 2
- **Owner agent:** fortran-gpu-dev (review: numerics-reviewer)
- **Depends on:** ADR-002 (accepted). Runs in parallel with 201 (see "Build flags").
- **Status:** done (2026-10-28; two review rounds, re-run by the manager)

## Goal
`Libraries/Modules/alf_la_mod.F90` provides the plain-array linear-algebra layer of ADR-002 §1 with two backends:
`HOST` (the exact BLAS/LAPACK calls ALF makes today; default when no `ALF_LA_*` macro is set) and `CUDA`
(cuBLAS/cuSOLVER on device-resident arrays). Every routine has a unit test that runs with both backends.
No upstream routine calls `alf_la` yet; `ALF.out` with `GPU=none` is unchanged.

## Context
- ADR-002 (interface rules, routine list, alignment, streams, error handling). Prototype to start from:
  `gpu_bench/src/bench_la.F90` (101/103: cuBLAS/cuSOLVER `bind(C)` interfaces, `use_device_addr`, 16-byte aligned
  alpha/beta – a complex(8) parameter crashes libcublas, lwork queries). Copy and adapt; don't `use` gpu_bench.
- **Routine set = what the GPU path will call.** Derive it from the hot-path call sites and put a table
  "routine → ALF call sites (file:line)" in the module header comment and the report. Starting list (verify):
  - `la_init`, `la_finalize`, `la_sync`, `la_backend_name`
  - `la_gemm` (ZGEMM), `la_gemm_cr`/`la_gemm_rc` (complex × real via DGEMM on the interleaved view, 101),
    `la_gemm_batched` (strided), `la_geru` (ZGERU)
  - `la_geqrf`, `la_ungqr`, `la_ormqr` (ZUNMQR), each with an `_lwork` query
  - `la_getrf`, `la_getrs` (CGR, `Prog/cgr1_mod.F90:124,141,164`), `la_trsm` (ZTRSM, `cgr1_mod.F90:571,619`);
    `la_getri` only if a call site needs it
  - `la_qrpp`: column-norm pre-pivot + `geqrf`, returning a ZGEQP3-compatible `IPVT` and `tau` – the 104 semantics,
    reference = the `#else` (variant 1) branch of `QR_PREPIVOT` in `Prog/cgr1_mod.F90` (≈ line 312 ff.) and its
    twin in `Prog/cgr2_2_mod.F90`. Column norms on the device as a kernel; the sort may be on the host (N
    numbers D→H, permutation H→D – a documented transfer).
  - `la_geqp3_host`: STABLOG's CGR/CGR2_2 path (108): `target update from` the matrix, host ZGEQP3,
    `target update to` factors, `tau`, `IPVT` (documented transfer). Under HOST it is plain ZGEQP3.
  - Row/column permutation (ZLAPMR, `cgr1_mod.F90:574`) and other small element-wise steps are **kernels for
    203–206**, not `alf_la` routines – list them in the report as "not in alf_la".
- Interface rules (ADR-002): explicit-shape arrays + leading dimensions in LAPACK order; caller guarantees device
  residency under CUDA; `alf_la` gets device addresses itself (`!$omp target data use_device_addr`); one in-order
  CUDA stream created in `la_init` and set on both handles; `la_sync` is the only blocking call; `info` checked,
  failure → `Terminate_on_error` (`Libraries/Modules/runtime_error_mod.F90:81`). Device `info`/workspace:
  preallocated once in `la_init` or grown on demand, never per call.
- `ALF_LA_ONEMKL`: not implemented here (302). That branch must fail at compile time with a clear `#error`.

### Build flags (until 201 lands)
- HOST: unchanged profiles, e.g. `. ./configure.sh GNU devel noMPI NO-INTERACTIVE` and IntelLLVM, PGI.
- CUDA: `. ./configure.sh PGI noMPI NO-INTERACTIVE`, then export by hand exactly what 201 will set:
  `ALF_FLAGS_MODULES`/`ALF_FLAGS_PROG`: replace `-mp` by `-mp=gpu -gpu=cc80`, append `-DALF_GPU -DALF_LA_CUDA`;
  `ALF_LIB`: append `-mp=gpu -gpu=cc80 -cudalib=cublas,cusolver`. Put this in a helper script.

## Files allowed
- `Libraries/Modules/alf_la_mod.F90` (new), `Libraries/Modules/Makefile` (add to `SRCS` only)
- `testsuite/alf_la.tests/` (new: `CMakeLists.txt` + test programs), `testsuite/CMakeLists.txt` (one
  `add_subdirectory` line)
- `docs/gpu/scripts/nividia/` (helper scripts)
Anything else: stop and report.

## Acceptance
1. **Tests** (`testsuite/alf_la.tests/`, registered with ctest): one program per routine family. Each builds its
   inputs on the host (seeded, reproducible), maps them (`target enter data`) under CUDA, calls `alf_la`, copies back
   and compares with a direct host BLAS/LAPACK call in the same program. Sizes N ∈ {1, 7, 64, 144, 400} (plus a
   rectangular and an `lda > m` case per routine; k ∈ {1, 32} for batched). Each test **fails (nonzero exit)** on a
   violated check – printing alone is not enough.
   - HOST backend: results **bit-identical** to the direct call.
   - CUDA backend: max relative difference (Frobenius, per output) ≤ 1e-12 for gemm/geru/trsm/getrs; for the QR
     family compare invariants, not raw factors: ‖A·P − Q·R‖/‖A‖ ≤ 1e-12, ‖QᴴQ − I‖ ≤ 1e-12, |R_jj| sequences
     ≤ 1e-12 relative, `IPVT` identical to the HOST `la_qrpp` (use matrices with well-separated column norms;
     document this). `la_geqp3_host`: bit-identical to direct ZGEQP3. Report the measured maxima per test.
2. ctest runs on nividia with: GNU devel noMPI (HOST), IntelLLVM devel noMPI (HOST), PGI noMPI (HOST), PGI noMPI +
   CUDA flags (CUDA). All existing tests still pass in each (39/39 + the new ones).
3. `GPU=none`: `make -j8 all` for GNU MPI: `ALF.out` does not contain `alf_la` symbols (`nm`), and is
   byte-identical to a build of the same tree without the new module **or** differs only in build metadata –
   show what you compared.
4. A small timing table (CUDA, N = 400, 1024): per routine µs per call vs 101's numbers, to catch a missed
   synchronisation or per-call allocation.
- Report: diff summary, routine/call-site table, test table (backend × test × max rel diff), timings.

## nividia rules
Test only; develop on the laptop and rsync to `~/alf-gpu/work/t202/` (`TMPDIR=~/alf-gpu/work/t202/tmp`). Envs:
`~/alf-gpu/env/{gnu,intel,nvhpc}.sh`. ≤ 4 cores, the GPU only while tests run (check `nvidia-smi` first). Runs
≤ 15 min. Clean up build trees at the end. Ticket 201 runs concurrently in its own work dir.

## Out of scope
Calling `alf_la` from ALF routines (203–206); the ONEMKL backend and SYCL shim (302); device kernels for
permutations/scaling; MPI; `configure.sh` (201).

## Resolution (manager, 2026-10-28)
- Delivered: `Libraries/Modules/alf_la_mod.F90` (HOST + CUDA; ONEMKL is an `#error` until 302), 13 ctest programs in
  `testsuite/alf_la.tests/`, helper scripts `docs/gpu/scripts/nividia/t202_*.sh`. Routine → ALF call-site table in
  the module header. Not in `alf_la` (kernels for 203–206): ZLAPMR/ZLAPMT, ZSCAL scalings, D separation, det/phase.
  Not provided: `la_gemm_rc`, `la_getri` (no hot-path call site).
- numerics-reviewer, round 1: NaN-blind checks under fast-math, `t_order` did not test WAR/WAW, NULL-stream choice
  unjustified, `la_qrpp` not compared with variant 1's real path (ZGEQP3 with fixed columns), plus test gaps. Round 2:
  `maxval`/`max` drop NaN, and `t_order` had no negative control. All fixed.
- **Stream finding** (ADR-002 §1 amended, DISCREPANCIES): nvfortran 26.9 runs synchronous target regions on the
  per-thread default stream (`acc_get_cuda_stream(acc_async_sync)`), `nowait` regions on a separate non-blocking
  stream; `!$omp interop` is unimplemented. `alf_la` uses the `acc_async_sync` stream. No `nowait` around `alf_la`
  without `la_sync`; single host thread, one GPU per process.
- Manager re-run on a clean `git archive` tree + the 202 files (nividia):
  - ctest 52/52 with GNU devel, IntelLLVM devel, PGI devel (HOST) and `PGI noMPI GPU=nvidia` (CUDA, A30);
  - `t_order` alone: PASS (8 checks, max 0); with `-DALF_LA_TEST_BADSTREAM`: FAIL 6/8 (RAW 0.41/0.86, WAR 1.37/1.08,
    WAW 0.070/0.070), so the canary works;
  - largest CUDA deviations: gemm 1.4e-15, trsm 3.0e-15, LU 8.2e-16, QR invariants 2.8e-14 (HOST 2.7e-14, rounding);
    HOST bit-identical to direct LAPACK, `la_qrpp` HOST bit-identical to ZGEQP3(jpvt=1) on the pre-permuted matrix;
    `la_geqp3_host` bit-identical on both backends;
  - `GPU=none` GNU MPI `ALF.out`: no `alf_la` symbols (module only in `modules_90.a`).
- Timing (CUDA, agent): close to 101 at N=1024; `la_qrpp`+ungqr at N=400 3.9 ms vs 2.8 ms (host sort of the norms,
  fixed ≈1 ms); `la_geqp3_host` 15 ms (N=400) / 199 ms (N=1024) incl. transfers. Revisit after a sweep profile (210).
- Open: ONEMKL backend and its ordering test (302); multi-GPU binding (405).
