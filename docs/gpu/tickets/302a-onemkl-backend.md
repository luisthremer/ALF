# 302a – ONEMKL backend of `alf_la` on the Arc B580

- **Phase / lane:** 3 (Intel lane)
- **Owner agent:** fortran-gpu-dev (Sonnet)
- **Review:** numerics-reviewer (code, acceptance file, scripts) · **Verification:** acceptance-verifier
- **Depends on:** 202 (interface, frozen); a consistent oneAPI on `intel` (matching ifx and MKL release); the owner's
  approval of `gpu-port-process-lanes`
- **Status:** draft (2026-10-29), not started

## Goal
Every public `alf_la` routine has an `ALF_LA_ONEMKL` implementation that runs on the B580, passes the `alf_la.tests`
routine tests and the ordering canary `t_order`, and is documented per routine as GPU / host / hybrid (measured).

## Context
- `Libraries/Modules/alf_la_mod.F90`: `ALF_LA_ONEMKL` today is a guard only (`LA_GUARD` in every public routine,
  `la_init` stops). Public interface (frozen): `la_init, la_finalize, la_sync, la_backend_name, la_reserve_work,
  la_gemm, la_gemm_cr, la_gemm_batched, la_geru, la_trsm, la_geqrf(_lwork), la_ungqr(_lwork), la_ormqr(_lwork),
  la_getrf(_lwork), la_getrs, la_qrpp, la_geqp3_host, la_debug_poison_last_info`. There is no `la_getri`
  (PLAN/ADR-002 corrected 2026-10-29).
- ADR-001 §2, ADR-002 §1: `!$omp dispatch` for BLAS (one call per dispatch, inside `target data use_device_addr`);
  a small C++ SYCL USM shim for LAPACK pieces that dispatch would run on the host. **Open question (ADR-002 §1):** the
  queue — OpenMP `interop` `sycl::queue` (ifx 2025.3 `omp_lib` lacks `omp_get_interop_ptr`, declare it `bind(C)`;
  103 prototype, not verified on a GPU) or a Level Zero `make_queue` fallback. `t_order` decides.
- oneMKL 2026.2 headers on `intel` (2026-10-29, read-only grep): offload interfaces exist for zgeqrf, zunmqr, zungqr,
  zgetrf, zgetrs, zgetri, ztrsm, zgemm, zgesvd; **none for zgeqp3** → STABLOG keeps host ZGEQP3 (`la_geqp3_host`) as on
  NVIDIA. Docs say getrf (n ≤ 256 on the CPU), geqrf and unmqr are hybrid under OpenMP offload; measure it.
- Prototype: `gpu_bench/src/bench_la.F90` (`-DBENCH_LA_ONEMKL`, `-DBENCH_LA_SHIM`), `gpu_bench/src/bench_sycl_shim.cpp`,
  `gpu_bench/README.md:69–92`. Runtime rule: ifx, `libomptarget`, `libsycl` and `libmkl_sycl*` from the **same**
  oneAPI release (2025.3 with 2026.2 fails at link time with SYCL symbols, 207).
- B580: single tile, FP64 present (rate unmeasured); results are correctness evidence, not PVC performance.
  `ONEAPI_DEVICE_SELECTOR=level_zero:0`; never plain `sycl-ls`.

## Files allowed
- `Libraries/Modules/alf_la_mod.F90` (`ALF_LA_ONEMKL` code paths only; HOST and CUDA unchanged)
- `Libraries/Modules/alf_la_sycl_shim.cpp` (new), `Libraries/Modules/Makefile` (C++ rule)
- `configure.sh` (intel block: shim and link flags only)
- `testsuite/alf_la.tests/` (enable routine tests and `t_order` under ONEMKL)
- `testsuite/gpu_reference/acceptance/302a.json`, `docs/gpu/scripts/intel/t302a_*`
- Not allowed: the ONEMKL guards in `gpu_stab_mod`, `gpu_state_mod`, `gpu_wrap_mod` (that is 302b).

## Acceptance (pre-registered in `302a.json` before implementation; see `302a.draft.json`)
1. `alf_la.tests` under `GPU=intel` on the B580: every routine test passes with the same tolerances as HOST/CUDA
   (unchanged test assertions); positive control: the poisoned-info and wrong-ld fault tests fail as on CUDA.
2. `t_order` (RAW/WAR/WAW without host waits) passes 100/100 repetitions; positive control: the test with the ordering
   deliberately broken (no sync, separate queue) fails.
3. Per-routine table, measured: where it runs (GPU / host / hybrid; method: Level Zero/`LIBOMPTARGET_DEBUG` trace or
   timing vs a host-pinned run), time at N = 144, 400, 1024.
4. `GPU=none` (GNU and ifx CPU): preprocessed `alf_la_mod.F90` identical to the base for every non-ONEMKL macro set;
   ctest count unchanged.
5. `GPU=nvidia` unaffected: `alf_la.tests` under CUDA on `nividia` bit-identical results to the base build.

## Out of scope
Activating 204–206 device paths on the B580 (302b), the τ path (302c), PVC (301, 304), performance work.
