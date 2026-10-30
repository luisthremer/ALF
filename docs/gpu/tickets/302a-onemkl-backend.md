# 302a – ONEMKL backend of `alf_la` on the Arc B580

- **Phase / lane:** 3 (Intel lane)
- **Owner agent:** fortran-gpu-dev (Sonnet)
- **Review:** numerics-reviewer (code, acceptance file, scripts) · **Verification:** acceptance-verifier
- **Depends on:** 202 (interface, frozen); a consistent oneAPI on `intel` (met 2026-10-30: ifx 2026.2 + oneMKL 2026.2,
  DECISIONS.md row 13); the process diff (merged)
- **Base:** `gpu-port` 5f062f29 (206 state + harness + pin change; neither 207 nor 212 is in this base)
- **Branch:** `gpu-port-302a-onemkl`
- **Acceptance files:** `testsuite/gpu_reference/acceptance/302a.json` (Intel part, `intel`) and `302a.nvidia.json`
  (`nividia` part, `acc_driver.py --part nvidia`), pre-registered in their own commit before any implementation (sha
  under Approvals; DECISIONS.md pinned at that commit)
- **Status:** open (2026-10-30; owner: "start 302a per decision 11")

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
  queue — OpenMP `interop` `sycl::queue` (ifx 2025.3 `omp_lib` lacked `omp_get_interop_ptr`; check ifx 2026.2, else declare it `bind(C)`;
  103 prototype, not verified on a GPU) or a Level Zero `make_queue` fallback. `t_order` decides.
- oneMKL 2026.2 headers on `intel` (2026-10-29, read-only grep): offload interfaces exist for zgeqrf, zunmqr, zungqr,
  zgetrf, zgetrs, zgetri, ztrsm, zgemm, zgesvd; **none for zgeqp3** → STABLOG keeps host ZGEQP3 (`la_geqp3_host`) as on
  NVIDIA. Docs say getrf (n ≤ 256 on the CPU), geqrf and unmqr are hybrid under OpenMP offload; measure it.
- Prototype: `gpu_bench/src/bench_la.F90` (`-DBENCH_LA_ONEMKL`, `-DBENCH_LA_SHIM`), `gpu_bench/src/bench_sycl_shim.cpp`,
  `gpu_bench/README.md:69–92`. Runtime rule: ifx, `libomptarget`, `libsycl` and `libmkl_sycl*` from the **same**
  oneAPI release (2025.3 with 2026.2 fails at link time with SYCL symbols, 207). On `intel`: `. ~/alf-gpu/env/intel.sh`
  (pinned components, ifx 2026.2.0, MKL 2026.2; MACHINES.md).
- B580: single tile, FP64 present (rate unmeasured); results are correctness evidence, not PVC performance.
  `ONEAPI_DEVICE_SELECTOR=level_zero:0`; never plain `sycl-ls`.

## Files allowed
- `Libraries/Modules/alf_la_mod.F90` (`ALF_LA_ONEMKL` code paths only; HOST and CUDA unchanged)
- `Libraries/Modules/alf_la_sycl_shim.cpp` (new), `Libraries/Modules/Makefile` (C++ rule)
- `configure.sh` (intel block: shim and link flags only)
- `testsuite/alf_la.tests/CMakeLists.txt`, `la_test_util.F90`, `t_order.F90` (the only existing test files that may
  change: `allowed_test_changes`; test assertions and tolerances stay unchanged); new files in `testsuite/alf_la.tests/`
- `docs/gpu/scripts/intel/t302a_*` (development only)
- Not allowed: `testsuite/gpu_reference/acceptance/**` (frozen), HOST and CUDA code paths, `GPU=none` behaviour.
- Not allowed: the ONEMKL guards in `gpu_stab_mod`, `gpu_state_mod`, `gpu_wrap_mod` (that is 302b).

## Acceptance (pre-registered in `302a.json` + `302a.nvidia.json` before implementation; see the two drafts)
Two parts, per the machine split: the B580 tests run on `intel` (`302a.json`), GPU=none identity (GNU and ifx) and the
unchanged CUDA backend on `nividia` (`302a.nvidia.json`). Controls need three test-only switches that 302a adds to
`testsuite/alf_la.tests`: `ALF_LA_TEST_BREAK_ORDER=1` (ordering deliberately broken), `ALF_LA_TEST_POISON_INFO=1`
(poisoned device info), `ALF_LA_TEST_WRONG_LD=1` (wrong leading dimension); each must make its test fail.
Exact rules the acceptance files depend on:
- **`GPU=none` identity including `__LINE__`:** the preprocessed HOST form of `alf_la_mod.F90` (9 macro sets) must be
  identical to the base **including the integers `__LINE__` expands to** (`LA_FAIL` passes `__LINE__`; the harness masks
  only `"X.F90",N` markers). ONEMKL code inserted above HOST code must not shift HOST lines: keep the line count of each
  `ALF_LA_ONEMKL` block unchanged (e.g. one `#include "alf_la_onemkl_<routine>.inc"` line replacing the one `LA_GUARD`
  line; new `.inc` files are allowed) or restore the numbering with a `#line` directive inside the ONEMKL-only region.
  Expected counts: 9 comparisons, 0 differing lines, 3 build files (`configure.sh`, `Libraries/Modules/Makefile`,
  `testsuite/alf_la.tests/CMakeLists.txt`), 1 other file (`alf_la_sycl_shim.cpp`); `.inc` files under
  `Libraries/Modules/` count as "other" too, so **use at most the shim as a non-Fortran file, or put ONEMKL bodies in
  `.F90`-named include files** (`alf_la_onemkl_*.F90` are added CPU-path Fortran and count as added_or_deleted:
  not allowed). Simplest: keep the ONEMKL bodies inside `alf_la_mod.F90` and compensate with `#line`.
- **Test hooks:** `ALF_LA_TEST_BREAK_ORDER`, `ALF_LA_TEST_POISON_INFO`, `ALF_LA_TEST_WRONG_LD` are read only in
  `ALF_LA_ONEMKL` code or in `testsuite/alf_la.tests`; when armed, alf_la prints ` alf_la test hook armed: <NAME>` on
  stdout. (302b extends the 212 production refusal to `ALF_LA_TEST_*`.)
- The ONEMKL ctest registers the same 13 tests as HOST/CUDA (`t_checks t_gemm t_geru t_trsm t_lu t_qr t_qrpp t_geqp3
  t_order t_fail t_fail_exit t_info_recycle t_info_sync`); `t_onemkl_guard` goes away; every ONEMKL test runs with
  `OMP_TARGET_OFFLOAD=MANDATORY`.

1. `alf_la.tests` under `GPU=intel` on the B580: every routine test passes with the same tolerances as HOST/CUDA
   (unchanged test assertions); positive control: the poisoned-info and wrong-ld fault tests fail as on CUDA.
2. `t_order` (RAW/WAR/WAW without host waits) passes 100/100 repetitions; positive control: the test with the ordering
   deliberately broken (no sync, separate queue) fails.
3. Per-routine table, measured: where it runs (GPU / host / hybrid; method: Level Zero/`LIBOMPTARGET_DEBUG` trace or
   timing vs a host-pinned run), time at N = 144, 400, 1024.
4. `GPU=none` (GNU and ifx CPU): preprocessed `alf_la_mod.F90` identical to the base for every non-ONEMKL macro set;
   ctest count unchanged.
5. `GPU=nvidia` unaffected: `alf_la.tests` under CUDA on `nividia` bit-identical results to the base build.

## Sequencing with 212
212 (NVIDIA lane, open) also edits `alf_la_mod.F90` (CUDA paths, `k_colnorms`) and `testsuite/alf_la.tests/CMakeLists.txt`
(`t_determinism`). Different blocks, but the same files: whichever ticket lands on `gpu-port` second is rebased onto the
first, its acceptance base moves, and its T2 runs on the rebased candidate (a new pre-registration commit for the base
change only, with the owner's approval at that time). The `alf_la` interface stays frozen.

## Approvals
- Pre-registration commit: `469a2255` (2026-10-30) on `gpu-port-302a-onemkl`; base `5f062f29`. Any change of
  `302a.json` / `302a.nvidia.json` after this needs the owner's quoted approval and fresh seeds (VERIFICATION.md §3),
  unless made before any implementation result exists and pre-registered again in a new commit recorded here.

## Evidence
`intel:~/alf-gpu/work/t302a/` and `nividia:~/alf-gpu/work/t302a/`: development runs (T0) under `dev/`; candidates in
`cand/agent/<sha12>/` (T1) and `cand/mgr/<sha12>/` (T2). Never delete raw outputs or binaries before the owner's gate.

## Handoff notes
(implementer, at ~250k tokens or at a phase boundary)

## Out of scope
Activating 204–206 device paths on the B580 (302b), the τ path (302c), PVC (301, 304), performance work.
