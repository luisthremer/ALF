# 302a – ONEMKL backend of `alf_la` on the Arc B580

- **Phase / lane:** 3 (Intel lane)
- **Owner agent:** fortran-gpu-dev (Sonnet)
- **Review:** numerics-reviewer (code, acceptance file, scripts) · **Verification:** acceptance-verifier
- **Depends on:** 202 (interface, frozen); a consistent oneAPI on `intel` (met 2026-10-30: ifx 2026.2 + oneMKL 2026.2,
  DECISIONS.md row 13); the process diff (merged)
- **Base:** `gpu-port` 7b3c0dc6 (206 state + harness incl. the ctest-only driver fix + decision 13; neither 207 nor 212 is in this base)
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
- `Libraries/Modules/alf_la_mod.F90` (`ALF_LA_ONEMKL` code paths, plus the `#line` directives of the `__LINE__` rule;
  HOST and CUDA otherwise unchanged)
- `Libraries/Modules/alf_la_sycl_shim.cpp` (new; the only non-Fortran source file), `Libraries/Modules/Makefile` (C++
  rule)
- `configure.sh` (intel block: shim and link flags only)
- Existing test files that may change (`allowed_test_changes`), only as stated under "Exact rules":
  `testsuite/alf_la.tests/CMakeLists.txt`, `la_test_util.F90`, `t_order.F90`, `t_info.F90`, `t_fail.F90`
- `docs/gpu/scripts/intel/t302a_*` (development only)
- Not allowed: `testsuite/gpu_reference/acceptance/**` (frozen), HOST and CUDA code paths, `GPU=none` behaviour, new
  Fortran files outside `testsuite/`, `.inc` or other new files under `Libraries/` besides the shim.
- Not allowed: the ONEMKL guards in `gpu_stab_mod`, `gpu_state_mod`, `gpu_wrap_mod` (that is 302b).

## Acceptance (pre-registered in `302a.json` + `302a.nvidia.json` before implementation)
Two parts, per the machine split: the B580 tests run on `intel` (`302a.json`), GPU=none identity (GNU and ifx) and the
unchanged CUDA backend on `nividia` (`302a.nvidia.json`, `acc_driver.py --part nvidia`).

Exact rules the acceptance files depend on (any deviation, e.g. a design without the shim or a different set of
changed build files, means re-registration before T1):
- **`GPU=none` and CUDA identity including `__LINE__`:** the preprocessed HOST form of `alf_la_mod.F90` (9 macro sets of
  `acc_gpunone.sh`) must be identical to the base **including the integers `__LINE__` expands to** (`LA_FAIL` passes
  `__LINE__`; the harness masks only `"X.F90",N` markers). The CUDA `alf_la.tests` log contains `alf_la_mod.F90 (NNN):`
  from `LA_FAIL` and is compared base vs candidate, so **CUDA line numbers stay unchanged too**. cpp ignores directives
  inside skipped groups, so a `#line` inside an ONEMKL-only branch does nothing for HOST/CUDA. What works (checked by
  review (b) with gfortran -E -cpp -P): a `#line N` as the first line of each HOST `#else` (and CUDA) branch that follows
  an ONEMKL branch, and an unconditional `#line N` right after the `#endif` of each ONEMKL-only region, N = the line the
  next line has in the base file. Every `LA_GUARD` line and every comment outside ONEMKL-only regions stays
  byte-identical (comments survive `-P`; an indented macro line expanding to whitespace is a differing line for
  `diff -B`), so the stale "not implemented yet (ticket 302)" comments stay. T0 checks that ifx fpp and nvfortran honour
  `#line` (a forced CUDA/ONEMKL `LA_FAIL` shows a sensible line). Expected counts: 9 comparisons, 0 differing lines,
  0 added/deleted, 3 build files (`configure.sh`, `Libraries/Modules/Makefile`, `testsuite/alf_la.tests/CMakeLists.txt`),
  1 other file (`alf_la_sycl_shim.cpp`).
- **Tests stay the tests:** under ONEMKL the same 13 tests (`t_checks t_gemm t_geru t_trsm t_lu t_qr t_qrpp t_geqp3
  t_order t_fail t_fail_exit t_info_recycle t_info_sync`) build from the same `t_*.F90` files with the same
  registrations as HOST/CUDA; `t_onemkl_guard` goes away; every ONEMKL test runs with `OMP_TARGET_OFFLOAD=MANDATORY`.
  In `la_test_util.F90` only the `#ifdef ALF_LA_CUDA` lines of the `dev_*` helpers may change (to cover ONEMKL) and new
  hook helpers may be added; `TOL_CUDA`, `check*` and `test_end` stay byte-identical. In `t_order.F90` only the backend
  `#ifdef` lines may change. In `t_info.F90` and `t_fail.F90` only the final "reached the end" line changes, to
  `write(*,'(A)') 'FAILED t_info: reached the end without failing'` (resp. `t_fail`), so that a run reaching the end is
  a check failure for the harness (output-only change; the owner is told). Reviewer (a) diffs
  `base..cand -- testsuite/alf_la.tests`.
- **Backend visible:** `la_backend_name()` returns `onemkl` (the log line `TEST t_gemm  backend=onemkl` is gated);
  `la_init` prints ` alf_la: ONEMKL device: <SYCL device name>` (gated: must contain `B580`) and stops if the OpenMP
  default device and the SYCL queue's device are not the same Level Zero device.
- **Test hooks** (read only in `ALF_LA_ONEMKL` code or in `testsuite/alf_la.tests`; when armed, alf_la prints
  ` alf_la test hook armed: <NAME>`; the verifier checks the armed line and that no gate run has an `ALF_LA_TEST_*`
  variable; 302b extends the 212 production refusal to `ALF_LA_TEST_*`). Each control must end as FAIL (check failures
  with `FAILED <name>` lines), never as a crash (exit 2 = INPUT_ERROR):
  - `ALF_LA_TEST_BREAK_ORDER=1`: `la_gemm` submits to a second queue; `la_sync` and `la_finalize` still wait on both
    (like BADSTREAM on CUDA), so nothing is freed while kernels run. T0 reports the failure fraction over 100
    repetitions; zero failures means the control is not killed: stop and report.
  - `ALF_LA_TEST_POISON_INFO=1`: the ONEMKL backend never reads the info slots back; `t_info_recycle` and `t_info_sync`
    must fail (FAILED lines via the changed end line).
  - `ALF_LA_TEST_WRONG_LD=1`: `la_gemm` calls with `opa = 'N'` and `lda > m` pass `lda' = m` (valid, in bounds);
    `t_gemm` must fail its checks.
- **No crash markers:** the C++ shim catches every exception (oneMKL throws on illegal arguments, e.g. `t_fail`'s m < n
  `ungqr`) and returns it as an info/error code, so `terminate called` never appears; no runtime message may contain
  `omptarget` (the harness counts it as a crash). T0 greps the logs for both.
- **Time:** one `t_order` repetition ≤ 8 s on the B580 (100 repetitions within 15 min), measured in T0; otherwise stop
  and report (a different count needs the owner's approval).

1. `alf_la.tests` under `GPU=intel` on the B580: all 13 routine tests pass with unchanged assertions and tolerances;
   controls: poisoned info and wrong ld fail.
2. `t_order` (RAW/WAR/WAW without host waits) passes 100/100 repetitions; control: broken ordering fails.
3. Per-routine table, measured: where it runs (GPU / host / hybrid; method: Level Zero/`LIBOMPTARGET_DEBUG` trace or
   timing vs a host-pinned run), time at N = 144, 400, 1024. **Not a gate:** judged by the owner at the gate; `la_gemm`,
   `la_gemm_cr`, `la_gemm_batched`, `la_geru` and `la_trsm` must be shown to run on the GPU.
4. `GPU=none`: preprocessed `alf_la_mod.F90` identical to the base (gfortran, 9 macro sets, incl. `__LINE__`); GNU ctest
   count unchanged (52); the ifx CPU build (ifx 2025.3.3 on `nividia`, a stand-in, not the 2026.2 pin) gives
   bit-identical output to the base on `hub6x6_cb`. Not gated: ifx-only preprocessor branches and an ifx ctest.
5. `GPU=nvidia` unaffected: `alf_la.tests` under CUDA on `nividia` gives output identical to the base build (timings
   dropped), including the `LA_FAIL` line numbers.

Toolchain record (decision 13, no package hold): the `intel` package set of 2026-10-30 is
`intel:~/alf-gpu/backup/intel_versions_2026-10-30.txt` (sha256 `fd00ecf4…d387`; ifx 2026.2.0-863, MKL core
2026.2.0-881, `intel-opencl-icd`/`libze-intel-gpu1` 26.35.39758.17, `libze1` 1.32.0, `libigc2` 2.41.15). The manager
compares it before every T1/T2 on `intel` and stops on any difference.

## Sequencing with 212
212 (NVIDIA lane, open) also edits `alf_la_mod.F90` (CUDA paths, `k_colnorms`) and `testsuite/alf_la.tests/CMakeLists.txt`
(`t_determinism`). Different blocks, but the same files: whichever ticket lands on `gpu-port` second is rebased onto the
first, its acceptance base moves, and its T2 runs on the rebased candidate (a new pre-registration commit for the base
change only, with the owner's approval at that time). The `alf_la` interface stays frozen.

## Approvals
- First pre-registration `469a2255` (base `5f062f29`) was superseded before any implementation result, after review
  (b) (M1–M7, S1–S8). Current pre-registration commit: `0707248f` (2026-10-30) on `gpu-port-302a-onemkl`; base
  `7b3c0dc6`. Any change of
  `302a.json` / `302a.nvidia.json` after this needs the owner's quoted approval and fresh seeds (VERIFICATION.md §3),
  unless made before any implementation result exists and pre-registered again in a new commit recorded here.

## Evidence
`intel:~/alf-gpu/work/t302a/` and `nividia:~/alf-gpu/work/t302a/`: development runs (T0) under `dev/`; candidates in
`cand/agent/<sha12>/` (T1) and `cand/mgr/<sha12>/` (T2). Never delete raw outputs or binaries before the owner's gate.

## Handoff notes
(implementer, at ~250k tokens or at a phase boundary)

## Out of scope
Activating 204–206 device paths on the B580 (302b), the τ path (302c), PVC (301, 304), performance work.
