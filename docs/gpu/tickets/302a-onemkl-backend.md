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
- **Time:** no time limit (decision 16, 2026-10-30, replacing decision 15's per-repetition limit): run times are not
  acceptance criteria of this ticket. History: first registered as ≤ 8 s per `t_order` repetition (100 within 15 min);
  decision 15 allowed ~19-min runs with ≤ 12 s per repetition; the second T0 (presence guard) had mean 11.50 s, min
  9.80 s, max 12.98 s, 22 of 100 above 12 s; implementer 3 traced it to a fixed ~0.45 ms per oneMKL `zgeqrf` call (8292
  8×8 calls per pass in the info-slot loop), see Handoff notes and PLAN.md ticket 306. A run that shows no progress for
  15 min is investigated (CLAUDE.md). Expect ~20 min per 100-repetition `t_order` run.

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
- Third pre-registration (after T0 results, with the owner's quoted approval, DECISIONS.md row 15, 2026-10-30): only
  the `t_order` time rule above changes; `302a.json` and `302a.nvidia.json` are unchanged. Commit `ef30771b`; the owner's
  per-repetition reading ("Each rep") was added in `cf31bc8b`. Decision 16 (no time limit) was pre-registered in
  `84dade95`, the current pre-registration.

## Evidence
`intel:~/alf-gpu/work/t302a/` and `nividia:~/alf-gpu/work/t302a/`: development runs (T0) under `dev/`; candidates in
`cand/agent/<sha12>/` (T1) and `cand/mgr/<sha12>/` (T2). Never delete raw outputs or binaries before the owner's gate.

## Handoff notes
(implementer, at ~250k tokens or at a phase boundary)

**Implementer 2 (2026-10-30, T0 of the candidate with the presence guard; uncommitted tree, T1 not yet run).**

- **Host-side presence guard** (ONEMKL module `alf_la_om_mod` only; HOST/CUDA text and line numbers untouched, ppcheck 0 lines).
  Before every `!$omp dispatch` (`om_gemm`, `om_gemm_cr`, `om_gemm_batched`, `om_geru`, `om_trsm`) and before every shim call
  that receives device addresses (`om_geqrf`, `om_ungqr`, `om_ormqr`, `om_getrf`, `om_getrs`; also `om_qrpp`, whose target regions
  would otherwise map implicitly) the first and the last element touched by every array argument are checked on the host with
  `omp_target_is_present` (default device). A miss calls `om_fail` (same format as the other failures): ` alf_la: la_gemm:
  argument A is not present on the device, nothing was queued`, then `Terminate_on_error` (error stop). The guard runs before the
  BREAK_ORDER branch too. Checked extents follow the BLAS/LAPACK argument shapes (`ext2`, `extv`, batched strides); empty
  extents are skipped. Internal scratch (`dws`) is checked as well.
- **Guard test (dev script, not a ctest):** `docs/gpu/scripts/intel/t302a_guard.{sh,F90}`. (1) `present`: `omp_target_is_present`
  alone: mapped A = 1, mapped A(n,n) = 1, unmapped B = 0, A after `target exit data delete` = 0. (2) code read: `guard_gemm` is the first
  statement after the opcode translation in `om_gemm`, before the hook branch and the dispatch; `om_fail` ends in `error stop`.
  (3) `pos`: all mapped, C(1,1) = reference. (4) ONE negative call (`neg A`, unmapped A through `la_gemm`, load checked, no other
  user): exit code 1 with the message above, no `omptarget`, no crash, machine unaffected. Output: `intel:~/alf-gpu/work/t302a/dev/guard_out/{present,pos,neg_A}.out`.
- **BREAK_ORDER mechanism:** `la_gemm` under the hook calls `alfsh_zgemm_unordered`, which submits the GEMM to a SECOND sycl queue
  (same context/device, created lazily) and returns without waiting. Two queues of one device ran in submission order on the B580
  (0 of 5 runs failed), which hid the missing ordering, so every unordered GEMM also depends on a SYCL host task (the gate) that
  spins until `alfsh_wait` (called by `la_sync` / `la_finalize`) releases it and then waits for the second queue. The GEMMs
  therefore run late, after the caller's following target regions, exactly the effect of a missing dependency; nothing is
  freed while they run (the wait is in the sync/finalize path). It is a faithful ordering control: the only change versus the
  correct path is that the dependency between the caller's stream and the GEMM is cut. Result: 100/100 repetitions FAIL.
- **Per-routine table** (measured, `intel:~/alf-gpu/work/t302a/dev/t0g/table/table.txt`, generated by `t302a_table.py` from the raw
  `table_*.txt`; method: MKL_VERBOSE trace, device time vs one-host-core pinned run): GPU: `la_gemm`, `la_gemm_cr`,
  `la_gemm_batched`, `la_geru`, `la_trsm`, `la_getrf`, `la_getrs`; hybrid (device call, host work): `la_geqrf`, `la_ungqr`,
  `la_ormqr`, `la_qrpp` (its geqrf); host (ZGEQP3 + D<->H/H<->D copies): `la_geqp3_host`.
- **Dev scripts** (all in `docs/gpu/scripts/intel/`, `t302a_*`): `sync` (rsync to `dev/src` on a machine), `build`, `ctest`, `reps`
  (t_order repetitions), `table`/`table.py`, `guard`, `gnu`, `ifxid` + `ifxrun` (+ `dev/h5c.py` on nividia compares all datasets),
  `cuda`, `ppcheck`. Pitfall found: sourcing `~/alf-gpu/env/intel.sh` resets `"$@"`: scripts save the arguments first.
- **Raw outputs:** `intel:~/alf-gpu/work/t302a/dev/t0g/` (A1 `a1_ctest.log`, controls `ctl_*.log`, `order_pass/`, `order_broken/`
  with per-repetition logs and `reps.txt`, `table/`), `dev/guard_out/`; older T0: `dev/t0/`. nividia: `dev/gnu_run2.log`,
  `dev/ifxrun2.log`, `dev/ifxid_h5c2.out`, `dev/cuda_*`.
- **Open:** `t_order` repetitions vary 9.8 to 13.0 s (mean 11.50, 22 of 100 above 12 s); see the report to the manager.

**Implementer 3 (2026-10-30, time profile of `t_order`; STOPPED for a design decision, nothing implemented; tree unchanged except dev scripts).**

- **Profile** (`t302a_tprof.{F90,sh}`, a timing copy of `t_order`; `t302a_micro.*`; `t302a_qrbench.{cpp,sh}`; build `intel:~/alf-gpu/work/t302a/dev/build_g1`):
  one repetition = 11.0 s in main: both passes together `la_init` 0.12 s, JIT of the first gemm/target regions ~0.05 s, first geqrf 0.74 s,
  RAW+WAR+WAW loops 1.0 s per pass (2.0 s), and the **"info slot traffic" loop of 2*LA_NSLOT+100 = 8292 `la_geqrf(8,8)` calls per pass
  = 3.9 s and 4.1 s (7.99 s, 73 percent)**. Process start-up and runtime init are small (`t_checks`, which never inits: 0.17 s;
  `t_gemm` 3.1 s). Per tiny `la_geqrf` 0.45 to 0.8 ms; the pieces (`t302a_micro`): `omp_target_is_present` guard 0.05 us each,
  `target data use_device_addr` 2 us, scratchpad query 0.4 us, **oneMKL `zgeqrf` itself 400 to 450 us**. Cause: oneMKL's hybrid geqrf, a
  fixed cost per call independent of the queue (in-order or not, waiting or not: 8x8 447/469 us, 64x64 540/466 us, 144x144 786/784 us),
  of `MKL_NUM_THREADS=1`, `KMP_BLOCKTIME=0`. The 3 s spread between repetitions is the jitter of these ~16600 calls (0.45 to 0.8 ms each);
  guard, shim queue set-up and JIT do not show up. `geqrf_batch` (strided, batch 1) is not an alternative: wrong results (max diff 0.4 at n=8).
- **Measured option (needs the owner's/manager's decision, because it moves `la_geqrf` for small sizes from hybrid to host with two
  USM copies):** in the shim, `m,n <= 64`: copy A to the host, MKL host `zgeqrf`, copy A and tau back. Micro (`t302a_qrbench`): 8x8 19 to 25 us
  vs 600 us (oneMKL), 64x64 37 us vs 570 to 710 us, 144x144 400 us vs 790 to 920 us, 400x400 2.2 to 2.4 ms vs 7.7 to 7.9 ms. End to end
  (dev experiment, reverted, build `dev/build_exp1`, not evidence): `t_order` 2.91 s and 2.74 s, `t_qr` 1.25 s. Other options: ahead-of-time
  compilation does not help (JIT is about 1 s of 11); a custom tiny-QR SYCL kernel (more code and numerics to verify).

**Implementer 4 (2026-10-30, review (a) findings of candidate 1; uncommitted tree).** Details of what changed and of the T0 runs are in the report to the manager; durable notes:

- **What `t_order` shows.** `ALF_LA_TEST_BREAK_ORDER` is a maximal-reordering control: every GEMM is held back by a gate until `la_sync` /
  `la_finalize`, so it shows that the test detects gross reordering, not partial overlap or a short race window. A PASS of `t_order` on the
  correct path shows host-synchronous ordering (every shim call waits, dispatched BLAS calls are enqueued on the OpenMP runtime's queue and
  the caller's target regions follow), not queue-level dependencies: the backend declares none and relies on the host-synchronous design.
- **BLAS dispatch and ADR-001 §2.** The BLAS calls use `!$omp dispatch` with the host-side presence guard (`need_*`, `omp_target_is_present`),
  not `target data use_device_addr` around each call as the ADR §2 wording says; the dispatch maps the (present) arrays itself. Deviation
  recorded for the ADR; `use_device_addr` is used only around shim calls.
- **`-I${MKLROOT}/include`** reaches ANA/PROG as well as `Libraries/Modules` (`configure.sh`'s intel block adds it to `ALF_FLAGS_MODULES`, `_ANA` and `_PROG` alike). Harmless
  (no MKL module is used there); left as is.
- **`om_fail` tag** is now `alf_la_mod.F90[onemkl]` (checked: no ctest or acceptance regex reads the file name of ONEMKL failures; the HOST/CUDA
  `LA_FAIL` text is unchanged).

- **Singular `la_getrf` on the B580** (`t302a_singular.{F90,sh}`, `intel:~/alf-gpu/work/t302a/dev/t0h/singular.out`): an 8x8 matrix with an
  exactly zero column runs through oneMKL `zgetrf` without exception or crash; `la_getrf`/`la_sync` return normally (info > 0 tolerated as on
  HOST), U(3,3) = 0. No `terminate`, no `omptarget`.

## Out of scope
Activating 204–206 device paths on the B580 (302b), the τ path (302c), PVC (301, 304), performance work.
