# Phase 4 (device-only ZGEQRF in the SYCL shim) - status 2026-11-12

## Gate A: PASSED (final kernel)
- Code: `Libraries/Modules/alf_la_sycl_devqr.hpp` + `alf_la_sycl_devqr.cpp` (-> `libalfdevqr.so`, dlopen'ed by `alf_la_sycl_shim.cpp` when
  `ALFSH_DEVQR=1`; the SYCL device code cannot live in the plain host shim object that ifx links; `configure.sh` got `-ldl`). Default off: production unchanged.
  Panel = ONE work-group, every panel column in the registers of G sub-groups (16 lanes, large GRF), ZGEQR2 (LAPACK ZLARFG semantics incl. safmin
  rescaling, look-ahead of the norm) + ZLARFT in the same kernel; trailing update ZLARFB = 3 oneMKL device GEMMs on an explicit V.
- `devqr_test.cpp` vs host MKL zgeqrf (N=36/144/576/1024, non-square, tiny; random, column scales 1e-12..1e12 sorted/shuffled, 1e-150..1e150):
  dR <= 2.1e-14, dV <= 5.9e-14, dtau <= 2.2e-14, orth ||Q^H Q-I|| <= host; planted tau(0)*=1+1e-12 FAILS (dtau 1e-12). Cond-1e24 matrices: R equal, reflectors differ (noise): informational.
- Host MKL zgeqrf leaves the last 1-row reflector as identity (tau=0, R(k,k) complex); reference LAPACK and device oneMKL (what ALF uses) make it real; devqr = the latter.
- Standalone SYCL programs need SYCL_UR_USE_LEVEL_ZERO_V2=1 here and then (sanity.cpp) a device->host copy completes ms AFTER wait() returns, and
  event waits return early (standalone timings = submission time only; hybrid oneMKL zgeqrf gives wrong results). All timing below is inside ALF.

## Gate B: FAILED (stopped here)
Per call inside ALF (shim timers, 1 tile, delta of two runs): oneMKL hybrid 9.1 ms (8 thr) / 26.3 ms (1 thr) at N=576, 26.3 ms (8 thr) at N=1024.
Device QR (final kernel): N=576 nb=8/12/16/32: 9.8 / 9.3 / 8.7 / 9.5 ms (1 thread, nb=16: 8.7 ms: independent of host threads); N=1024 nb=16/32: 19.9 / 21.3 ms.
Gate (<= 0.5x of oneMKL 8 thr) = 4.7 / 14.5 ms: NOT met (1.9x / 1.4x too slow); gain vs oneMKL only with few host threads (3x at 1 thread) and at N=1024 (1.3x).
Cause (probes inside ALF, `ALFSH_LAUNCHBENCH=1`): launches are cheap (3 us kernel, 8 us small GEMM), the panel chain is the limit: one work-group = one Xe core,
>= 7 us per column step (barrier, ZLARFG fp64 division chain, SLM broadcast of v to all sub-groups = SLM bandwidth), 15-36 us per step for 40 register slots
(576 rows) where the kernel spills; 576 steps -> 4-8 ms before the 3 trailing GEMMs. Global-memory panel (v1): 38 us/step (memory latency), 20 ms/call.
Not tried: panel spread over several work-groups with grid-wide sync, TSQR-like panel (changes the algorithm), 2-queue look-ahead.
Checks with the final kernel (devqr on): 3a same-seed vs CPU build: L6 Mz=F 2.2e-15, L6 Mz=T 6.0e-11, L12 Mz=F 1.6e-14, confout_0 and acceptance identical;
3c planted defect (tau(1)*=1+1e-12 in a throwaway lib): FAIL 4.1e-13 / 6.0e-13. 3b Precision Green (info): L24 mean 1.4-1.9e-14 / max 1.4-4.5e-11 (baseline 1.46e-14 / 3.07e-11), L32 mean 3.9-4.0e-14 / max 4.4-6.1e-11 (baseline 9.9e-14 / 1.2e-10).
