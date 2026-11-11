# Phase 3 (async update phase) status, after step 1: GATE A FAILED, stopped (no source change)

Microbenchmark `asyncbench.F90` (jobs ab1/ab2/ab3 = SuperMUC 560333-560335, 1 tile, ifx 2025.3 AOT pvc, 8 threads): chain of 32 dependent
small kernels (1024 elements) + a dispatched zgemm 1024x32x32 and a copy kernel every 8th, 300 repetitions, 12 us/op unit.
- (a) synchronous (today): 11.9-12.7 us/op (9 runs).
- (b) `nowait depend(inout: tok)` on every target region and `!$omp dispatch nowait depend(...)`, one taskwait at the end: 12.0-12.7 us/op,
  i.e. 97-103 % of (a) (needs <= 50 %). Result bitwise identical to (a) (checksum equal in all runs, order preserved). The "issue only" time
  (before the taskwait) equals the total: the encountering thread runs each nowait target task to completion; with real host work between
  launches (sin loop, +12 us/op at hw=1000) the totals are sync 24.1 vs async 23.7 us/op: no overlap of host work with device work either.
- (b') nowait without depend: 34 us/op (slower) and checksum differs from run to run (races): not usable.
- (c) knobs, each x {sync, async}, all 11.9-12.7 us/op: LIBOMPTARGET_LEVEL_ZERO_COMMAND_MODE=sync|async|async_ordered, COMMAND_BATCH=copy|compute|all|all,64
  (+ COMMAND_MODE=async), NO_SYCL_FLUSH=1, USE_MULTIPLE_COMPUTE_QUEUES=1, LIBOMP_NUM_HIDDEN_HELPER_THREADS=1, OMP_NUM_THREADS=1: no effect.
  USE_IMMEDIATE_COMMAND_LIST=0: 28-29 us/op (2.4x worse, with or without COMMAND_MODE=async).
Conclusion: with ifx 2025.3 / Level Zero RTL the OpenMP `nowait` mechanism does not defer anything; B3 as specified (nowait + dependencies) gains nothing.
The previous latbench value (empty target nowait + taskwait 7.6 us vs sync 6.0 us) already pointed that way. No source changes were made.
