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

# Phase 3b (SYCL queue route), steps 1 and 2 (estimate); no production source change yet
GATE A' PASSED. syclbench (OpenMP-offload Fortran, sycl::queue from `omp interop` as om_init, SYCL kernels + oneMKL SYCL gemm from a C++ shared library,
queue is in-order; same chain/checksum as asyncbench; SuperMUC 560361): OpenMP sync 12.3 us/op; SYCL wait after every op 9.2; SYCL async, one wait at the end
4.3 us/op (issue-only 4.1) = 35 % of OpenMP sync; checksums bitwise identical in all variants (also to the OpenMP result). Note: SYCL kernels of our own
only run if the C++ part is linked by the sycl driver (shared lib); an ifx-linked .o gives "No kernel named ... found".
Step 2 (measured, not guessed). L=24, one tile, per sweep: update phase 4.2-4.4 s = 8640 blocks (nd=32; 18 per slice x 480 slices) x ~500 us. Host-side
per-vertex work is negligible (0.18 s/sweep, all site decisions use the host block Gb; NO device read per site). Unavoidable host syncs per sweep: exactly ONE
per block (the D2H of the block gather gbd in blk_begin) = 8640; none per site, none per slice. Timers (P3_TIMERS): blk_begin 2.1-2.8 s/sweep, blk_end 1.3 s,
slice level (hopping wraps) 0.28 s.
KEY FINDING: blk_begin is not launch latency: the 32 KB `target update from(gbd)` costs 380 us (blkbench E1; 16 KB: 10 us, 512 KB: 281 us), the D2H to
PAGEABLE host memory falls off a cliff above the default staging buffer. `LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576` (runtime env var, no source
change) makes it 22 us; in ALF (pb2, 560370, 2 alternating pairs, 3 sweeps): sweep 8.8 -> 6.6 s, update 4.3 -> 1.9 s/sweep (blk_begin 2.6 -> 0.21 s);
stab/cgr/obs unchanged; 1 MB, 4 MB, 16 MB equal; confout_0.h5 and data.h5 identical (h5diff) to the unset run.
Estimate for SYCL async on top of that knob (blkbench: OMP-sync block 188 us with knob, SYCL async 1 wait/block ~140 us): update 1.9 -> ~1.5 s/sweep,
i.e. ~0.4 s/sweep (~6 % of the 6.6 s sweep); the rest of the block is device execution (scale 22, flush gather 2x15, gemm 2x20 + 2x26 us measured with a
wait after every op) that async submission cannot hide because the next block's gather depends on the flush (host decisions in between are ~20 us).
Without the knob, SYCL async would only gain 570 -> 500 us/block (D2H dominates).
