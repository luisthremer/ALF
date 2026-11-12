# Phase 5 (faster device LAPACK on Max 1550) - status 2026-11-12 ~13:10 (STOPPED: ssh supermuc2 hangs, owner must log in)

All numbers C = measured on one PVC tile (job 560715, a dev-loop daemon, standalone harness `qr2_test.cpp`: true device time of a chain of calls closed
by a kernel that sets a host-USM flag, host spins on it; the old device QR measures 7.2 ms at N=576 in this harness, 7.8 ms inside ALF, so the harness
is a good proxy). NOT yet measured inside ALF (Gate B/C): the ALF build + jobs could not be started (ssh hang, see below).

## Cause found (task 1, C)
- Phase-4 panel: 6.9-7.2 us per column even at M=128. Ablations of the new panel (compile-time macros) at M=64, 16 columns: scalar chain (sqrt + 3 divisions)
  ~1.2 us/col, lane reductions ~0.7, the step's synchronisation chain (barrier + the rest) ~2.8; a work-group barrier ITSELF costs 0.16 us (bar_test.cpp),
  an SLM flag handoff 0.6 us (worse). The old kernel's extra cost was the ZLARFT T recurrence on sub-group 0 + Z dot products of idle sub-groups.
- Panel ~90 % of the old device geqrf (integral over the 36 panels); the 3 trailing GEMMs are small.
- SLM bandwidth bounds the kernels: an 8-byte-per-lane load at 16-byte stride (re, im separately) occupies a 256-byte pass of the SLM pipe; one 16-byte (double2)
  load per slot, the reflector kept in registers for the second pass, 2 columns per sub-group, G sub-groups per column (row split) fixed most of it.
  Global loads under branches serialise (~1 us each): branch-free clamp + select.  Work-groups are limited to 32 sub-groups (512 work-items) in 256-GRF mode
  (G=3 with jb=16 fails to launch).  sycl::group_barrier(wg) has a global fence; it.barrier(local_space) is cheaper (0.16 vs 0.23 us).

## What was built (files, uncommitted except as noted)
- Libraries/Modules/alf_la_sycl_qr2.hpp (NEW): panel2_t (T-free ZGEQR2 panel, one work-group, G sub-groups per column, ZLARFG semantics + safmin rescaling,
  look-ahead norm, fast path for the scalar chain), apply2_t (reflector application, sub-group(s) per NV columns, SLM-staged double-buffered reflector chunks,
  no T, no GEMM; used for the QR trailing update and for ZUNMQR; side R by conj-transposed rows), geqrf2() (nb=16, panels M<=256 use 32), unmqr2().
- alf_la_sycl_devqr.cpp (entry points alfdevqr_zgeqrf2, alfdevqr_zunmqr; committed by the lead in ac3b7f8a as work in progress), alf_la_sycl_shim.cpp
  (ALFSH_DEVQR=2 -> new geqrf; ALFSH_DEVUNMQR=1 -> device unmqr with fallback to oneMKL when declined; ALFSH_DEVQR=1 unchanged = phase 4; ALFSH_QRNB default 16 at level 2),
  Libraries/Modules/Makefile (libalfdevqr.so: depends on the new header, uses ALF_SYCL_BLK_FLAGS = AOT for PVC when set).
- p5/: qr2_test.cpp (Gate A + timing), bar_test.cpp (barrier micro-benchmark), daemon.sh/q.sh/sync.sh/env.sh/t_*.sh/r*.sh (dev loop on the SuperMUC),
  bld.sh, p5lib.sh, build_comb.sh (ALF builds + job helpers, prepared, NOT run).

## Gate A (standalone, vs host MKL zgeqrf/zunmqr; C)  - PASSED for the final kernels except the planted-defect runs (not yet run: see open)
- geqrf: 36/144/576/1024 square, types random / column scales 1e-12..1e12 sorted+shuffled / 1e-150..1e150 / zero+real columns; non-square (576x144, 144x576,
  1000x37, 37x1000, 1300x90, 1536x64), tiny (1x1, 5x5, 17x33); nb 8/16/32: dR <= 2.2e-14, dV <= 6.6e-14, dtau <= 2.8e-14, orthogonality <= host (1.5e-15).
- unmqr L/R x N/C, (36..1024)^2 and rectangular/partial-k cases: reldiff vs host zunmqr <= 2.3e-15. Q^H Q - I <= 1.3e-15, ||QR-A||/||A|| <= 2.1e-15 (n=64/576/1024).
- Planted defects in the test (DEVQR_PLANT=1 tau(0)*=1+1e-12, 2 wrong conj/trans, 3 reflector scaled): implemented, NOT yet run.

## Timing (standalone, 1 tile, per call; C)
| N | old device QR (ALFSH_DEVQR=1) | new geqrf2 | new unmqr (L,N / R,C) | oneMKL in ALF (job 560667) |
| 576 | 7.2 ms | 3.9 ms (nb 16, 32 for M<=256) | 1.34 / 1.35 ms | geqrf 10.9, unmqr 2.6 |
| 1024 | 16.7 ms | 10.4 ms | 5.3 / 5.3 ms | geqrf 26, unmqr 8.3 |
Gate B target (inside ALF, ALFSH_TIMING, no ASYNCBLK): geqrf <= 0.6 x 7.8 = 4.7 ms (standalone 3.9: expected met), unmqr <= 0.5 x oneMKL = 1.3 ms at N=576
(standalone 1.34: marginal), 4.1 ms at N=1024 (5.3: NOT met in standalone).

## Open / next steps
1. ssh supermuc2 hangs (timeout 124, twice; other sessions' polling ssh also hang): the ALF build (bld.sh), Gate B/C jobs, planted-defect and node
   Precision-Green runs were not started. Sources in p5/src_alf on the SuperMUC may be incomplete (rsync of the tree was interrupted).
2. unmqr at N=576 needs ~3% and at N=1024 ~30% more: the per-reflector step is latency bound (~1.7-2.3 us/step at nq 576); idea: block 2-4 reflectors per
   reduction round with a per-chunk Gram matrix; apply launch fixed cost ~20 us x 36 per geqrf.
3. geqrf: panel 2.5 ms of 3.9; overlapping the single-Xe-core panel with the wide trailing apply on a 2nd queue (not done), 2 columns per sub-group in the panel.
4. Task 4 (device ZGETRF) not started.

## UPDATE after ssh recovered (jobs 560715 daemon, 560721 build, 560722 g1, 560723 g2, 560724 g5, 560725 g6, 560726 g3a, 560727 g3b, 560730 g3c; scripts g*.sh, p5lib.sh here)
- Gate A planted defects (qr2_test, job 560715): PLANT=1 FAIL (dtau 1e-12), PLANT=2 FAIL (reldiff 1.4), PLANT=3 FAIL (1e-12); clean: 0 failed.
- Default build (no switches) vs comb build, bitwise L6F/L12F: PASS (0.0) (560722).
- Gate B in ALF (ALFSH_TIMING, no ASYNCBLK, marginal 4 vs 2 sweeps; 560722), ms/call: L24 geqrf V3(devqr1) 5.48, V5(devqr2) 3.91, V6(+unmqr_dev) 3.98; unmqr oneMKL 2.41/2.61, unmqr_dev 1.37;
  getrf 3.7-4.0. L32 geqrf V3 16.95, V5 10.36, V6 10.34; unmqr oneMKL 7.63/7.44, unmqr_dev 5.30.
- Gate C same seed: L6F 2.3e-15, L12F 1.6e-14 PASS; L6T 1.358e-10 FAIL (tol 1e-10); noise study (g4, g5): oneMKL itself fails 2 of 5 seed sets (1.13e-10, 1.84e-10).
- Single tile PG L24 2.6e-11, L32 2.7e-11; node PG (32 ranks, CCS 4) max 6.5e-10 PASS (V4 5.8e-10).
- Planted: tau 1e-12: same seed FAIL (4.2e-13 / 2.0e-12), node PG 2.0e-8 FAIL; wrong conj: FAIL (run aborts). Races: single tau slot / single exchange slot / single v buffer NOT caught (PG identical to clean);
  apply chunk-end barrier removed and single staging buffer: caught in the node config (abort "Unstable matrix").
- Timing alternating V4 vs V7 (async+staging; devqr1 vs devqr2+devunmqr): L24 5.85-5.91 -> 4.95-4.97 s/sweep (x3); L32 18.70/18.44 -> 16.41/16.42 (x2).
  Node CPU-time-diff, sweeps/s per card: 32 ranks 4/tile CCS 4 OMP 3: V4 0.812/0.815, V7 0.841/0.846; 8 ranks 1/tile OMP 14: V4 0.315/0.302, V7 0.354/0.361.
