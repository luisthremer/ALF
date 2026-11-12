# AMD (MI210, gfx90a) notes

Toolchain: gfortran 14 + GCC amdgcn offload, hipBLAS/hipSOLVER via `alf_la_hip_*` (see configure.sh `GPU=AMD`).

## Faster libgomp GCN plugin (optional, ~3x on the real run)
`plugin-gcn-pool-prof.patch` applies to GCC 14.2.0 `libgomp/plugin/plugin-gcn.c` (after changing the wait state at
`run_kernel` to `HSA_WAIT_STATE_ACTIVE`, which alone gave nothing):
- pools the per-launch kernarg block and completion signal (empty target region 253 us -> 37 us),
- pools small device allocations,
- `GCN_PROF=1` prints per-kernel wall time at exit.
Build only `libgomp-plugin-gcn.la`, put its `.libs` first in `LD_LIBRARY_PATH`.

## Runtime settings
`GCN_NUM_THREADS=4` (default 16): parallel-region/reduction startup on GCN scales with thread count; 4 was best at L=16.

## Measured (L=16, 2 sweeps, 1 process): 66 s CPU stock -> 21 s (plugin + threads + obs kernels)

## HIP kernels (2026-11-11)
`Libraries/Modules/alf_hip_kernels.hip.cpp` (hipcc) + `alf_hipk_mod.F90` replace the hot small OpenMP-target kernels
(copies, block gather/scale, flush gather, scatter, correlations, scalar/local observables, stabilisation helpers, slice
multiplication) under `#if defined(ALF_LA_HIP)`; the OpenMP versions remain for other backends. A HIP launch+sync is ~15 us,
a libgomp-GCN parallel target region ~150 us. Remaining GPU time is dominated by rocSOLVER's many tiny QR kernels.

Single process, 2 sweeps, CPU time: L=16 66 s -> 6 s; L=24 95 s -> 19 s (A30 reference ~10 s/sweep with 4 ranks).
4 concurrent L=24 processes: 42 s wall for 2 sweeps each. L=28 ~16.6 s/sweep, L=32 ~28.8 s/sweep (single process).

## QR in HIP, deferred synchronisation (2026-11-11, later)
Profile of the L=24 sweep (rocprofv3 kernel trace): 2.07M kernel dispatches per 2 sweeps, 93% of them from rocSOLVER's
`zgeqrf` (2392 kernels, 8.8 ms per call at n=576) and `zunmqr` (2138 kernels, 8.5 ms); together 44% of the sweep.
- `Libraries/Modules/alf_hip_qr.hip.cpp`: blocked Householder QR (zlarfg/LAPACK output format) with a register-resident
  panel kernel (NB=16, 256 threads; NB=8, 512 threads for m > 1024) and a wave-per-columns reflector-application kernel used
  for the trailing update and for the whole `unmqr` (no T factor, no GEMM). m <= 2304, else (or `ALF_HIP_QR=0`) hipSOLVER.
  n=576: geqrf 8.8 -> 4.4 ms, unmqr 8.5 -> 1.7 ms. Standalone test `qr_test.cpp` (vs hipSOLVER; planted defects are caught).
- Deferred synchronisation: HIP wrappers no longer `hipDeviceSynchronize` after every launch; the patched libgomp plugin
  (`plugin-gcn-presync-hook.patch`, applied on top of the pool patch) calls `alf_hip_presync()` before every OpenMP kernel
  launch and device copy/free. Needs that plugin (else every call syncs); `ALF_HIP_ASYNC=0` forces the old behaviour.
  Output is bit-identical to the synchronous mode.
- All remaining stabilisation kernels (ident, permutations, ct, colinv, phases, ...), colnorms and permute_cols are HIP now.

Steady-state seconds per sweep (difference of 2- and 4-sweep runs, MI210, same parameters as the A30 production run):

| | A30 (solo, N_delay=32) | A30 4 ranks + MPS | MI210 solo | MI210 4 concurrent (no MPS) |
|---|---|---|---|---|
| L=24 | 5.25 | 10.5 per rank | 4.41 | 8.9 per rank |
| L=28 | - | - | 8.84 | 18.8 per rank |
| L=32 | 18.7 | - | 17.5 | - |
