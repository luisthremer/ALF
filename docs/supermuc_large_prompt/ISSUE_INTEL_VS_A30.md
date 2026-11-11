# Intel Max 1550 vs NVIDIA A30: ALF finite-T DQMC, point-defect Hubbard model

Team session 2026-11-11, branch `production-mpi-intel-large-prompt-team` (worktree `ALF_wt_large_prompt_team`).
"C" = confirmed by our own job on SuperMUC-NG Phase 2; "E" = estimate. Scripts and raw logs: `phase2/team/`
(remote `$SCRATCH/large_prompt/team/`). WBS: `phase2/team/wbs_team.json`.

## 1. Problem

On paper a Max 1550 is several times an A30, yet one ALF chain runs slower on it. Reference parameters `ref_L24/`
(L=24, N=576, N_FL=2, Beta=24, Dtau=0.1, U=5.3, submatrix N_delay=32, Nwrap=10). A30 targets (TUNING_DEFECT.md,
REPORT_2026-11-10.md, not re-measured): one chain 5.25 s/sweep (L=24), 18.2 s/sweep (L=32); 0.37 sweeps/s per GPU
with 4 chains, 0.43 with 8.

## 2. Before / after against the target

| metric | A30 | Max 1550 before | Max 1550 after | status |
|---|---|---|---|---|
| L=24, one chain, one tile [s/sweep] | 5.25 | 9.7 (9.0-10.5; 10.1-14.6 in other jobs) (C) | 9.3 (8.5-9.9) (C) | not met, ~1.8x |
| L=24, one chain, one card (COMPOSITE) [s/sweep] | 5.25 | 17.3-19.4 (C) | not re-measured | worse than one tile |
| L=32, one chain, one tile [s/sweep] | 18.2 | 29-32 (C) | 27.4 (1 run) (C) | not met, ~1.5x |
| sweeps/s per card, L=24 (32 ranks/node = 4 per tile, ZEX_NUMBER_OF_CCS) | 0.43 | 0.51 (C) | 0.52 (C) | **met** |

The per-card throughput target was already met before any change (it contradicts the premise "slower"); only the
single-chain speed is behind. Before/after timings are from alternating runs in the same job with the in-loop sweep
timer (no setup); run-to-run spread of the L=24 sweep is +-0.7 s, almost all of it in the update phase.

## 3. Profile (L=24, one tile, 8 host threads, per sweep, before the change)

| phase | s/sweep (C) | share | where it runs |
|---|---|---|---|
| updates (`wrapgrup`+`wrapgrdo`, submatrix) | 5.1-5.6 | ~55 % | device, latency bound |
| CGR (incl. QR/LU) | 1.74 | 17 % | device + host (oneMKL hybrid) |
| stabilisation (`stab_wrapur/ul`, QR) | 1.39 | 14 % | device + host (oneMKL hybrid) |
| observables (device equal-time) | 1.18 | 12 % | device |

- oneMKL LAPACK inside CGR + stabilisation: 2.43 s/sweep (C): 192 geqrf (9.4 ms/call with 8 host threads, 27 ms with
  1), 192 unmqr (2.8 ms), 96 getrf (3.9 ms). L=32: 9.8 s of 31.4 s (geqrf 29 ms/call).
- What runs on the host and why: oneMKL's device `geqrf` takes `geqrf_wy_lookahead_hybrid`; the panel (`host::geqrf_wy2`:
  MKL zgeqrf, zlarft, ztrmm) and `host::larfb_wy` (host zgemm) run on the CPU (VTune call stacks, C). The oneMKL reference
  documents geqrf on GPU as hybrid; there is no switch. getrf is host-assisted as well. Everything else is on the device
  (host-side `la_qrpp` sort and `Pivot_phase` < 0.1 %).
- Why the update phase is slow: per sweep ~75,000 synchronous kernel launches (60,000 in `gpu_lowrank_mod`), ~40,000
  dispatched oneMKL BLAS calls (each with an interop sync), ~374,000 small H->D copies and ~168,800 "first-private array"
  allocations (descriptors of module allocatables) (C, LIBOMPTARGET_DEBUG counts). Measured unit costs: empty sync kernel
  6 us, small kernel 9 us, small dispatched zgemm 18 us, 16-byte copy 0.8 us (C); with immediate command lists off the
  non-LAPACK part grows by 4.8 s/sweep (C). Launch/sync overhead is ~1.7-2.5 s/sweep (E); the rest of the update time is
  device time that could not be measured (no working GPU tracer, see 6).
- Device utilisation of one chain is low: 4 chains per tile give 2.6x the throughput of 1 chain (C), so one chain uses
  well under half of a tile.

## 4. Cause

Two parts: (a) the update phase is bound by per-launch/per-sync latency, which is higher on PVC/Level Zero with ifx
OpenMP offload than on the A30 (where the same code is also latency bound, per TUNING_DEFECT.md); (b) oneMKL's QR/LU on
the GPU is hybrid and does its panels on the host (2.4 s/sweep at L=24, 9.8 s at L=32; on the A30 QR/LU is ~9 % of kernel
time).

## 5. What was tried and what changed

| step | result |
|---|---|
| Phase 1: oneMKL batched `geqrf_batch`/`getrf_batch` (batch 1) via the shim | no gain: geqrf 9.7/28.9 ms vs 10.1/28.4 ms (same path, same Precision Green); getrf 2.6-4x slower (C). Reverted, patch `phase2/team/phase1_batch_switch.patch` |
| Phase 2 (commit de3f7306): explicit-shape dummies instead of module allocatables in the submatrix kernels; fused `k_blk_scale` | bitwise-identical output; first-private arrays 168,800 -> 28,544, launches ~75,000 -> 49,530 per sweep (C); L=24 -0.4 s (4 alternating pairs) to -1.2 s (prototyper, 3 pairs), L=32 -3.8 s (1 pair) (C) |
| runtime knobs | immediate command lists already default (off: +4.8 s); 4-16 host threads best for the hybrid geqrf (C) |

Numerics (C): same-seed vs same-compiler ifx CPU build: Mz=F 2.2e-15 (L6) / 1.5e-14 (L12), Mz=T 5.5e-11, configurations
and acceptance identical, ctest 53/53, Precision Green unchanged (L24 mean 1.46e-14, max 3.07e-11). Controls: planted G
perturbation 1e-12 is caught by the CPU comparison only marginally (1.1e-13 vs tol 1e-13); `ALF_GPU_NS=0` is NOT caught at
L=12 (8e-15) - the CPU comparison is blind to small orthogonality losses of Q; the new GPU-vs-old-GPU bitwise check
catches a 1e-12 factor perturbation (1.9e-11). CUDA path: compiles with nvfortran (A30 host), not run (GPU busy with the
production job); `map(present,...)` is Intel-only because nvfortran 26.9 rejects it.

## 6. Tools that do not work on this system (C)

VTune 2025.1 and 2025.8 `gpu-offload` collect no GPU data; `LIBOMPTARGET_PROFILE` prints nothing; the PTI view tracer
segfaults with the Level Zero tracing layer and records nothing without it; standalone SYCL programs need
`SYCL_UR_USE_LEVEL_ZERO_V2=1` and then time oneMKL differently from ALF (not representative). What works: shim timers
(`ALFSH_TIMING=1`), phase timers (`-DALF_TEAM_TIMERS`), `LIBOMPTARGET_DEBUG=1` counts, VTune host hotspots.

## 7. Remaining options (owner decisions)

| option | gain (E) L=24 / L=32 | effort | risk |
|---|---|---|---|
| B3: asynchronous update phase (nowait + dependencies, Intel path), remove most per-kernel syncs | 1.5-2.5 s / 3-5 s | days | ordering bugs; needs the bitwise GPU-vs-GPU check |
| A3: own SYCL blocked Householder QR (LAPACK algorithm, fused panel kernel) in the shim | ~1 s / ~4.7 s | days | rounding-level change vs oneMKL; same-seed check may move |
| B1 for the rest: same dummy-argument treatment for `gpu_obs_mod` and `gpu_wrap_mod` kernels (remaining 28,544 first-private arrays) | ~0.1-0.2 s | hours | low |
| run parameter (owner): `Nwrap`, `N_delay` | earlier session: -36 % with Nwrap 20, N_delay 64 (not re-checked) | none | moves numerics; proposal only |

Even with B3 + A3 the single L=24 chain is estimated at ~6 s (E), still above 5.25 s; the A30 has lower per-launch cost
and no host-side QR. For production throughput the Max 1550 is already ahead with 4 ranks per tile.

## 8. Budget

GPU: ~3.0 node-h ~ 12 GPU-h (E; sacct unavailable, from job logs), partition `preempt` (owner OK), one job preempted
(560275, partial data used). Usage: weekly 30 % -> ~36 % (~6 % of the 10 % cap); 2 Sonnet subagents (prototyper,
numerics reviewer).
