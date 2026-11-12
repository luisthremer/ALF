# Intel Max 1550 vs NVIDIA A30: ALF finite-T DQMC, point-defect Hubbard model

Team session 2026-11-11/12.
- Branches: `production-mpi-intel-large-prompt-team` (Phases 1-2, staging buffer); `-team-async` (Phase 3/3b); `-team-sycqr` (Phase 4); and `-team-combined` (all of them; worktree `ALF_wt_lp_team_combined`).
- "C" = confirmed by our own job on SuperMUC-NG Phase 2; "E" = estimate.
- Scripts and raw logs: `phase2/team/` (remote `$SCRATCH/large_prompt/team/`). WBS: `phase2/team/wbs_team.json`.

## 1. Problem

On paper a Max 1550 is several times an A30, yet one ALF chain ran ~2x slower on it.
- Reference parameters `ref_L24/`: L=24, N=576, N_FL=2, Beta=24, Dtau=0.1, U=5.3, submatrix N_delay=32, Nwrap=10.
- A30 targets (TUNING_DEFECT.md, REPORT_2026-11-10.md, not re-measured):
  - one chain: 5.25 s/sweep (L=24), 18.2 s/sweep (L=32);
  - 0.37 sweeps/s per GPU with 4 chains, 0.43 with 8.

## 2. Before / after against the target (combined branch, one binary, run-time switches)

| metric | A30 | Max 1550 before | after: staging + async + device QR (V4) | status |
|---|---|---|---|---|
| L=24, one chain, one tile [s/sweep] | 5.25 | 10.5-10.7 (C) | 6.12-6.17 (C) | 1.17x the A30, not met |
| L=32, one chain, one tile [s/sweep] | 18.2 | 30.3-31.7 (C) | 19.1-19.6 (C) | 1.06x the A30, nearly met |
| sweeps/s per card, L=24 (32 ranks/node, 4 per tile, CCS 4) | 0.43 | 0.549 (C) | 0.737-0.785 (C) | **met, 1.7-1.8x the A30** |

Steps, each a run-time switch of the same binary (jobs 560386-88, 560581-83; alternating runs, in-loop timers):

| variant | switches | L=24 | L=32 | per card | numerics |
|---|---|---|---|---|---|
| V0 | none | 10.5-10.7 | 30.3-31.7 | 0.549 | reference |
| V1 | `LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576` | 7.05-7.08 | 23.2 | 0.568 | bitwise identical |
| V2 | V1 + `ALFSH_ASYNCBLK=1` | 6.53-6.59 | 22.0-22.4 | 0.655-0.688 | bitwise identical |
| V4 | V2 + `ALFSH_DEVQR=1` | 6.12-6.18 | 19.1-19.6 | 0.737-0.785 | rounding level, all checks pass (sec. 5) |

L=24 V4 phase split (C): update 1.44-1.47, stabilisation 1.40, CGR 1.83-1.86, observables 1.19 s/sweep.

## 3. Profile and cause (before the changes)

| phase at L=24 | s/sweep (C) | where it runs |
|---|---|---|
| updates (submatrix) | 5.1-5.6 (~55 %) | device, latency bound |
| CGR (incl. QR/LU) | 1.74 | device + host (oneMKL hybrid) |
| stabilisation (QR) | 1.39 | device + host (oneMKL hybrid) |
| observables | 1.18 | device |

- **Update phase.** Each sweep made ~75,000 synchronous kernel launches, ~40,000 dispatched BLAS calls and ~374,000 small H->D copies (C, LIBOMPTARGET_DEBUG counts).
  - Its single biggest cost: the 32 KB D->H copy of the gathered block, 8640 per sweep. It took 380 us with the default staging buffer and 22 us with a 1 MB staging buffer (C). That alone is -3.5 s/sweep at L=24.
- **oneMKL geqrf/getrf.** They are hybrid. The geqrf panel (`host::geqrf_wy2`) and `larfb_wy` run on the host CPU (VTune call stacks, C). Before, this cost 2.4 s/sweep at L=24 and 9.8 s at L=32.

## 4. What changed

1. **Phase 1: oneMKL batched geqrf/getrf.** No gain. Reverted (C).
2. **Phase 2 (de3f7306): explicit-shape dummies and a fused kernel in the submatrix update.**
   - First-private arrays dropped from 168,800 to 28,544 per sweep, launches from ~75,000 to 49,530.
   - L=24 -0.4 to -1.2 s, bitwise identical (C).
3. **Phase 3: OpenMP `nowait`/`depend`.** Gate A failed (97-103 %): ifx runs these target tasks undeferred (C). Dropped.
4. **Phase 3b: staging buffer, 1 MB** (environment variable, in the job files, b2343b98).
   - The largest single lever: L=24 10.6 -> 7.1 s (C).
5. **Phase 3b: SYCL-async block path (`ALFSH_ASYNCBLK=1`).**
   - Block kernels are submitted on the OpenMP interop in-order queue without a wait: 4.3 us per op vs 12.3 us synchronous (C).
   - L=24 -0.5 s, L=32 -1 s, per card +15-20 % (C).
6. **Phase 4: own device Householder QR (`ALFSH_DEVQR=1`, `libalfdevqr.so`, dlopen'ed).**
   - It replaces oneMKL's hybrid geqrf: a register-resident panel kernel (ZGEQR2+ZLARFT in one work-group), then a blocked update.
   - Falls back to oneMKL for M > 1536 or when the device limits are not met.
   - L=24 -0.4 s, L=32 -2.6 s, per card +10-15 % (C).

## 5. Numerics (all C)

- Same seed vs the same-compiler ifx CPU build with V4 (and V3):
  - Mz=F: L6 2.2e-15, L12 1.6e-14 (tol 1e-13);
  - Mz=T: 6.0e-11 (tol 1e-10).
- V1 and V2 are bitwise identical to the branch head.
- Precision Green:
  - single tile: L=24 1.8e-14 / 4.5e-11, L=32 5.1e-14 / 4.4e-11;
  - node configuration (32 ranks, CCS 4): V2 max 3.7e-10, V4 max 5.8e-10.
- **Device-QR race, found and fixed.** With 4 ranks per tile under CCS, the first device QR gave wrong results (Precision Green max 0.55-1.96).
  - Cause: the tau buffer of ZLARFT was reused two panel steps early.
  - Fixed in 29005a0f: 3 slots, plus a device-limit check.
  - The pre-fix library was kept as a planted defect. In the node configuration it FAILS the check (max 0.40), and the fixed one passes (job 560582). The race only shows under concurrent ranks per tile, so this node check is now part of the device-QR gates.
- Other controls:
  - A planted 1e-12 G perturbation is caught by the CPU comparison only marginally (1.1e-13 vs tol 1e-13).
  - `ALF_GPU_NS=0` is NOT caught at L=12. The bitwise GPU-vs-GPU check catches a 1e-12 factor perturbation.
- CUDA path: compiles with nvfortran on the A30 host. `map(present,...)` is Intel-only.
- Open: a fresh numerics review of the device QR before production use (WBS).

## 6. Tools on this system (C)

- What does not work:
  - VTune 2025.1/2025.8 `gpu-offload` collects no GPU data;
  - `LIBOMPTARGET_PROFILE` prints nothing;
  - PTI segfaults;
  - standalone SYCL programs need `SYCL_UR_USE_LEVEL_ZERO_V2=1` and then time differently from ALF.
- What works: shim timers (`ALFSH_TIMING=1`; with the async path they include queued work), phase timers (`-DALF_TEAM_TIMERS`), `LIBOMPTARGET_DEBUG=1` counts, and VTune host hotspots.

## 7. AMD comparison and next steps (Phase 5)

- **The AMD MI210 session** (branch `gpu-delayed-submatrix-defect-production-mpi-amd`) reached L=24 4.41 s/sweep and L=32 17.4 (their measurement). Their levers:
  - their own HIP Householder QR for geqrf AND unmqr (wave-per-column reflector application, no T factor);
  - syncs deferred to the next kernel launch (patched libgomp plugin, ROCm-specific).
  - Their observables take 0.12 s of device time at L=24, against 1.19 s phase time here.
- **Phase 5 plan** (WBS), to close the remaining ~0.9 s at L=24:
  - a clean per-call profile;
  - device ZUNMQR (AMD-style apply);
  - device ZGETRF;
  - the observables phase (dummy-argument treatment and fusion);
  - then the same checks as above.
- **Production recommendation now:** the staging buffer (V1) is free and bitwise identical; V2 is also bitwise identical. Keep V4 off in production until the fresh numerics review is done.

## 8. Side finding: faulty GPU tile in production (2026-11-12)

- Sporadic production aborts ("G difference > 10") were traced to tile 1 of node i20r02c06s02 by the production session: the same input aborts there and runs clean with the tile hidden or on another node.
- Our stress test (`diag/t_stress.F90`) repeats the dev_cgr chain at N=1024 with a bitwise check every iteration. Its positive control catches a planted error on all 8 tiles.
- LRZ ticket draft: `diag/LRZ_TICKET_i20r02c06s02.md`.

## 9. Budget

- GPU: about 3.0 node-h for Phases 1-2, plus about 10 node-h for Phases 3/3b/4, the combined evaluation, and the gates (E; sacct unavailable, from job logs). Partition `preempt` (owner OK).
- Usage: weekly 30 % at the start -> 41 % (2026-11-12 ~10:30), cap 25 %.
