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

| metric | A30 | Max 1550 before | V2 (bitwise identical) | V4 | V7 (Phase 5) | status |
|---|---|---|---|---|---|---|
| L=24, one chain, one tile [s/sweep] | 5.25 | 10.5-10.7 (C) | 6.53-6.59 (C) | 5.85-6.17 (C) | **4.95-4.97 (C)** | met with V7 |
| L=32, one chain, one tile [s/sweep] | 18.2 | 30.3-31.7 (C) | 22.0-22.4 (C) | 18.4-19.6 (C) | **16.41-16.42 (C)** | met with V7 |
| sweeps/s per card, L=24 (32 ranks/node, 4 per tile, CCS 4) | 0.43 | 0.549 (C) | 0.655-0.688 (C) | 0.74-0.82 (C) | **0.84-0.88 (C)** | **met, ~2x the A30** |

V4 and V7 were measured in different jobs; the V4 vs V7 alternation in one job gave L24 5.85-5.91 -> 4.95-4.97, L32 18.44-18.70 -> 16.41-16.42, per card 0.812-0.815 -> 0.841-0.846 (C).

Steps, each a run-time switch of the same binary (jobs 560386-88, 560581-83; alternating runs, in-loop timers):

| variant | switches | L=24 | L=32 | per card | numerics |
|---|---|---|---|---|---|
| V0 | none | 10.5-10.7 | 30.3-31.7 | 0.549 | reference |
| V1 | `LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576` | 7.05-7.08 | 23.2 | 0.568 | bitwise identical |
| V2 | V1 + `ALFSH_ASYNCBLK=1` | 6.53-6.59 | 22.0-22.4 | 0.655-0.688 | bitwise identical |
| V4 | V2 + `ALFSH_DEVQR=1` | 6.12-6.18 | 19.1-19.6 | 0.737-0.785 | rounding level, all checks pass (sec. 5) |
| V7 | V2 + `ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1` | 4.95-4.97 | 16.41-16.42 | 0.841-0.879 | rounding level; L6T open, race test partly blind (sec. 7) |

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

## 7. Phase 5: second device QR and device ZUNMQR (V7)

- **AMD MI210 session** (branch `gpu-delayed-submatrix-defect-production-mpi-amd`): L=24 4.41 s/sweep, L=32 17.4 (their measurement), with their own HIP Householder QR for geqrf AND unmqr (wave-per-column apply, no T factor) and deferred syncs (patched libgomp, ROCm-specific). Phase 5 took the apply idea over.
- **Clean profile before Phase 5** (job 560667, sync variants, marginal 4 vs 2 sweeps, per sweep): L24 geqrf device 1.50 s (oneMKL 2.09), unmqr 0.50, getrf 0.41; L32 geqrf 5.20 (8.17), unmqr 1.95, getrf 1.24 (C).
- **Cause in the Phase 4 panel** (C, microbenchmarks): ~7 us per column, dominated by the synchronisation chain and the ZLARFT T recurrence; SLM bandwidth limits the kernels.
- **What was built** (`alf_la_sycl_qr2.hpp`, default off): a T-free ZGEQR2 panel in one work-group, a reflector-apply kernel used for the QR trailing update and for ZUNMQR (`ALFSH_DEVQR=2`, `ALFSH_DEVUNMQR=1`, fallback to oneMKL when the device limits are not met).
- **Per call inside ALF** (C, job 560722): L24 geqrf 5.48 -> 3.91 ms, unmqr 2.4-2.6 -> 1.37 ms; L32 geqrf 16.95 -> 10.36 ms, unmqr 7.4-7.6 -> 5.30 ms. Gate B (geqrf <= 0.6x, unmqr <= 0.5x) missed narrowly (0.61-0.71x and 0.52-0.71x), but the sweep targets are met.
- **Checks (C):**
  - Standalone vs host MKL: R, V, tau <= 7e-14; unmqr <= 2.3e-15; planted tau / conj / scale defects all FAIL.
  - Default build (no switches) bitwise identical to before at L6F/L12F.
  - Same seed vs CPU: L6F 2.3e-15, L12F 1.6e-14 PASS; **L6T 1.36e-10 FAIL (tol 1e-10)**. The criterion is noise-limited: oneMKL itself fails 2 of 5 seed sets (1.13e-10, 1.84e-10). A new pre-registered L6T rule needs the owner's approval.
  - Precision Green: single tile L24 2.6e-11, L32 2.7e-11; node configuration (32 ranks, CCS 4) max 6.5e-10 PASS. Planted tau 1e-12 fails both same seed and node PG.
  - Three narrow planted races (single tau slot, single exchange slot, single v buffer) are NOT caught by node PG; the gross ones are.
- **Race test** (`p5/race/NOTES.md`, jobs 560754-560807, ~5.4 node-h, C): stress test (`diag/t_stress.F90`, bitwise vs iteration 1) with V7 kernels, 32 ranks, 4 per tile, CCS 4.
  - Real kernels: **0 bad results** in 0.85 M iterations (N=1024) and 3.39 M (N=576) on two nodes; oneMKL and DEVQR=1 controls also 0.
  - Plants: one-buffer caught at once; exchange slot caught rarely (4 in 337k); v buffer caught only with an amplifier setting (`ALFDEVQR_G=1`); **tau slot never caught** (2.2 M iterations).
  - Meaning: races at >= 1e-5 per iteration are excluded (95 %: < 3.5e-6 at N=1024, < 9e-7 at N=576); a tau-type race (the Phase 4 bug class) would be missed. Code review of the tau-slot logic not done.
- **Throughput vs walkers per tile** (job 560772, V7, L24, sweeps/s per card, C): 1/tile 0.365, 2/tile 0.458, **4/tile 0.879**, 8/tile 0.261. V7 helps one chain much more than throughput (+3.7 % at 4/tile, +16 % at 1/tile).
- **Not done:** device ZGETRF, the observables phase.

## 8. Phase 6: batched multi-walker rebuild (assessed, no-go recommended)

- **Work per walker-sweep at L24** (job 560773, MKL_VERBOSE, C for BLAS): 8.8 TFLOP in oneMKL (5376 large zgemm 8.2 TF, 34,560 rank-k updates 0.49 TF, ztrsm 0.07) + ~0.4 TF device QR and own kernels (E) = ~9.2 TFLOP. A tile at 4 walkers runs ~4 TFLOP/s.
- **Ceiling benchmark** (`p6/ceil_bench.F90`, jobs 560780-560801, C): zgemm 576^3 12.9 TFLOP/s alone, ~11-12 under node load; batched up to 16.7 alone, 12.4 loaded. The rank-k updates are memory-bandwidth bound under load; QR, LU and observables are latency bound.
- **Ceiling (E):** formula on gemm speed 3.1-4.1x per card; grounded 1.5-2.3x; realistic 1.3-1.6x.
- **Design** (`p6/DESIGN.md`): walker index on all device arrays, lockstep walkers, per-walker RNG; ~62 person-days (E).
- **Prototype** (batched hopping wrap, kept as `p6/proto_walker_batch.patch`, not in the source):
  - oneMKL batched/wide zgemm is NOT bitwise equal to the single call at N=144/400/576 (rounding-level differences that depend on the batch size, C), so the bitwise-per-walker gate fails; the bit-safe mode gives no gain.
  - Batching the wrap saves ~2 % of a sweep (C).
- **Phase split, one chain L24 V7** (C): updates 1.64, stabilisation 0.90, CGR 1.23, observables 1.31 s/sweep; under CCS 4, 70-90 % of each phase is latency (E, fit).
- **Recommendation:** no full rebuild. If more throughput is wanted: try launch-latency settings on the existing binary first (free, untested), then a bit-safe pilot with the walker index only in our own kernels (updates, observables; 15-20 person-days, E).

## 9. Production recommendation (2026-11-12)

- Production runs branch `production-mpi-intel-unstable-guard` (same base 5667acda; it differs only in `Prog/control_mod.F90`; a test merge with this branch has no conflicts).
- **V1 now:** `LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576` in the job script. Environment only, no rebuild; bitwise identical in our tests. Confirm with one short same-seed run of the production binary with and without it (bitwise).
- **V2 next:** merge this branch into the production branch, rebuild, set `ALFSH_ASYNCBLK=1`; confirm bitwise identity against the current production build (same seed, short run) before switching. Gain about +15-20 % per card over V1.
- **Not yet V4/V7:** V4 waits for the final numerics review; V7 for the L6T rule and the tau-slot review. Their extra throughput at 4 walkers per tile is +7-20 % over V2.

## 10. Side finding: faulty GPU tile in production (2026-11-12)

- Sporadic production aborts ("G difference > 10") were traced to tile 1 of node i20r02c06s02 by the production session: the same input aborts there and runs clean with the tile hidden or on another node.
- Our stress test (`diag/t_stress.F90`) repeats the dev_cgr chain at N=1024 with a bitwise check every iteration. Its positive control catches a planted error on all 8 tiles.
- LRZ ticket sent by the owner (2026-11-12 ~12:20); closed from our side. Control node 560639: 100,711 iterations on 8 tiles, 0 differences (C). The suspect-node stress runs were cancelled before running (owner).

## 11. Budget

- GPU: about 51 node-h in total (E; sacct unavailable, from job logs), of the 250 node-h (1000 GPU-h) team limit. Of this: Phases 1-4 and gates ~13, faulty-node diagnosis and stress tests, Phase 5, race test ~5.4, Phase 6 ~2.5. Partition `preempt` (owner OK).
- Usage: weekly 30 % at the start -> 41 % (2026-11-12 ~10:30), cap 25 %.
