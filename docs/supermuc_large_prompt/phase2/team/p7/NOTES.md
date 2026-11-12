# Phase 7: kernel-level GPU measurement (where does the GPU time go?)

2026-11-12, SuperMUC-NG Phase 2, ifx/oneMKL 2025.3.0, binary `$SCRATCH/large_prompt/team/p5/build_aottt/Prog/ALF.out` (unchanged), reference input `ref_L24` (L=24, N=576, N_FL=2, submatrix N_delay=32, Nwrap=10).
"C" = measured by a job of this session (job id given), "E" = estimate. Measure only: no source under Prog/ or Libraries/ was edited; the numbers come from run-time switches, the existing build, and
four small LD_PRELOAD measurement libraries (`src/`, built into `$SCRATCH/large_prompt/team/p7/lib/`, not part of the source tree).
Nodes: 1 node-job each; all runs `-p preempt --qos=preempt -x i20r02c06s02`. Layouts: **L1** = one rank on tile 0 (OMP 8); **L32** = 32 ranks, 4 per tile, CCS 4, OMP 3.
V7 = `LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 ALFSH_ASYNCBLK=1 ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1`, V2 = same without the two DEVQR switches.

## 1. Tools (what worked)

| tool | result |
|---|---|
| `unitrace` / `onetrace` | **not installed**: no binary in any oneAPI tree under /lrz/sys/intel (2025.0.1 ... 2026.2.0), no module. Only the PTI *library* 0.11 (oneapi_2025.0.1 / 2025.1.0 / 2025.2.1). Not built from source (no authority to download; a local clone of pti-gpu made during the first minutes was deleted unused). |
| any Level Zero tracing-layer tool (`ZE_ENABLE_TRACING_LAYER=1`: PTI, unitrace, onetrace) | **crashes ALF at start-up** (SIGSEGV, jobs 560821 P3, 560825 gdb). Cause (C, gdb backtrace): with the layer on, the OpenMP plugin's SYCL interop init (`__tgt_sycl_init_interop` -> `make_platform<level_zero>`) throws `UR_RESULT_ERROR_UNINITIALIZED` and the C++ unwinder jumps to a stack address. A trivial OpenMP-offload program runs fine with the layer on (560823). `SYCL_UR_USE_LEVEL_ZERO_V2=1` makes ALF abort (SIGABRT, 560826/560829). |
| Intel Advisor 2025.4 `--collect=survey --profile-gpu` | **failed** (560823/560825/560826): "Roofline collection is not supported for MPI run"; without mpiexec the app returned non-zero and the data file is "corrupted". No GPU data. |
| VTune | not retried (known not to work). |
| XVE_ACTIVE/XVE_STALL, occupancy metrics (metric query) | **not available**: no unitrace, and `/proc/sys/dev/i915/perf_stream_paranoid` is not readable on the compute node. Occupancy below is therefore **launch-geometry based (E)**, not hardware counters. |
| `LIBOMPTARGET_PLUGIN_PROFILE=T` | **works** (note: different from `LIBOMPTARGET_PROFILE`): per OpenMP kernel count, host time and *device time from Level Zero event timestamps*, plus DataAlloc/DataRead/DataWrite count and time. Printed per rank at exit (`Thread 0`). Covers only kernels launched by libomptarget. |
| `LIBOMPTARGET_INFO=-1` | **works**: per launch "Number of teams / Team sizes" (launch size). Raw log is 3.6 M lines per sweep, parsed on the node (`jobs/parse_info.py`, `data/kinfo_V7.txt`). |
| `MKL_VERBOSE=1` | counts only: the SYCL-API lines carry `host:N/A device:N/A` (no time). |
| own LD_PRELOAD libraries (`src/`) | **work**, bitwise neutral (Precision Green identical, 560845 pass-through run). Method: after a wait on the application queue, run the same call on a separate in-order *profiling* queue of the same context/device, wait, read `command_end - command_start` of the event. Fully serialised, so device times are per-kernel (isolation) times. `alf_la_sycl_blk_p7.cpp` = block path (gather/scale/fgather/copy/rank-k gemm), `p7off.cpp` = the oneMKL OpenMP-dispatch variants `mkl_blas_{zgemm,ztrsm,zgemm_batch_strided}_omp_offload_lp64` (with call-site return addresses, resolved with `nm`), `p7mkl.cpp` = exported SYCL gemm/trsm (superseded; dispatched calls do not pass through it). `ALFDEVQR_PROF=1` (existing switch) gives the device-QR panel/apply split. |
| the 32-rank serialised run (job 560846) | **hung** (ranks idle after ~200 QR calls), cancelled after 7 min. So L32 has only the OpenMP-kernel profile (C) + phase timers (C); SYCL/BLAS device times at L32 are not measured. |

## 2. Per-kernel table, V7, L1 (one chain), per sweep

Sources: omp kernels = PLUGIN_PROFILE 3 sweeps minus 1 sweep, /2 (job 560831); block ops and dispatched BLAS = serialised profile (560845, 560847); QR/LU = `ALFSH_TIMING` marginal (560831). Launch size: omp = teams x team-size (INFO, 560831), SYCL = global range from the source.
Occupancy proxy (E): Max 1550 tile = 64 Xe cores x 8 XVE x 8 threads = 4096 hardware thread slots (SIMD32 for omp, SIMD16 for the QR kernels).

Top 20 by device time (`data/report_tables.md` has the same, plus per phase):

| # | phase | kernel / op | calls/sweep | device s/sweep | mean us | launch size | fill of tile (E) |
|--|--|--|--|--|--|--|--|
| 1 | S+C | device QR geqrf call (sync; stages: panel 60 %, apply 40 %) | 192 | 0.822 | 4281 | panel: 1 work-group x 512 WI; apply: ~70 WG x 256 WI | panel 0.8 %, apply ~27 % |
| 2 | C | getrf oneMKL call (hybrid host/device, sync) | 96 | 0.379 | 3943 | oneMKL | partly on the host |
| 3 | U | zgemm 576^3 (hopping wraps, dispatch) | 2880 | 0.315 | 109 | oneMKL | full (14 TFLOP/s) |
| 4 | U | block rank-k gemm 576x576xk (SYCL) | 34560 | 0.312 | 9.0 | oneMKL, k=10-14 | memory bound (~1.2 TB/s) |
| 5 | O | corr l437 (omp) | 480 | 0.300 | 626 | **1 x 1** | 0.02 % |
| 6 | S+C | device QR unmqr call | 192 | 0.262 | 1367 | ~72 WG x 256 WI | ~27 % |
| 7 | O | zgemm 576^3 (symm_apply) | 1920 | 0.206 | 107 | oneMKL | full |
| 8 | O | corr l473 | 480 | 0.166 | 346 | **1 x 1** | 0.02 % |
| 9 | O | corr l450 | 480 | 0.153 | 319 | 512 x 1024 | full |
| 10 | U | block scale (SYCL) | 8640 | 0.093 | 10.8 | 2xbx576, mean 36864 WI | ~28 % |
| 11 | O | corr l418 | 480 | 0.090 | 187 | 512 x 1024 | full |
| 12 | U | block fgather (SYCL) | 17280 | 0.084 | 4.8 | kx576, mean 6045 WI | ~5 % |
| 13 | O | corr l386 | 480 | 0.068 | 142 | 512 x 1024 | full |
| 14 | C | k_cgr_phase l1232 | 96 | 0.059 | 610 | 1 x 64 | 0.05 % |
| 15 | C | k_prec l1327 | 96 | 0.051 | 529 | 4 x 1024 | ~3 % |
| 16 | S+C | k_colnorms l845 | 192 | 0.044 | 231 | 288 x 2 | ~1 % |
| 17 | U | block gather (SYCL) | 8640 | 0.043 | 5.0 | 2x32x32 = 2048 WI | ~2 % |
| 18 | C | ztrsm 576 LUNN | 50 | 0.039 | 771 | oneMKL | ~30 % (E) |
| 19 | C | k_det_phase l1207 | 96 | 0.034 | 351 | 1 x 64 | 0.05 % |
| 20 | U | block copies (H->D uploads, D->H gather) | 43200 | 0.029 | 0.7 | DMA | - |

Phases: U = updates (WRAPGR), S = stabilisation (wrapur/wrapul), C = CGR, O = observables. The two `corr` kernels with 1x1 are the `acc0` loops over `no = 1..Norb` (Norb = 1), a serial loop over N=576 sites in ONE thread (Prog/gpu_obs_mod.F90, lines ~437 and ~473).
Kernels that exist but are small are in `data/report_tables.md` / `data/j1/V7_pp3/rank0.log`.

L32 (32 ranks, 4 per tile), omp kernels only (C, job 560841, mean over 32 ranks, slow-down vs L1):
corr l450 3.86x, l386 3.51x, l418 3.32x, k_cpy (obs) 2.62x, k_cpy_n (wrap) 2.39x, vslice_r 2.19x; **the 1x1 and 1x64 kernels 1.05-1.07x**. I.e. the big kernels see the tile as saturated (~4 ranks share it: x3.3-3.9), the tiny ones do not feel the sharing.
omp-kernel device time per rank and sweep: O 1.72 s (L1: 0.84), C 0.19, U 0.09, S 0.04, shared 0.08.

## 3. Device busy vs wall per phase

Wall = TEAM_TIMERS of the unperturbed run (V7_n3, job 560831 and 560841). Busy = sum of device times of the phase from the tables above (event timestamps, isolation times), QR/getrf = duration of the synchronous library call (they are device bound inside, see below), shared kernels and the QR/unmqr calls split 50/50 between S and C by call count (96 + 96 per sweep).

**L1, V7 (s/sweep)**

| phase | wall | device busy | busy/wall | idle (gaps) | what the busy time is |
|--|--|--|--|--|--|
| U updates | 1.444 | 0.914 | 63 % | 0.53 (37 %) | rank-k gemm 0.31, hopping zgemm 0.32, block small kernels 0.25, omp wrap 0.04 |
| S stabilisation | 0.803 | ~0.64 | ~80 % | ~0.16 (20 %) | QR 0.41 + unmqr 0.13 (calls) + gemm 0.04 + kernels 0.06 |
| C CGR | 1.210 | 0.83 .. 1.21 | 69-100 % | 0 .. 0.38 | QR 0.41 + unmqr 0.13 + getrf 0.38 (hybrid: part on the host, unknown) + kernels 0.19 + trsm/gemm 0.11 |
| O observables | 1.180 | 1.047 | 89 % | 0.13 (11 %) | corr (1x1) 0.47, corr (512x1024) 0.31, zgemm 0.21, small 0.07 |
| sweep | 4.907 | 3.4 .. 3.8 | **70-78 %** | 1.1 .. 1.5 (22-30 %) | |

Occupancy-weighted ("full-tile-equivalent" seconds, weights from the launch geometry, E): U 0.70, S 0.15, C 0.21, O 0.55, total **1.6 s of 4.9 s = ~33 %**. The device is mostly *busy*, but with kernels that fill only a small part of the tile; only gemm, rank-k gemm and the big corr kernels fill it.

**V2, L1**: same U and O (1.48, 1.18), S 1.57 and C 2.03 s: the oneMKL hybrid geqrf (panel on the host CPU, 2.13 s/sweep in 192 calls of 11.1 ms) leaves the device idle during the host panel (C: ALFSH_TIMING, job 560831).

**L32, V7 (per rank, s/sweep)**: wall U 2.66-2.75, S 1.69, C 1.87, O 2.73, sweep 9.3-9.4 s (C, 560841/560848). Device busy per rank is not measurable (see tools). Indirect evidence (C+E):
- the big kernels are x3.3-3.9 slower than alone, the slice-level hopping gemm phase of U is 0.275 -> 0.919 s (x3.3, P3 timer): **tile saturated** in gemm and corr kernels;
- obs at L32: the full-fill kernels of 4 ranks need 4 x (0.31 corr + 0.21 gemm) = 2.1 s plus ~0.3 s small ones = 2.4 s of the 2.73 s wall: **observables are throughput bound at L32** (E, from isolation times);
- full-tile-equivalent demand of 4 ranks = 4 x 1.6 s = 6.4 s in 9.3 s wall = **~69 % tile utilisation (E)**; remaining headroom ~1.45x (matches the Phase-6 grounded ceiling 1.5-2.3x / realistic 1.3-1.6x).

## 4. Host round trips per sweep, V7, L1 (C unless marked)

| item | per sweep | cost |
|--|--|--|
| OpenMP target-kernel launches (each synchronous: the host waits) | 14,970 (U 3872, S 864, C 1152, shared 1882, O 7200) | host time in the launch call 3.1 us mean (PLUGIN_PROFILE); a synchronous launch+wait costs ~12 us (earlier bench, C, previous session) -> ~0.18 s (E) |
| D->H copies, OpenMP (`DataRead`) | 2,304 | 2.7 us host each |
| D->H copy of the block gather (+ wait), SYCL path | 8,640 (one per block) | **22 us each** (`blk_begin` 0.192 s / 8640), 0.19 s = 13 % of U wall; 30 us at L32 (0.235 s) |
| synchronous library calls (device QR, unmqr, getrf) | 480 (192 + 192 + 96) | 4.3 / 1.4 / 3.9 ms |
| H->D copies, OpenMP (`DataWrite`) | 45,856 | 5.3 us host each = **0.243 s** (device time 3.4 ms in total) |
| H->D copies, SYCL block path | 34,560 | **9.0 us submit each = 0.31 s** |
| SYCL block ops submitted (async): gather 8,640, fgather 17,280, scale 8,640, rank-k gemm 34,560, copies 43,200 | 112,320 | submit cost: gather 4.4, fgather 4.5, scale 4.4, **gemm 9.7**, copy 9.0 us (`ALFP7_SUBMIT`, 560845) = **0.88 s/sweep = 61 % of the update wall** (`blk_end` 0.78 s + `blk_begin` 0.21 s) |
| dispatched oneMKL calls (zgemm 576^3, ztrsm) | ~5,100 + 100 | 109 us device each; host submit+complete 124 us |

Update phase per block (8,640 blocks/sweep, 167 us wall each, C): host submission of ~13 async ops 91 us (`blk_end`), D->H round trip 22 us, slice-level hopping wraps 32 us (0.275 s/8640), host Metropolis/rest ~20 us. The device executes the block's ops in ~65 us (gather 5 + fgather 2x4.8 + scale 10.8 + 4 gemm x 9 + copies). So the update phase is **host-submission bound**: the device waits for the host (37 % idle), not the other way round. At L32 `blk_end` grows 0.784 -> 1.334 s (x1.7, C): host/driver contention for submissions.

## 5. Standard updates (Task 5): the GPU path supports them

`Sequential_scheme = "standard"` and `"delayed"` are supported on the device (`gpu_lr_mode` = LR_STD / LR_DELAYED / LR_SUBMATRIX in `gpu_lowrank_mod.F90`; standard = the per-flip device rank-1 update of ticket 205 in `gpu_update_mod.F90`). Same seeds, V7 switches, L1, 2 sweeps (C, job 560831):

| scheme | N_delay | s/sweep | upd s/sweep | stab / cgr / obs |
|--|--|--|--|--|
| standard | - | **18.32** | 14.81 | 0.81 / 1.24 / 1.19 |
| delayed | 32 | 19.66 | 16.15 | 0.80 / 1.24 / 1.19 |
| submatrix | 16 | 5.72 | 2.23 | 0.80 / 1.24 / 1.18 |
| submatrix | 32 (reference) | 5.01 (3-sweep timer: 4.91) | 1.51 | 0.81 / 1.24 / 1.18 |
| submatrix | 64 | 4.75 | 1.26 | 0.80 / 1.23 / 1.19 |
| submatrix | 128 (job 560848, 3 sweeps) | 5.14 | 1.63 | 0.82 / 1.22 / 1.19 |

Per flip (138,240 flips/sweep): standard 107 us, delayed 117 us, submatrix N_delay=32 10.9 us. Submatrix is 3.7x faster than standard in the update phase (9.8x), whole sweep 3.7x. Precision Green of all variants is at the usual 1e-11 level.
N_delay scan with V7 (alternated, 2 runs each, 3 sweeps; C, 560848): L1 32: 4.92, 64: 4.73-4.74, 128: 5.14; L32 32: 9.32-9.33, 64: 9.08-9.10 (-2.6 %), 128: 9.50.

## 6. Conclusion

**Single chain (L1).** The GPU is not mainly idle between launches: device busy 70-78 % of the wall (C+E), idle 22-30 %. What is wrong is the *fill*: only ~1/3 of the tile's capacity is used (E), because most busy time is spent in kernels that fill almost nothing:
- O (1.18 s): 0.47 s (40 %) in two **one-thread** kernels (`corr` l437, l473); busy 89 %.
- S+C: QR/LU calls = 1.46 s/sweep (30 % of the sweep, 3.9-4.3 ms per call): the device-QR panel is one work-group (0.8 % of the thread slots) and is 60 % of the geqrf time, apply/unmqr fill ~27 %, getrf is hybrid. Latency chains *inside* kernels, not launch gaps (the QR call takes 3.98 ms unperturbed and 4.27 ms with a host wait after every stage, i.e. waits between stages are not what costs the time).
- U: 63 % busy, 37 % idle, and the update is **host-submission bound** (0.88 s of host submit cost for 112k async ops vs 0.56 s of device time in the block path): that is the only phase where gaps (latency) dominate.

**Node layout (L32).** Idle time of one rank is filled by the other three: the tile reaches ~69 % of its capacity (E). The saturated parts are obs (corr + zgemm), the hopping-wrap zgemm and the rank-k gemm (bandwidth bound), i.e. **throughput/bandwidth**, while QR panels, 1x1 kernels, block gather/scale and the host submission remain latency/under-fill that sharing hides. The "70-90 % latency" Phase-6 fit is therefore right for the *under-fill + host* part at L1 but overstates it at L32; the remaining headroom at L32 is ~1.4-1.6x at most (E).

**Levers (E, not tried):**
1. **Observables, two serial kernels** (`corr` l437/l473, 1 thread, 0.47 s/sweep = 9.5 % at L1): parallel reduction over I. L1 gain ~ -0.4 s/sweep (-8 %); at L32 they overlap with other ranks (1.05x slow-down), so gain only ~1-3 % there. The big corr kernels (0.31 s alone, x3.5 at L32, only 576 outputs per launch on 512x1024 threads) are 3.5x slower than their data volume (16 MB) implies: a better mapping could cut the saturated obs phase at L32 (2.73 s of 9.3 s) by ~1/4 -> +5-8 % throughput.
2. **Update block path, fewer and cheaper host submissions**: pack the 4 H->D uploads per block into one (-0.23 s), fuse the 4 rank-k gemms of a block into one batched/own kernel (-0.2 s) or reuse a SYCL graph / command list: -0.4 .. -0.5 s/sweep at L1 (-9 %), more at L32 (blk_end is x1.7 slower there): +5-8 % throughput. `N_delay = 64` is the free part of this: -3.8 % (L1), -2.6 % (L32), measured (C).
3. **QR/LU call chains (1.46 s = 30 % at L1)**: 1-work-group panel (0.8 % fill) and hybrid getrf; only helps L1 (the other ranks hide it at L32) -> not a throughput lever, but the single-chain lever (device getrf, 2-3 work-groups per panel).
Throughput headroom beyond these is bounded by gemm/rank-k bandwidth saturation (~69 % used at L32, E).

## 7. Budget, jobs, files

Jobs (all single-node, preempt): 560821, 560823, 560825, 560826, 560829 (probes, ~1 min each), **560831** (L1: V7/V2 timing, PLUGIN_PROFILE, serialised block profile, INFO, schemes; 7.8 min), **560841** (L32 timing + PLUGIN_PROFILE; 4.7 min), 560843 / 560844 (MKL interposer, superseded; 3.7 min), **560845** (submit cost + dispatched BLAS profile; 2.0 min), 560846 (call sites + L32 serialised; hung, cancelled, 7.1 min), **560847** (call sites; 1.0 min), **560848** (N_delay scan; 5.5 min).
Total ~37 min = **~0.62 node-h** (from the start/end `date -Is` lines of the logs; limit 2 node-h).

Remote: `$SCRATCH/large_prompt/team/p7/` = `jobs/` (scripts, logs `*.sh.<id>.log/.err`), `src/` (measurement libraries), `lib/` (built .so), `work_*/runs/*` (raw rank logs of every run), `data/`.
Local: `docs/supermuc_large_prompt/phase2/team/p7/` = `jobs/`, `src/`, `data/` (rank logs and job logs of every run; `data/report_tables.md`, `data/kinfo_V7.txt`, `data/ALF.out.nm.txt`), `ana.py`, `tables.py`, `report.py` (tables from the raw data), this file.
Reproduce: `jobs/p7lib.sh` (helpers), `j1_single.sh` (L1), `j2_node.sh` (L32), `j4_off.sh` / `j6_site.sh` (serialised profile), `j7_nd.sh` (N_delay); `jobs/a_probe*.sh` = tool probes.
