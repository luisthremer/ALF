# Phase 8 / blk: fewer and cheaper host submissions in the GPU submatrix block path (ALFSH_BLKPACK)

2026-11-13, SuperMUC-NG Phase 2, ifx/icpx/oneMKL 2025.3.0. "C" = measured by a job of this session (job id), "E" = estimate.
V7 = `LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 ALFSH_ASYNCBLK=1 ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1`, V2 = same without the two DEVQR switches.
Reference L=24 (ref_L24), N_delay 32 and 64.  Layouts: L1 = one rank on tile 0 (OMP 8); T1 = 8 ranks, one per tile (OMP 14); L32 = 32 ranks, 4 per tile, CCS 4, OMP 3.

## Result in one paragraph

`ALFSH_BLKPACK=1` is bitwise identical to the current code (G1 at L6F/L6T/L12F/L24F, V7 and V2, CMP_TOL=0, planted defects fail) and cuts the block path from 13 to 8
async SYCL operations per block and the host-to-device copies from 4 to 1, **but it only helps where the host submission is the bottleneck: one chain, N_delay = 32
(update phase -10.7 %, sweep -2.7 %).  At the production setting N_delay = 64 it gives nothing at L1 / T1 and costs +1.5 % sweep time at the node layout (32 ranks, 4/tile)**
(in-loop timers, C).  `ALFSH_BLKPACK=2` (own fused kernel for the first flush GEMM, NOT bitwise, deviation <= 2.5e-15 relative at L24 same seed) is -17.5 % update at L1/N_delay 32
but again +/-0 .. +1.3 % sweep at the node layout with N_delay 64.  The expectation of the ticket (+5-8 % per card at the production layout) is **not met**: at 4 ranks per tile with
N_delay 64 the update phase is limited by the device (shared tile), not by the host submission cost; removing host work does not shorten it.  Recommendation: do not switch BLKPACK on in production
at N_delay = 64; it is useful (bitwise) only for a one-chain-per-tile run at N_delay <= 32.

## 1. What changed (files, behind `ALFSH_BLKPACK`, default off = the code before)

Only `Prog/gpu_lowrank_mod.F90` (+226/-5) and `Libraries/Modules/alf_la_sycl_blk.cpp` (+148).  Makefile, shim, obs and qr2 untouched.  Intel-only parts are inside the existing `ALF_LA_ONEMKL` guards; the
only edit outside a guard is in `blk_end` (vmin/vmax/nfm computed before the upload instead of after it: same values, used by all paths).  `GPU=none` is not affected (the module is `ALF_GPU` only).

Per block (8,640 blocks/sweep at N_delay 32) the async path before: copy bsd (H->D), gather, copy gbd (D->H), [wait], copy cbuf, copy bS_d, scale, then per flavour fgather, gemm1, gemm2, copy nsigma
= 13 ops, 4 H->D copies.  With BLKPACK=1:

| step | old | BLKPACK=1 | bitwise |
|---|---|---|---|
| block sites to the device | H->D copy of `bsd` | passed **by value** in the gather kernel (512 B struct, N_delay <= 64) | yes |
| gather -> host | gather to device + D->H copy | gather kernel writes the host-USM landing buffer directly | yes |
| wrap factors za/zb, block sites, pending sites S, packed M, new field values nsigma(vmin:vmax,t) | 4 H->D copies (cbuf 35 KB at nd 32, bS_d, nsigma, bsd) | packed into ONE host-USM staging buffer (compact layout, only the used part, every M segment 64-byte aligned), ONE copy | yes |
| scale + field values to nsigma_d | scale kernel, separate nsigma copy | the same scale kernel (same 3-D launch, same expressions) also copies the nns field values (work items fl=0, a=1, j<=nns) | yes |
| fgather per flavour | 2 launches, shared Qt buffer | 1 launch for both flavours, own Qt buffer for flavour 2 (same values) | yes |
| gemm1/gemm2 | 4 oneMKL zgemm | the same 4 oneMKL zgemm calls (same shapes/ops/scalars/ld; B of gemm1 at a 64-byte aligned offset) | yes |
| submission | `queue::parallel_for/memcpy` (event) | `sycl::ext::oneapi::experimental::parallel_for/memcpy` (no event), `ALFSB_NOENQ=1` = old call | yes |

Result: 8 ops, 1 H->D copy per block.  `ALFSH_BLKPACK=1` without `ALFSH_ASYNCBLK=1` does nothing; N_delay > 64 prints a message and falls back.

`ALFSH_BLKPACK=2` (separate sub-switch, **NOT bitwise**): additionally `alfsb_fgemm1` fuses the first flush GEMM (Y = Qt M^T) into the fgather kernel (one work item per (fl,i,m), k-sum in the order j = 0..kk-1 of
my kernel instead of oneMKL's zgemm; Qt is not stored): 6 ops per block.

Debug env (planted defects, the checks must fail): `ALFSH_BLKPACK_PLANT=1` wrong offset of the M segment of flavour 2, `2` wrong offset of the field values, `3` no wait after the gather, `4` M of the wrong
flavour, `5` last field value not scattered.  Ablations: `ALFSB_GDEV=1` (gather via device buffer + copy), `ALFSB_NOENQ=1`.

## 2. Build

`git archive HEAD` (7e8b6032) + the two changed files in `$SCRATCH/large_prompt/team/p8/blk/src_alf` (other agents' edits not included; md5 of both files equals the built copies); `jobs/c_build.sh` /
`jobs/build_comb.sh` (procedure of p5/c_build5.sh, new variant `aottxt` = AOT for PVC + `-DALF_TEAM_TIMERS` + text observables, one binary for all gates; `txtnew` = JIT text build as the p5 txtnew).
Final build job 561089 (2.6 min).  Both compile without errors.

## 3. Gates

### G1 bitwise (final binary, job 561091 = V7, 561092 = V2; same binary, same seed, fresh start, ALFSH_BLKPACK unset vs 1, CMP_TOL=0, all observable files + confout_0 + acceptance): all PASS, worst rel dev 0.000e+00
- V7: L6F (N_delay 8 and 4, the latter exercises k=1 -> kk=2), L6T, L12F (32), L24F (32 and 64).  V2: L6F (8, 4), L6T, L12F, L24F (32).  The log checks that the B run printed the BLKPACK line and the A run did not.
- Default build vs the previous binary (p5 txtnew), no ALFSH switches, L6F and L12F: PASS 0.0 (561092); same with ALFSH_ASYNCBLK=1 (new blk_end structure, BLKPACK unset): PASS 0.0.
- The same G1 passed on the first version of the code (jobs 560970/560973; the scale kernel was a 1-D launch then); the scale kernel, the plants 4/5 and BLKPACK=2 were changed after, so the final binary was fully re-run.
- BLKPACK=2 vs unset, V7, same seed (561091): worst rel dev 3.7e-16 (L6F), 1.4e-15 (L12F), 2.5e-15 (L24F), acceptance identical; not bitwise, as announced.

### G2 planted defects (561092, V7, B = BLKPACK=1 + plant, each at L6F and L12F): every one **FAILS** (10/10)
plant 1 (M offset +1 element, flavour 2): "Green function contains NaN", abort rc 1; plant 2 (field-value offset +1): "Unstable matrix" rc 3; plant 3 (no wait after the gather): rc 3; plant 4 (M of the wrong flavour):
rc 1 / rc 3; plant 5 (last field value of the block not uploaded): rc 3.  All abort within the first stabilisation (3-5 s), i.e. the observables are missing, G1 reports FAIL.  The check can tell: the correct code passes bitwise
where every one of these fails (no near-miss defect that survives was found; the check is not weaker than bitwise because it IS bitwise).

### G3 node Precision Green (32 ranks, 4/tile, CCS 4, OMP 3, V7, L24, 3 sweeps; job 561098)
rank 0 info: mean 1.963e-14, max 6.536e-10 for BLKPACK=0 and 1 (identical to all digits), N_delay 64; N_delay 32: 1.961e-14 / 6.536e-10 for both; BLKPACK=2: max 6.57e-10 / 6.56e-10 (same level, not bitwise).
Old vs new run directories (all 32 confout + merged observables) compared with compare.py tol 0: PASS (N_delay 64 and 32).  (Earlier job 560998: 6.536e-10 as well.)

### G4 timing (old and new alternating in one job; metric = in-loop TEAM_TIMERS, mean over ranks, s/sweep, 3 sweeps; 0 = BLKPACK unset)
The CPU-time-difference method (p6/m_tput.sh) was tried first (jobs 560998, 561035): between identical runs it scatters by +-10 % (e.g. 7.73 vs 9.42 s/sweep for the same setting, 1 vs 3 sweeps, node drift), useless for effects of 1-3 %.
The in-loop timers reproduce to 0.5 %.

**One chain on tile 0 (L1), V7, L24** (561099 for nd 32; 561035 for nd 64):

| N_delay | BLKPACK | update s/sweep | sweep s/sweep | block begin+end (host) us/block |
|---|---|---|---|---|
| 32 | 0 | 1.46 1.49 1.50 1.49 (mean 1.49) | 4.89 4.95 4.95 4.96 | 116.6 |
| 32 | 1 | 1.32 1.34 1.33 1.33 (mean 1.33; -10.7 %) | 4.77 4.81 4.80 4.80 (-2.7 %) | 99.3 |
| 32 | 2 | 1.22 1.23 (-17.5 %) | 4.68 4.68 (-5.1 %) | 87.2 |
| 64 | 0 | 1.23 1.25 | 4.66 4.70 | 149.2 |
| 64 | 1 | 1.25 1.25 | 4.68 4.70 | 151.6 |
| 64 | 2 | 1.21 1.22 | 4.66 4.66 | 143.5 |

(Earlier job 560980 with the first code version: nd 32 0: 1.47-1.50, 1: 1.34, nd 64: 1.26/1.25 -- same picture.)  With BLKPACK the host no longer waits at blk_end (0.79 -> 0.43 s/sweep) but waits at blk_begin for the
device (0.21 -> 0.43): the chain is device-latency bound after the change (gather, scale, fgather, 4 gemm in series on the in-order queue).  At N_delay 64 the block is longer on the device (more k) and the update phase
was already device bound: no change.

**One rank per tile (T1: 8 ranks, OMP 14), V7, L24, N_delay 64** (561035): update 1.427 (0) / 1.417 (1) / 1.376 (2), sweep 5.157 / 5.162 / 5.120 s: -0.1 %, -0.7 %.

**Node layout (32 ranks, 4/tile, CCS 4, OMP 3), V7, L24, in-loop means over 32 ranks** (561098; ABBA order 0 1 2 2 1 0):

| N_delay | BLKPACK | update | stab | cgr | obs | sweep | sweeps/s/card (8 chains per card / sweep) |
|---|---|---|---|---|---|---|---|
| 64 | 0 | 2.464 2.463 | 1.708 1.704 | 1.874 1.871 | 2.747 2.738 | **9.103 9.092** | 0.879 |
| 64 | 1 | 2.536 2.518 | 1.731 1.729 | 1.886 1.893 | 2.797 2.775 | 9.258 9.222 (+1.6 %) | 0.866 (-1.5 %) |
| 64 | 2 | 2.507 2.474 | 1.737 1.734 | 1.886 1.897 | 2.808 2.776 | 9.246 9.188 (+1.3 %) | 0.868 (-1.3 %) |
| 32 | 0 | 2.722 2.737 | 1.697 1.702 | 1.872 1.866 | 2.727 2.745 | 9.321 9.354 | 0.858 |
| 32 | 1 | 2.579 2.572 | 1.709 1.709 | 1.875 1.880 | 2.757 2.746 | 9.225 9.214 (-1.4 %) | 0.869 (+1.3 %) |
| 32 | 2 | 2.477 2.468 | 1.715 1.711 | 1.872 1.872 | 2.766 2.758 | 9.137 9.116 (-2.3 %) | 0.877 (+2.2 %) |

N_delay 64 (production) without BLKPACK is the best of all (0.879 per card).  At 4 ranks per tile the update phase is bound by the shared tile (device time of gemm/rank-k, E, see p7 section 3), so fewer host
submissions only help at the small block size, and at N_delay 64 the new path is 1-1.5 % slower in update AND (unexplained, E: the obs/stab phases are 1-2 % slower too, probably less overlap between ranks'
phases once the update phase is shorter/different) in the rest.  Differences of 1-2 % are at the edge of what the in-loop timers resolve (spread between repeats <= 0.7 %), so "no gain" is the safe reading for the node layout at N_delay 64.

**Ablations of BLKPACK=1 (L1, nd 32, 561099):** update 1.32-1.34 (default), 1.29 with `ALFSB_GDEV=1` (gather via a device buffer and a copy: the direct host-USM write is not faster), 1.31 with `ALFSB_NOENQ=1`
(the enqueue functions give nothing measurable), 1.29 with both.  Not tested at the node layout (budget).

## 4. Why the expectation was not met (E, supported by the numbers above)
- The profile of p7 (host submit 91 us vs device 65 us per block at L1, nd 32) was right for that case; the fix works there (host blk_end 0.79 -> 0.43 s/sweep, update -10 %), the remaining per-block time (99 us) is the
  in-order chain of 8 small device ops (each ~5-10 us plus launch gaps), i.e. device latency; BLKPACK=2 shows the chain still shortens with fewer ops (-0.1 s more).
- At N_delay 64 and at 4 ranks per tile the block path is device/throughput bound; host cost reductions cannot show up.  The rank-k gemms (4 oneMKL calls, 9 us each, 40 % of the chain) are the only large item left and cannot be
  fused bitwise (p6: oneMKL batched/wide gemm is not bitwise equal to single calls).
- SYCL command graph (checked, not used): the probe `jobs/graph_probe.cpp` (jobs 560966/560969) died at queue creation on the compute node (`UR_RESULT_ERROR_UNINITIALIZED`, also under mpiexec) before reaching the graph code, so
  feasibility inside the OpenMP interop queue is **not tested**; it would in any case need per-block kernel arguments and gemm shapes (k) that change every block, i.e. graph node updates through oneMKL calls, which oneMKL does not offer.  Judged not worth the remaining budget.
- Not done: the 45,856 small OpenMP `target update to` copies (5.3 us each, 0.24 s/sweep) are not in the block path of gpu_lowrank_mod (the ab path has none); they come from other modules (files I may not edit).

## 5. Budget, jobs, files
Jobs (all single-node, preempt, all logs in `logs/`): builds 560963, 560986, 561089; probes 560966, 560969 (graph); first versions of the gates 560970 (g12), 560973 (g1b), 560980/560997 (g4a, 560997 preempted at 12:02),
560998 (g4n, hit its time limit, CPU-time metric); final gates 561091 (G1 V7 + BLKPACK=2), 561092 (G1 V2 + default-vs-previous + G2); timing 561035 (L1 nd 64 + T1), 561098 (node, in-loop), 561099 (L1 nd 32 + ablations).
Node time from the start/end lines: about 3.9 node-h of the 4 (builds 0.12, gates 1.35, timing and its repeat after the CPU-time method failed 2.4).  No job was started after the budget was reached.

Local: `docs/supermuc_large_prompt/phase2/team/p8/blk/` = this file, `jobs/` (build_comb.sh, c_build.sh, g1.sh = final G1/G2, g12.sh/g1b.sh = first versions, g4a*.sh/g4n.sh = first timing versions with the CPU-time metric,
t4n.sh/t4l.sh = final timing, p8lib.sh, wrap_p.sh, inloop.py, graph_probe.*, sync.sh), `logs/` (all job logs, `p3_summary.txt` = block timers).  Remote: `$SCRATCH/large_prompt/team/p8/blk/` (src_alf, build_aottxt, build_txtnew, work_*, jobs).
Sources changed (not committed): `Prog/gpu_lowrank_mod.F90`, `Libraries/Modules/alf_la_sycl_blk.cpp`.
