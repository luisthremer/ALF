# Phase 6 design note: B walkers in one MPI rank (batched multi-walker GPU path)

Status 2026-11-12. "C" = measured by our own job on SuperMUC-NG Phase 2 (jobs and files named); "E" = estimate.
Measurements: `ceil_bench.F90` (jobs 560781, 560787, ...), `m_tput.sh`/`m_flop.sh` (560772/3), prototype `m_proto.sh`/`m_gate2.sh`.
Raw outputs: `data/`. Summary script: `ceil_sum.py`.

## 0. What the measurements say (read this first)

1. **Gate 1 (ceiling >= 1.5x) passes on the task's formula, and the formula is not informative.** One tile does 12.8 TFLOP/s on a single 576^3 zgemm when it is the only busy tile, 16-18.6 with a batch
   of walkers (C, job 560781/560787); with all eight tiles busy: 10-11 and 12.4; at 1024^3: 18.5 and 20.0 alone, 14-15 busy. Floor = 9.2 TFLOP / best rate = 0.55 s (alone) to 0.74 s (node loaded) per
   walker-sweep, card ceiling 2 tiles / floor = 3.6 to 2.7 sweeps/s = 4.1x to 3.1x today's 0.879. But the formula puts all 9.2 TFLOP at gemm speed. They are not: the rank-k flush gemms
   (0.49 TFLOP, 34,560 calls, 366 GB of G traffic per walker-sweep) are memory-bandwidth bound at production load (C, job 560801: 4 processes per tile move ~0.8-1.1 TB/s of the ~1.6 TB/s per tile, E, data sheet;
   batching them gains 1.0-1.3x there), QR/LU/trsm and the observables are latency bound.
2. **A grounded ceiling (E).** Tile time per walker-sweep today is 2.28 s (2 tiles / 0.879 per card). Floor of the throughput-bound part: big zgemm 8.2 TFLOP / 12.4 = 0.66 s + flush 366 GB / 1.2 TB/s = 0.30 s
   + QR/LU/trsm/observables about 0.2 s (E) = about 1.2 s, so at most 1.9x per card even if every launch and sync were free. A second estimate: fitting the phase times of one chain with and without CCS 4 (C, job 560792: sweep 5.32 s on a
   full tile, 8.31 s on one of four compute slices) as time = latency + throughput x (1 or 4) gives a throughput part of ~1.0 s per walker-sweep (E, exact if the slice slows the throughput part by 4; a factor 3 gives 1.5 s)
   and a ceiling of 1.5-2.3x. A realistic implementation reaches part of that: **1.3-1.6x per card (E)**.
3. **Batching the big gemms alone buys 2-7 %** (C + E). The prototype (section 9) measures the hopping wrap of WRAPGRUP/WRAPGRDO: 0.34 ms per call unbatched on a loaded node (1 process per tile) and 0.27 ms per walker at B = 4
   (-22 %), 0.225 at B = 16 (-33 %); with CCS 4 (one compute slice) 0.845 -> 0.733 ms (-13 %). That part is 6-10 % of a sweep, so -2 %; adding the other big gemms (stabilisation wraps, CGR, decompose) -> 2-7 % (E).
4. **oneMKL's batched and wide gemm do not give the bits of the single gemm, at almost every size** (C, job 560787, bit comparison against the single-call result on the same inputs). Equal to the single gemm:
   N = 1024 for every B; N = 256 only for B = 2; N = 36 for B <= 4. Different in 99.7 % of the elements: N = 144, 400, 576 for every B >= 2, N = 256 for B >= 4. The result also depends on B (no size rule). In the real sweep
   (Gate 2, job 560791) batched mode 1 vs the unbatched default differs by 1.2e-15 (L12F) and 2.2e-15 (L24F) in the worst observable of a short run, so "walker b of a batched run is bitwise identical to the same walker alone" holds
   only when the gemms stay one call per walker (mode 3: bitwise PASS at L6F, L12F, L24F, but no gain: 0.348 vs 0.353 ms per call).
5. **The throughput lever is latency, not FLOPs.** Single chain (V7, 1 rank per tile, no CCS, C, job 560772/560792): 5.33 s/sweep = updates 1.64 + stabilisation wraps 0.90 + CGR 1.23 + observables 1.31 + rest 0.25.
   The CCS fit puts 70-90 % of every phase into latency (updates 83 %, wraps 77 %, CGR 89 %, observables 71 %, E). The observables take 1.3 s here against 0.12 s of device time on the AMD (E, their number).
   A walker index in these kernels lets one launch serve B walkers; that, not the batched gemm, is what a rebuild would buy.

## 1. Execution model

One rank owns B walkers (chains). All walkers advance in lockstep through the same (sweep, direction, time slice, block of 32 sites, stabilisation block):
- host orchestration (main.F90, Wrapgr_mod, gpu_update_mod, gpu_stab_mod) stays one program; every host-only step becomes `do w = 1, B`;
- every device step that is today one launch / library call becomes one launch with a walker index (own kernels) or B back-to-back library calls (oneMKL, see 6);
- the walkers are independent chains with their own fields, RNG, Phase, observables, bins and files; they share only read-only tables
  (hopping tables, `Hd`, V_exp when the parameters are equal) and the device scratch dimensioned for B.

Within one slice the site updates of one walker are sequential (each decision needs G after the previous accept). In lockstep, step s of all B walkers is
issued together: one batched gather / one D->H of B scalars, B host decisions (cheap, 8 % of the update phase today, C: 0.35 of 4.12 s per two sweeps),
one batched H->D of the accept flags, one launch applying the accepted rank-1 steps of all walkers (walkers that rejected get a no-op via a per-walker flag).

## 2. Device arrays: a walker index

The walker index is the LAST index, so the data of one walker is contiguous and `G_d(1,1,nf,w)` is a legal argument of every existing routine (the
single-walker code stays callable per walker, which is what the self-tests and the B = 1 reference use).

| array (module) | today | with walkers | size per walker, L=24 / L=32 (E) |
|---|---|---|---|
| G_d (state) | (N,N,N_FL) | (N,N,N_FL,B) | 10.6 / 33.6 MB |
| udw_U, udw_V, udw_D / udw_L (state) | (N,N,2,N_FL_eff) | + B | 21 / 67 MB (x N_FL_eff) |
| udv_U, udv_V, udv_D / udv_L (state) | (N,N,NSTM,N_FL_eff) | + B | 255 / 1070 MB (x N_FL_eff) |
| nsigma_d, nsigma_t (state) | (N_op,Ltrot), (N_op) | + B | 2.2 / 4.9 MB |
| scratch scr_c, scr_vec, scr_r, scr_work, scr_ipiv (state) | one set | + B | 2x5.3 + ... |
| sA, sB, sC, st_Gn, st_Gt, st_tau, st_Dl, st_c, st_r, st_nrm, st_Lt, st_ipvt, st_perm, st_ipiv (stab) | one set | + B | 5x5.3 MB + vectors |
| X_d, Y_d, gbuf, S_h, slot_h, M_h/M_d, jc, Gb, bidx, bs_h ... (lowrank) | per flavour | + B (the pending-column counters `jc` are per walker) | small |
| Gsv, Gw, Ts, Asy, Bsy, Tmp, acc, tacc (obs) | one set | + B | ~ 4 G |
| hop_*, V_meta, V_P, vr_ptr, vr_ent, Hd (state, wrap) | one | shared | 42 MB |
| V_exp | one | shared if the coupling is equal for all walkers of the rank, else (.., B) | 0.1 MB |

Total about 0.3-0.6 GB per walker at L=24 and 1.2-2.3 GB at L=32 (E): B = 8 fits a 64 GB tile with room to spare. Ltau > 0 (gpu_tau_mod, 2627 lines) is
refused with B > 1 (not needed at Ltau = 0, which is the production setting).

## 3. Which phases batch (per walker-sweep, L=24)

| phase | today | batched form | bitwise vs alone | gain (E) |
|---|---|---|---|---|
| G wrap (WRAPGRUP/DO hopping part): 2 gemm + copy per flavour and slice | 2 la_gemm (N^3 with the shared Hd) + k_cpy_n, td kernels | wide gemm (n = B N) or strided batch for the Hd product (left), strided batch sB = 0 (right), one copy kernel, td kernels with walker index | gemm: NOT at N=576 (C, Gate 2 F); kernels: yes by construction | -22 % of this part = about -2 % of a walker-sweep (C, section 9) |
| Op_V wrap `k_vwrap_all`, `k_site_wrap` | one kernel per slice and flavour | grid over (r, c, w), field read from nsigma_d(n,t,w) | yes (same expression per element) | launch share |
| site update: gather, G(p,p), accept (rank-1) | host block mode: blk_begin/blk_end device sequence per block of 32 sites | block gather/scale/flush kernels with walker index; host decisions in a loop over w; flush gemms per walker (see 6) | yes if the flush gemm keeps its shape (k per walker) | the update phase is 1.63 s/sweep (C), 156 us per block end |
| stabilisation wraps (wrapur/ul): per slice 1 N^3 gemm + Vslice kernel | per walker | as the G wrap: left gemm, wide or strided, Vslice kernel with w | as above | 0.86 s/sweep (C) |
| decompose (QR with column-norm pivot, trsm, scaling) | libalfdevqr: panel kernel (1 work-group) + trailing apply, per matrix; pivot vector D->H for `Pivot_phase` | batch dimension in the panel and apply kernels (work-group index = matrix); one D->H for B pivot vectors; scaling kernels with w | yes (same arithmetic per matrix) | part of stab 0.8 + CGR 1.25 |
| CGR (2 gemm, LU, QR, ormqr, trsm, gemm) | per walker | gemm with distinct operands: strided batch (DD, 15-18 TFLOP/s at B = 16, C; bits not compared) or one call per walker; LU/trsm: per walker calls (getrf via oneMKL, hybrid) | gemm: see 6; rest yes | CGR is launch bound (1.25 s/sweep, C) |
| observables (gpu_obs_mod): 1.26 s/sweep | many small kernels per measurement | walker index in all kernels, accumulators acc(.., w) | yes | largest single latency lever (E: AMD needed 0.12 s) |
| precision checks, phase, Control_* | host / small kernels | per walker | yes | small |

## 4. Monte Carlo: accept/reject and RNG per walker

- **RNG.** `Random_Wrap` is the Fortran intrinsic RANDOM_NUMBER with a process-global state (Ranset / Ranget = RANDOM_SEED PUT/GET). Per-walker streams:
  keep `rng_state(K,B)`; the walker loop does `RANDOM_SEED(PUT=state(:,w))` before the host steps of walker w and `GET` after (2 calls, ~0.1 us each, once per walker per
  host step, not per draw). Draw order per walker is unchanged: the draws are `ranf_wrap()` in Wrapgr / the update (one per site step plus the T0 test), none in device code.
  Walker w of rank r is seeded like rank r*B + w of today's run (same seeds file line), so "walker alone" = today's single-walker run of that rank.
- **Decision.** Metropolis accept/reject stays on the host per walker (needs G(p,p) from the device: one batched D->H of B scalars per step instead of B copies). The field value and the
  accept flag go back in one batched H->D. Nothing couples walkers, so lockstep issuing cannot change any walker's trajectory.
- **Determinism rule.** Every device kernel computes, per element, the same expression in the same order as the single-walker kernel (true by construction when the walker index
  only offsets the base address). The library gemms are the exception (6).

## 5. MPI, PARALLEL_PARAMS, I/O

- Logical walker id g = rank*B + w. `mpi_per_parameter_set`, Group_Comm, `igroup` (Temp_<g>/ directories under PARALLEL_PARAMS + TEMPERING) and the seeds file are indexed by g; the
  MPI_COMM_SPLIT colours are computed on g, so the number of groups is B x ranks. Tempering swaps between two groups on the same rank become a local exchange of the fields
  (host copy of nsigma + `gpu_upd_resync_fields` for the two slots), between ranks they stay MPI messages on the walker's own buffer.
- Info, confout and bins: B file sets per rank (Temp_<g>/info, confout_<g>, data), written by the rank at the end of each bin; HDF5 serial is called B times. Observable reductions
  (Obs_* accumulation) are per walker, so the statistics of a bin are exactly those of B ranks with one walker.
- Environment: one rank per tile with CCS 4 would give the process one of four compute slices (measured: a single process then gets 3.9 TFLOP/s on zgemm 576^3, a quarter of the tile
  11.1, C job 560781 r8 vs r8ccs). So the productive layout is 4 ranks per tile x B walkers each (CCS 4, as today), and B up to 4 raises walkers per tile to 16. 8 processes per tile is
  not an option (0.261 per card, C, job 560772).

## 6. Bitwise identity and the oneMKL gemm (the finding that shapes the plan)

C (job 560787 `data/b*.out`, 560781 `data/r*.out`): the results of `zgemm_batch_strided` (shared operand), of the wide gemm (n = B N) and of the tall gemm are compared bit by bit with the single
`zgemm` on the same inputs: equal only at N = 1024 (all B), N = 256 (B = 2), N = 36 (B <= 4); at N = 36 (B >= 8), 144, 256 (B >= 4), 400, 576 99.7 % of the elements differ (rounding level). The distinct-operand strided batch
(CGR products) was timed but not bit-compared. oneMKL picks a different kernel (different K blocking) for batched calls, and which one depends on N and B. Consequences:
1. **Gate "walker b batched == walker b alone" cannot hold with a batched oneMKL gemm at L = 24.** B-independence of the batched result is also size-dependent (wide(B) vs wide(2) equal at N = 576, strided(B) vs strided(1) not),
   so even a new "batched-mode" reference is fragile: it would have to be re-verified for every (N, B, oneMKL version).
2. **Bitwise-safe option (recommended if bits matter):** one oneMKL call per walker with the same shape as today, issued back to back and asynchronously on the in-order queue (the shim's async-block path shows the
   pattern: 4.3 us vs 12.3 us per op). It gains only the host wait (a part of the ~50-100 us launch+sync), no kernel-rate gain; all own kernels (walker index) stay bitwise identical by construction.
3. **Numerical-variant option:** batched oneMKL with the owner's approval of a rounding-level change relative to the legacy path (the default build unchanged), tested same-seed against the CPU at 1e-13 (L6F, L12F)
   and with B fixed per production run. Gain on the gemm part: +23 % rate at 4 processes per tile (C: 12.4 vs 10.1 TFLOP/s per tile), about +7 % of the total throughput (E).
4. **Own gemm kernel:** not recommended (our device QR needed Phases 4 and 5 to reach 3.9 ms at N = 576; a zgemm at 12 TFLOP/s is out of reach).
5. The rank-k flush gemms (k = pending count, 4 x 8640 per sweep, mean k ~ 5) have per-walker k: they stay one call per walker in every option.

## 7. Effort per module (E, person-days for one developer who knows the code; includes unit tests against the single-walker path)

| module (file, lines) | work | days |
|---|---|---|
| gpu_state_mod.F90 (1893) | walker dimension on ~25 arrays, B-aware setup/finalize, per-walker mirrors of the fields, `gpu_supported` refusals (Ltau, td dense, Global moves) | 3 |
| gpu_wrap_mod.F90 (2645) | walker index in the ~40 production kernels (hop groups, dense, el_left/right, vop, vwrap_all, td_apply); batched Hd products (wide/strided, measured in `ceil_bench`); keep self-tests for w=1 | 6 |
| gpu_lowrank_mod.F90 (1696) + gpu_update_mod.F90 (633) + alf_la_sycl_blk.cpp (102) | lockstep site step: batched gather / G(p,p) / accept, per-walker pending counters and flushes, block mode begin/end for B, async queue submission (SYCL kernels get a 4th range dimension) | 12 |
| gpu_stab_mod.F90 (1905) | per-walker loops first (simple), then batched decompose/CGR kernels with w, one D->H of B pivot vectors | 7 |
| alf_la_sycl_qr2.hpp (452), alf_la_sycl_devqr.cpp (95), alf_la_sycl_shim.cpp (304) | batch dimension in panel/apply kernels and entry points `alfdevqr_zgeqrf2_batched`, `_zunmqr_batched`; the Phase 4/5 race (tau, exchange, v buffer slots) has to be re-proven with batch work-groups: planted-defect and node tests again | 8 |
| alf_la_mod.F90 (2105) | `la_gemm_batched` exists (strided); add per-walker-pointer variant and B-aware workspace queries; B back-to-back async gemm helper | 3 |
| gpu_obs_mod.F90 (713) | walker index in all kernels, acc(.., w), per-walker flush | 4 |
| main.F90 (1460), Wrapgr_mod, Control, Random_Wrap | walker loop around host-only parts, per-walker derived types (udvl, udvr, udvst, GR, Phase, nsigma, obs), RNG swap, B file sets, tempering with local swaps | 9 |
| validation | walker b vs alone (bitwise per phase self-test and full run, L6/L12/L24), same-seed vs CPU, node Precision-Green, race test with B, throughput runs, production dry run | 10 |
| **total** | | **~62 days (about 12 weeks)** |

The prototype of this task (sections 8, 9) covers only the batched hopping wrap of WRAPGRUP/WRAPGRDO and costs about 3 of those days; the rest of the table is untouched.

## 8. What the prototype is (and is not)

`gpu_wrap_mod.F90` (`gpu_wb_*`), switched by `ALF_GPU_WB=B` (default off, `ALF_GPU_WB_MODE` 1 = strided batch both sides, 2 = wide gemm left + strided right, 3 = separate gemms as control).
It is not a multi-walker engine (ALF still runs one chain per rank). The hopping part of WRAPGRUP/WRAPGRDO (mmthr + mmthl_m1, resp. mmthl + mmthr_m1) is applied to B slots of a batch array; the
real chain's G sits in slot s (s rotates with every call so every slot position is exercised), the other slots hold scaled copies (stand-in walkers) and are discarded; the real slot is copied back.
The time-dependent hopping elements (`td_apply`) are applied per slot, not batched. `ALF_GPU_WB_TIME=1` prints per-call wrap times (baseline and batched); `ALF_GPU_WB_PLANT=1` is a planted defect
(1e-12 relative on one element) for the negative control of the bitwise gate.

## 9. Prototype results (L=24, V7 switches: staging buffer, ALFSH_ASYNCBLK=1, ALFSH_DEVQR=2, ALFSH_DEVUNMQR=1; ALF_GPU_DEBUG=1 is needed for any ALF_GPU_* variable, allow-list only)

**Gate 2 (job 560791, same txtnew binary, CMP_TOL=0, compare.py on all observables, confout and acceptance; `data/gate2.560791.txt`):**

| case | L6F | L12F | L24F (1 sweep) |
|---|---|---|---|
| A control: default vs default | PASS 0 | PASS 0 | PASS 0 |
| B mode 3 (B = 4 separate gemms, stand-in slots rotating) vs default | **PASS 0** | **PASS 0** | **PASS 0** |
| C planted defect (1e-12 on one element of G) vs default (must fail) | FAIL 4.8e-13 | FAIL 6.9e-13 | FAIL 5.0e-13 |
| D mode 1 (strided batch) B = 4 vs mode 1 B = 1 | PASS 0 | FAIL 1.2e-15 | FAIL 2.2e-15 |
| D mode 1 B = 16 vs mode 1 B = 1 | FAIL 4.2e-16 | FAIL 1.7e-15 | FAIL 2.2e-15 |
| E mode 2 (wide + strided) B = 4 vs B = 2 | PASS 0 | FAIL 2.5e-15 | PASS 0 |
| E mode 2 B = 8 vs B = 2 | FAIL 3.7e-16 | FAIL 1.4e-15 | PASS 0 |
| F mode 1 B = 4 vs default (the Gate 2 of the task) | PASS 0 | FAIL 1.2e-15 | FAIL 2.2e-15 |

Reading: **Gate 2 is met only by the bit-safe mode 3** (per-walker gemm calls; this also shows that the stand-in slots, the rotating slot position, the copies in and out do not leak into the real chain). It is not met by any batched oneMKL call
(F), and a "batched mode alone" reference does not rescue it (D). At L24 (N = 576) the wide/strided result is independent of B >= 2 (E), at N = 144 and with B = 8/16 it is not. The check can tell right from wrong: the planted defect fails
at 5e-13, the differences of D/E/F are 4e-16 to 2.5e-15 (rounding level). Control A and B pass with exactly 0, so CMP_TOL=0 is a meaningful bitwise test.

**Wrap time per call** (one flavour: 2 zgemm 576^3 + 2 copies + td kernels + for batched: stand-in fill and copies in/out; 1920 calls per 2 sweeps; `data/proto.560792.txt`):

| condition | baseline unbatched | B=1 | B=2 | B=4 | B=8 | B=16 | mode 2 B=4 | mode 3 B=4 |
|---|---|---|---|---|---|---|---|---|
| 1 process per tile, no CCS, all 8 tiles busy, **ms per walker** | 0.340 (mean of 8 tiles, 0.30-0.36) | 0.403 | 0.298 | 0.266 (8 tiles, 0.233-0.278) | 0.255 | 0.225 | 0.295 | 0.348 |
| same, change vs baseline | | +18 % | -12 % | **-22 %** | -25 % | **-33 %** | -13 % | +2 % |
| 1 process per tile, CCS 4 (one compute slice), ms per walker | 0.845 | | | 0.733 (-13 %) | | | | |
| 4 processes per tile (32 ranks, CCS 4), ms per call and rank (4 walkers share the tile) | 0.93 (0.88-0.97) | | | | | | | |

Of the batched time per call at B = 4, fill + copy in/out is 0.046 ms (0.012 per walker), td kernels 0.002 ms (they are not batched here). The Precision Green of every run passes (max 3.3e-11, node run 6.5e-10, same as the default).
Share of the hopping wrap in a sweep: 960 calls x 0.34 ms = 0.33 s of 5.32 s (6 %) alone, 960 x 0.93 ms = 0.89 s of 9.2 s (10 %) in the node configuration (C).

**Implication.** Batching only this part changes a walker-sweep by -2 % (alone: 22 % of 6 %), -3 % with B = 16. The default build with the prototype compiled in is bitwise the pre-prototype build (job 560805, section 10).

## 10. Default build unchanged, budget, recommendation

- Default build with the prototype compiled in (no switches) vs the pre-prototype build (p5 txtnew), same seed, V7 switches, CMP_TOL=0: L6F, L12F, L24F all PASS, worst deviation 0 (job 560805, `m_gate0.sh`).
- Budget: about 1.3 node-h (jobs 560780-560805 incl. one cancelled; sum of job run times 74 min), of 15.
- **Recommendation (E where marked): no-go on the full rebuild as it is specified now (12 weeks, ~62 person-days).** It was justified by the big-gemm ceiling (3-4x), but that ceiling is not reachable (flush gemms are bandwidth bound,
  QR/LU/observables latency bound), batching the big gemms buys 2-7 % and cannot be bitwise at L=24 with oneMKL batched calls, and the grounded ceiling of the whole idea is 1.5-2.3x with 1.3-1.6x realistic (E), for 12 weeks of work
  plus the re-proof of the device-QR race with batch work-groups.
- **If the owner wants the latency lever, do a pilot first (about 15-20 person-days, E), bit-safe (one oneMKL call per walker, walker index only in our own kernels):** observables (gpu_obs_mod, 1.3 s/sweep, 71 % latency) and the block
  begin/end of the update (gpu_lowrank_mod + alf_la_sycl_blk.cpp, 1.64 s/sweep, 83 % latency), B = 4 walkers per rank with CCS 4 (16 per tile). Decision criterion measured on the pilot: observables + update time per walker at B = 4 <= 0.45 of today's
  (that projects to ~1.4x per card, E) and bitwise walker-b-vs-alone at L6/L12/L24. If it fails, stop; if it passes, extend to stabilisation/CGR.
- Cheap things to try before any code (E, untested here): launch-latency settings on the existing binary (Level Zero immediate command lists, `SYCL_UR_USE_LEVEL_ZERO_V2=1`, `LIBOMPTARGET_LEVEL_ZERO_USE_IMMEDIATE_COMMAND_LIST`) since 70-90 % of every phase is latency.
- What this prototype did not measure: 4 processes per tile each with B = 4 (needs the engine), batched observables / update / stabilisation (not written), the td kernels batched, B > 16.

## 11. Appendix: Gate 1 raw numbers (C; ceil_bench via alf_la ONEMKL, host-synchronous calls, omp_get_wtime over >= 12 repetitions; per tile, TFLOP/s)

Columns: R1 = one process, tile alone (job 560781 r1); R8 = 1 process per tile on all 8 tiles, no CCS; R8c = same with ZEX_NUMBER_OF_CCS=4 (the process gets one of four compute slices);
R32 = 4 processes per tile, CCS 4, aggregate per tile (as production). `data/*.sum.txt`, raw `data/*.out`.

| N = 576 | R1 | R8 | R8c | R32 |
|---|---|---|---|---|
| single zgemm, same operands repeated | 12.86 | 11.09 | 3.91 | 11.94 |
| B separate back-to-back zgemm on B different matrices (B = 16) | 12.39 | 9.57 | 3.91 | 10.09 |
| gemm + copy kernel (the wrap's hopd_left) | 11.13 | 8.92 | 3.68 | 10.16 |
| strided batch, shared A, B = 4 / 16 | 16.09 / 16.44 | 12.13 / 12.22 | 4.53 / 4.90 | 12.37 / 12.17 |
| wide gemm n = B N, B = 4 / 16 | 15.83 / 16.72 | 12.03 / 12.42 | 4.54 / 4.92 | 12.42 / 11.75 |
| tall gemm (right side), B = 16 | 16.66 | 12.34 | 4.92 | 11.90 |
| strided batch, distinct A and B, B = 16 | 15.50 | 11.29 | 4.88 | 11.44 |
| wide gemm + copy kernel, B = 8 | 14.81 | 11.24 | 4.67 | 10.65 |
| **N = 1024**: single / wide B = 8 / strided B = 16 | 18.55 / 20.05 / 18.64 | 15.78 / 15.69 / 13.19 | 4.98 / 5.09 / 4.99 | 14.27 / 14.21 / 13.01 |

Floor = 9.2 TFLOP / best rate; card ceiling = 2 x rate / 9.2 TFLOP:

| basis | rate | floor (s per walker-sweep) | card ceiling (sweeps/s) | vs 0.879 |
|---|---|---|---|---|
| tile alone, best batched 576 | 16.7 | 0.55 | 3.63 | 4.1x |
| node loaded (R8/R32), best batched 576 | 12.4 | 0.74 | 2.70 | 3.1x |
| node loaded, single call (R32, back-to-back distinct operands) | 10.1 | 0.91 | 2.20 | 2.5x |
| grounded estimate (section 0, item 2, E) | | ~1.2 | 1.7 | 1.9x |

Rank-k flush gemm (C -= X Y^T, N = 576, `K_*` lines, job 560801), microseconds per gemm:

| k | R1 single / strided B=4 / B=16 | R8 single / B=4 / B=16 | R32 single / B=4 / B=16 |
|---|---|---|---|
| 8 | 21.7 / 9.6 / 6.1 | 22.4 / 9.8 / 6.5 | 37.4 / 28.2 / 53.5 |
| 32 | 26.9 / 13.9 / 8.7 | 27.6 / 15.0 / 10.5 | 57.9 / 49.4 / 72.9 |

On an idle or one-per-tile device batching helps 2-3x (launch bound); at production load (R32) it gives 1.0-1.3x or less (memory bandwidth bound: 4 x 10.6 MB per 37-58 us).
