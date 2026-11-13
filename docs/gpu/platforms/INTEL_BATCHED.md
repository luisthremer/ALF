# Intel Max 1550 vs H100: regime of each phase, upper bounds, batched multi-walker rewrite

Branch `intel-batched` (from `gpu-port-unified` 33a62652), session 2026-11-13 01:50-03:45.
"C" = measured (job id given), "E" = estimate or data-sheet number. Sources: `docs/gpu/scripts/intel_batched/`.
Earlier measurements used and not repeated: Phase 7 per-kernel profile (`docs/gpu/sessions/supermuc_large_prompt/phase2/team/p7/NOTES.md`, "p7"),
Phase 6 design study and gemm ceiling (`.../p6/DESIGN.md`, "p6").

<!-- STEP1-TABLES -->

## 2. Regime of each phase and upper bounds

### 2.1 What the hardware can do (E, data sheets; C where marked)

| | Max 1550, one tile | Max 1550 card (2 tiles) | H100 SXM (Helma) |
|--|--|--|--|
| FP64 peak (vector / matrix) | 26 TF/s | 52 TF/s | 34 / 67 TF/s (DMMA, used by cuBLAS ZGEMM) |
| HBM | 1.64 TB/s, 64 GB | 3.28 TB/s, 128 GB | 3.35 TB/s, 94 GB |
| last-level cache | 204 MB L2 | 408 MB | 50 MB |
| zgemm 576^3 measured | 12.8 alone, 10-11 node busy, 16-18.6 batched (C, p6 job 560781/787) | 2x | not measured here (E: 30-45 TF/s with DMMA) |
| zgemm 1024^3 measured | 18.5-20 alone, 14-15 node busy (C, p6) | 2x | |

So "by pure flops the Max 1550 should obliterate the H100" is **not** true for this code: in FP64 the H100's tensor-core rate (67) is above the whole card's
vector rate (52), the HBM bandwidth per device is the same (3.3 TB/s), and the card is two separate devices (tiles) that must each be filled by their own
work. Per *tile* the Max 1550 has 39 % of the H100's FP64 matrix rate and 49 % of its bandwidth. The two advantages of the Max 1550 are the large L2
(204 MB per tile: G of one walker, 5.3 MB per flavour, and the G of 16 walkers stay L2-resident) and the price/count of devices.

### 2.2 Required rate to match the H100 (L = 24)

Work per walker-sweep at L = 24: 9.2 TFLOP, of which 8.2 TFLOP in big zgemm (576^3), 0.49 TFLOP in rank-k flush gemms (34,560 calls, 366 GB of G traffic) (p6, E from operation counts).

| target | needed per tile | today per tile | needed gain |
|--|--|--|--|
| H100 single chain 3.51 s/sweep | one chain on one tile in 3.51 s (2.6 TF/s effective) | 4.04 s (C, V7 + LU-free) | 1.15x |
| H100 throughput 1.176 sweeps/s per GPU (8 ranks MPS) = 10.8 TF/s effective | card 1.176 -> tile 0.588 sweeps/s = 1.70 s tile time per walker-sweep = 5.4 TF/s effective | 0.42-0.44 sweeps/s (4 chains/tile, 9.3 s, C, older binary) = 2.3 s, 4.0 TF/s | 1.34-1.40x |
| H100 L = 32 throughput 0.232 sweeps/s per GPU | 0.116 sweeps/s per tile = 8.6 s tile time per walker-sweep (L32: ~69 TFLOP/sweep, E) = 8.0 TF/s | one chain/tile 16.4 s (C) = 0.061/tile; multi-chain L32 not measured | 1.9x vs one chain/tile |

The FLOP rates needed (5.4 TF/s per tile at L=24, 8 TF/s at L=32) are 40-55 % of what the tile's zgemm delivers at these sizes. **The H100 numbers are not
excluded by the hardware**: the gap is software (fill, latency), not peak rate or bandwidth. They are also not reachable by a small tweak: the throughput
target needs ~1.37x at L=24 and ~1.9x at L=32 relative to one chain per tile.

### 2.3 Regime per phase (L = 24, one chain = "L1", 4 chains/tile = "L32" in p7 notation)

Phase walls tonight (V7 + LU-free, one chain, C): updates 1.44, stabilisation 0.82, CGR 0.82, observables 0.74, rest 0.22 = 4.04 s.
Device busy from the p7 per-kernel profile (C, job 560831/845/847; CGR minus the removed getrf 0.38 s; observables minus the two serial corr kernels that
the observables time no longer contains: 1.18 -> 0.74 s).

| phase | wall (C) | device busy (E from p7) | regime, one chain | regime, 4 chains/tile |
|--|--|--|--|--|
| updates | 1.44 | 0.91 (rank-k 0.31, hopping zgemm 0.32, block small kernels 0.25) | **host-submission bound**: 112k async ops/sweep x 4-10 us submit = 0.88 s (C, p7), plus 8,640 D->H round trips of 22 us | submission contention (blk_end x1.7, C) + rank-k / hopping gemm saturating the tile (x3.3) |
| stabilisation | 0.82 | ~0.64 (device QR 0.41: panel = one work-group, 0.8 % fill; unmqr 0.13) | **under-fill / latency chain inside kernels** (QR panel) | hidden by the other chains |
| CGR | 0.82 | ~0.55-0.65 | same as stabilisation (QR) + small kernels 1x64 | hidden |
| observables | 0.74 | ~0.58 (corr 512x1024 kernels 0.31, zgemm 0.21) | throughput kernels with bad mapping (3.5x slower than their data volume) | **saturated** (2.4 of 2.73 s device demand, C+E, p7) |

None of the phases is HBM-bandwidth bound at the hardware limit: the most bandwidth-heavy op (rank-k flush, 10.6 MB per call) runs at 1.2 TB/s effective
on L2-resident data (p7: 9.0 us per call), which is below what the tile streams from L2 (Step 1 triad, table above).

### 2.4 Upper bounds (E)

| scenario | single chain L=24 | per card L=24 | basis |
|--|--|--|--|
| today | 4.04 s | 0.84-0.88 sweeps/s | C |
| current kernels, perfect asynchronous submission (no gaps, no host waits) | ~2.8-3.0 s (sum of device busy 2.7 + 0.2 rest) | ~0.95-1.0 (gaps are already filled by the other chains; only submission contention goes) | p7 busy times |
| + update block batched / fused (one launch per op for all walkers, one D->H per lockstep step) | — | +8-15 % (E: blk_end + blk_begin = 0.99 s/sweep host time at L1, x1.7 under sharing) | Step 1 block emulation |
| batched multi-walker code that fills the tile (all phases with a walker index, QR panels batched, batched gemm) | — | 1.3-1.6x of today = 1.1-1.4 sweeps/s (p6 grounded ceiling 1.5-2.3x) | p6 + this section |
| hardware floor (all FLOPs at measured zgemm rate, flush at L2 rate) | ~0.8 s | ~2.4 | 9.2 TFLOP / 12 TF/s per tile |

## 3. Verdict

<!-- VERDICT -->

## 4. Design of the rewrite (batched multi-walker, one process per tile or per compute slice)

Follows p6 section 1-5 (execution model, walker as last array index, per-walker RNG / files / MPI ids) with three changes that the evidence of this session
makes:
1. **Bitwise walker-b-vs-alone is not the acceptance test.** The owner's acceptance test is same-seed agreement with the same-compiler CPU build (Mz=F 1e-13,
   Mz=T 1e-10) and a planted defect that fails it. The batched oneMKL gemm differs from the single call at 1e-15 (p6 Gate 2 F: 2.2e-15 at L24F), two orders
   below the tolerance. So p6's no-go reason 4 ("batched gemm is not bitwise") does not block: batched gemm is allowed as long as the CPU same-seed gate
   passes (to be shown per kernel class, with the planted defect). Own element-wise kernels stay bitwise per walker by construction (shown by `ib_wtest`).
2. **The first target is the update phase**, because it is the only host-submission-bound phase at one chain and the one with submission contention at
   4 chains/tile; the multi-walker form removes (W-1)/W of the host submissions and D->H round trips.
3. **Layout:** W walkers in one process; `G_d(N,N,N_FL,W)`, per-walker strides passed to the SYCL kernels (`alfsb_*_w`); lockstep: all walkers are at the same
   (slice, block) so b, nd, nfa/nfb are shared; the pending count k differs per walker (kw(w) on the device, kmax for the launch range; the gemm uses the
   per-walker k with mode 0, or kmax with zero-padded columns with mode 1).

Kernels (alf_la_sycl_blk.cpp, done): `alfsb_gather_w`, `alfsb_scale_w`, `alfsb_fgather_w`, `alfsb_gemm_w` (mode 0: W oneMKL calls = bitwise; mode 1:
one `gemm_batch` strided).

Host threading: the Metropolis decisions of the W walkers of one block step are independent; at b = 32 sites x W walkers the host work per block is ~20 us
x W (p7); with W <= 8 run them in one host thread (loop over w), beyond that in `!$omp parallel do` over w (each walker's RNG state swapped in and out,
p6 section 4).

Switches (planned): `ALF_GPU_WALKERS=W` (default 1 = today's path, unchanged), `ALFSB_GEMM_BATCH=0|1` (mode of the flush gemm).

What stays: the single-walker path (W = 1) is today's code and remains the reference; QR / CGR / observables stay per walker in the first stage (called in a
loop over w on the same in-order queue, without host waits between walkers).

## 5. Implemented and measured so far

<!-- IMPLEMENTED -->

## 6. Handoff

<!-- HANDOFF -->
