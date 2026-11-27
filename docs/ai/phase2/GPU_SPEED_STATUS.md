# GPU speed status: H100 vs a full CPU node, gap analysis, ranked levers

Session gpu-speed, 2026-11-27/28 (CG-ALF T2 + P2.9.5/P2.9.6). Measured on Helma.
Raw outputs: Helma `~/gpuspeed/` (`gs_*.o<job>`, `runs_*/*/bins.log`, `info`).

## Setup

- **Model:** Hubbard, square L×L, beta = L, dtau = 0.1, U = 5.3, N_FL = 2, N_SUN = 2, Nwrap = 10, Ltau = 0,
  submatrix updates. ham_h0 = 0 unless stated (see lever 1 for h0 = 1).
- **N_delay:** 32 at L ≤ 24, 64 at L ≥ 32 (the scan below).
- **CPU node:** 2× AMD EPYC 9965 (384 cores, 755 GB).
  - Code-modularity-CPU 84d38562, ifx 2025.3 + threaded MKL + Intel MPI, `HELMA MPI`.
  - Every MPI rank is an independent Markov chain, so node throughput = ranks / (s per sweep per rank).
- **GPU:** 1× H100 94 GB of a Helma GPU node (32 cores per GPU).
  - Branch `rd/cgalf-s1` (S1 merged), nvhpc 25.7, `ALF_GPU ALF_LA_CUDA`, noMPI build.
  - Multi-chain runs are independent processes, each with its own seed and directory, under MPS with 25 % active threads.
- **Timing:** steady-state s/sweep from the times at which bins complete (NSweep = 1 per bin, 1-s polling, the first bin
  excluded). Setup is not included; on the CPU, setup at L = 32 takes about 64 min (KNOWN_ISSUES P1, `Op_V_is_symmetric`).

## 1. Headline: one H100 vs one CPU node

| L | N | CPU node, best layout: s/sweep per chain | CPU node: sweeps/s | H100, 1 chain: s/sweep | H100, best MPS: sweeps/s | per-GPU / per-node | per-chain speedup |
|---|---|---|---|---|---|---|---|
| 16 | 256 | 17.0 (384×1) | 22.6 | GPU16 | GPU16M | R16 | C16 |
| 24 | 576 | 91.5 (192×2) | 2.10 | 2.68 | 1.73 (8 ch.) | 0.82 | 34 |
| 32 | 1024 | CPU32 | CPU32N | 7.04 | 0.40 (8 ch.) | R32 | C32 |
| 48 | 2304 | not run | | GPU48 | | | |

- **Per-chain speedup** = CPU s/sweep per rank ÷ H100 single-chain s/sweep. At L = 24 it is 34x against a rank with two
  threads, and 89x against a single-core rank (238.5 s at 384×1).
- **Per node:** a Helma GPU node has 4 H100, so at L = 24 it delivers about 6.9 sweeps/s against 2.1 for a CPU node, i.e.
  3.3x. Per card, one H100 at L = 24 does slightly *less* work than one 384-core EPYC 9965 node.
- **Other CPU layouts:**
  - L = 16: 192×2 gives 13.3 s/rank and 14.4 sweeps/s.
  - L = 24: 384×1 gives 238.5 s/rank and 1.61 sweeps/s (289 GB RSS).
  - The h0 = 1 pinning field costs the CPU about 15 % (L = 16, 384×1: 19.6 s against 17.0 s).
- **Other vendors**, from their own reports and not re-measured here. Those runs used h0 = 1 and other branches, so they are
  not comparable 1:1:
  - MI210: L = 24 at 4.41 s/sweep single chain; 4 concurrent chains at 8.9 s each, i.e. 0.45 sweeps/s.
  - Max 1550 (`intel-batched`): L = 24 at 1.128 sweeps/s per card, i.e. 0.96x of an H100 measured the same way;
    L = 32 at 0.194 sweeps/s per card (0.84x).
- **Earlier single-core comparison** (`ANTRAG_SUPERMUC/common/data/gpu_vs_cpu.dat`, EPYC 7763 core at L = 24): 148.6 s/sweep,
  so an H100 is 42x faster than *one core*. The full-node baseline here replaces that figure for node-to-node statements.

## 2. Gap analysis (H100, single chain, nsys, h0 = 0)

Kernel time by class: three sweeps plus setup, `rd/cgalf-s1`.

| class | L = 24: s | L = 24: % | L = 32: % | where it stands vs E1 |
|---|---|---|---|---|
| QR panel (geqr2/larft/larfb) | 1.43 | 36.4 | 32.5 | ~0.1 TFLOP/s effective (est.) vs 44.5 TF ZGEMM: about 300x below |
| GEMM (cuBLAS) | 1.21 | 30.7 | 30.1 | 150k launches at 8 µs average: small or skinny GEMMs, k = 32 flush at 10-17 TF |
| stabilisation reductions (`k_cgr_phase`, `k_prec`, `k_colnorms`, `k_phase`) | 0.46 | 11.7 | 9.2 | one-thread serial loops, 0.74 ms per call at N = 576 |
| lowrank update kernels | 0.27 | 6.9 | 4.6 | latency-bound, 104k launches |
| observables | 0.23 | 5.9 | 5.7 | |
| wrap (checkerboard) | 0.18 | 4.7 | 14.5 | `k_cb_left_g` 58k × 24.5 µs at L = 32 |
| LU/TRSM + other | 0.15 | 3.7 | 3.4 | |
| **total kernel** | **3.93** | | 11.84 s | |

- **Host side at L = 24:**
  - 265k `cuStreamSynchronize`, 68 % of CUDA API time.
  - 118k H2D and 34k D2H copies.
  - 338k kernel launches.
  - L = 32: 380k syncs and 548k launches.
- **The GPU is busy only about 50 % of a single-chain sweep:** kernel time is about 1.3 s/sweep at L = 24 against 2.68 s
  wall, and about 3.9 s against 7.04 s at L = 32.
  - Main cause: one host round trip per vertex step. `gpu_lowrank_mod` `k_sub_up`/`k_sub_down` do a
    `target update from` so that the host can make the Metropolis decision. That is on the order of 10⁵ syncs per sweep.
  - MPS hides part of this idle time, because other chains fill the gaps. At L = 24, 8 chains give 4.6x the single-chain
    throughput; at L = 32 the gain is only 2.8x, because kernel time dominates there.
- **Why QR dominates although cgalf-research T1 found CholeskyQR2 worth only 1-3 %:** T1 replaced only the decompose QR.
  The CGR QR (`la_qrpp`, pivoted column norms on the host, then `geqrf`) carries most of the panel time.
  (Inferred from the call counts; not yet checked call by call.)

## 3. Ranked levers

Gain is in sweeps/s per GPU in the production layout (MPS multi-chain) unless stated.

| # | lever | expected gain | evidence | cost | line / phase |
|---|---|---|---|---|---|
| 1 | Device observables for ham_h0 ≠ 0 (the production pinning field, Jakob's parameter set). Today `gpu_obs_mod.F90:111` returns early and every measurement runs on the host Obser. | L = 24, h0 = 1: 13.3 → ~2.7 s/sweep, i.e. up to 5x for production runs | measured: E3 13.3 s against 2.6 s at h0 = 0; H0SLOW | small to medium: one more local observable (Obs_local) on the device | GPU line, P2.9.x |
| 2 | Nwrap 10 → 20 (cgalf-research T1, `rd/cgalf-stab` 7ad95cf6, `docs/cgalf_T1/README.md`) | −20 % sweep time at L = 24/32 | measured by T1; Precision Green max ~1e-10 (cold start; a thermalised check is running) | parameter only, after the thermalised check | CG-ALF T1 / production input |
| 3 | Remove the per-vertex host round trip: Metropolis decision on the device (pre-generated uniforms keep the chain bitwise identical), then CUDA graphs over a whole N_delay block | single chain up to −45 % (GPU idle about 50 %); in MPS production est. −10 to −25 % | estimated from nsys (85k syncs per sweep, 1.4 s idle per sweep at L = 24) | medium: `gpu_lowrank_mod`, RNG order | GPU line, P2.9.5 |
| 4 | Faster CGR QR: CholeskyQR2 or batched flavours (2 QRs per call) for the `la_qrpp` path, not only the decompose | QR is 33-36 % of kernel time; est. −10 to −20 % kernel time | measured share; gain estimated from T1's 16-26 % per QR plus batching | medium; Precision Green gate as in T1 | CG-ALF (extends T1) |
| 5 | Real arithmetic for real-HS Hubbard (I2 / P2.9.8) | est. −30 to −50 % kernel time (4x fewer flops in GEMM and QR) | estimated from E1: DGEMM per flop ≈ ZGEMM from N ≈ 1000 | high: type split through the stack | P2.9.8 |
| 6 | Parallel per-i factors in the stabilisation phase kernels (this session's prototype) | single chain −3 to −6 % at L = 24, 0 to −2 % at L = 32; MPS ±1 % | measured, see below | done, prototype | GPU line, P2.9.6 |
| — | Larger N_delay | **negative**: at L = 32 it is slower above 64 | measured, scan below | — | drop |
| — | N_FL = 1 when h0 = 0 (spin symmetry) | −50 % | estimated | model-specific; production uses h0 = 1 | not pursued |

**Ranking logic:** 1 applies only to h0 ≠ 0, but that is the production case. 2 is a free parameter change.
3 and 4 are the two code levers that act on the measured gap. 5 is the largest but also the most expensive.

### Prototype (lever 6)

Branch `rd/gpu-speed-proto` (from `rd/cgalf-s1`), commit bd3543af. It changes only `Prog/gpu_stab_mod.F90`:
- `k_phase` and `k_cgr_phase` compute their per-i factors in a parallel kernel, and only the product stays serial.
- `k_prec` reads row-wise.

Helma job 1471267, single chain, 2 repetitions each; 8-chain MPS A/B:

| | S1: s/sweep | proto: s/sweep | change |
|---|---|---|---|
| L = 24, 1 chain | 2.683 / 2.683 | 2.515 / 2.587 | −3.6 to −6.3 % |
| L = 32, 1 chain | 7.038 / 7.038 | 7.038 / 6.894 | 0 to −2.0 % |
| L = 24, 8 chains (per chain) | 4.640 | 4.616 | −0.5 % |
| L = 32, 8 chains (per chain) | 20.23 | 19.91 | −1.6 % |

- **Kernels (nsys, L = 24):**
  - stabilisation reductions 0.461 s → 0.381 s (−17 %);
  - `k_cgr_phase` 738 µs → 445 µs;
  - `k_prec` became slower, 249 µs → 309 µs, so its row-wise read should be reverted.
- **Correctness:** `confout_0` is bitwise identical, so the Markov chain is unchanged. `Ener_scal` differs only in the last
  bit (relative 3e-16), because the phase product is reordered.
- **Conclusion:** a small but real single-chain gain. Moving the serial product to the host (copying N numbers) should
  remove most of the remaining 445 µs; not tried.

### N_delay scan (H100, single chain, h0 = 0)

NDSCAN

## 4. Hours used

| job | partition | what | elapsed | hours |
|---|---|---|---|---|
| 1463685 gs_gpu | preempt H100 | first pass (multi-chain runs invalid: identical seeds) | 9.4 min | 0.16 GPU-h |
| 1463710 gs_gpu2 | preempt H100 | N_delay scan, coarse | 7.7 min | 0.13 GPU-h |
| 1471249 gs_gpu3 | preempt H100 | failed at start (taskset core ids) | 0.3 min | 0.005 GPU-h |
| 1471267 gs_proto | preempt H100 | prototype A/B + nsys | 10.4 min | 0.17 GPU-h |
| GPU4JOB | preempt H100 | single chain, N_delay, MPS layouts | GPU4T | GPU4H |
| 1463702 gs_cpuA | preempt_cpu | L16, L24 at 384×1 | 26 min | 0.44 node-h |
| 1463703 gs_cpuB | preempt_cpu | L24/L32 at 192×2 (L32: no bin after 64 min of setup) | 98 min | 1.64 node-h |
| 1471250 gs_cpuC | preempt_cpu | L16 with h0 = 0, L24 at 192×2 | 15 min | 0.26 node-h |
| CPUDJOB | preempt_cpu | L32 at 96×4 | CPUDT | CPUDH |
