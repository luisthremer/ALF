# CG-ALF ("Crazy good ALF") — plan (P2.9)

Started 2026-11-25 (evening), at the owner's request ("arrive at CG-ALF tonight and plan and run it").
Owner's goals (OWNER_BRIEF.md, point 15): largest system sizes, real GPU utilisation (FLOPs, not "busy"
time), FP32 where it is validated, a numbers module (compile-time real/complex/integer), multi-chain
throughput, new stabilisation and update algorithms; always model-agnostic, ALF stays modular.

## 0. Where we start (measured earlier, `gpu-port` docs)

| Fact | Source |
|---|---|
| H100, one chain, L = 24 / 32 (beta = L, submatrix, Nwrap 10): 3.51 / 9.82 s per sweep | `docs/gpu/platforms/helma_h100.md` |
| H100, 8 chains under MPS, L = 24: 6.80 s per sweep per chain = 1.18 sweeps/s per GPU | same |
| One chain per GPU is latency-bound (A30 L = 32: 44 % of kernel time in ~32 000 small flush GEMMs; 65 % of host API time in stream synchronisation) | `docs/gpu/SPEEDUPS.md` |
| Production (many chains under MPS) is throughput-bound: a kernel runs 99.97 % of the time — "busy" is not "efficient" | same |
| CPU profile (upstream): stabilisation 26-55 %, local updates 13-35 %, B products 6-42 % (dense hopping), measurements 5-17 % | `docs/gpu/PROFILE.md` |
| Production model: Hubbard, `N_FL = 2`, `Mz = .F.`, flux N_Phi = L^2/2, checkerboard: complex hopping and complex HS coupling | `/jakob_params` |

## 1. The metric

For every GPU (H100, A30, Max 1550, MI210) and size L:
- **useful FLOPs per sweep** F(L), counted (not estimated) per kernel class, with the algorithm as it is;
- **time per sweep** t(L), single chain and best multi-chain layout (sweeps/s per GPU);
- **efficiency** E = F / (t * P_peak), with P_peak the FP64 peak that the kernel class can use (tensor-core FP64
  for GEMM, vector FP64 otherwise), and the memory-bandwidth bound for O(N^2) kernels (roofline).
CG-ALF progress = sweeps/s per GPU at fixed physics accuracy, and E per kernel class. The "theoretical limit"
is t_min(L) = sum over kernel classes of max(F_k / P_k, bytes_k / BW).

## 2. Experiments, in order

| ID | What | Where | Output | State |
|---|---|---|---|---|
| E0 | FLOP and byte count per sweep per kernel class (counters in the linear-algebra layer, env switch, GPU path of `gpu-port`) | H100, L = 16/24/32 | table F_k(L), bytes_k(L) | planned |
| E1 | Library ceilings: ZGEMM/DGEMM/CGEMM/SGEMM (and batched) for N = 256 ... 4096 | H100 (cuBLAS), Max 1550 (oneMKL) | achieved TFLOP/s per type and size: the ceilings for every idea below | **tonight** |
| E2 | Timeline of one sweep: kernel time per class, gaps, host syncs (nsys), single chain and 8 chains | H100 | where the time goes on the GPU | tonight if budget allows |
| E3 | Efficiency E per class = E0 / E2, compared with E1 | — | the roofline table: which class is far from its ceiling | after E0, E2 |

## 3. Ideas, with the experiment that decides them

| # | Idea | Expected gain (to be measured) | Decided by | Numerics |
|---|---|---|---|---|
| I1 | **Walker batching**: many Markov chains in one process, every kernel batched over walkers (one launch for all chains) instead of MPS processes | removes launch latency and small-GEMM inefficiency; target: single-process throughput >= best MPS layout, then above | E1 (batched vs single GEMM), E2 | bitwise per chain (same operations) |
| I2 | **Real arithmetic** where the model allows it (real hopping, real HS coupling: Hubbard Mz without flux, many others): compile/run-time choice of real vs complex G | ~2-4x fewer FLOPs and 2x less memory for every N^3 kernel | E1 (DGEMM vs ZGEMM ceiling) | identical physics; roundoff differs |
| I3 | **Numbers module** (generalises I2; start from `precdef` / Florian Goth's work): kind parameters chosen at compile time per model class (real/complex, integer hopping where it pays) | enables I2, I4 cleanly | design review | — |
| I4 | **Mixed precision**: FP32 (or TF32/FP32-split on tensor cores) for wraps and updates between stabilisations, FP64 for stabilisation and measurements; precision checks decide Nwrap | up to 2x (FP32 vector) or more (tensor cores) on the GEMM share | E1 (SGEMM/CGEMM ceilings), then precision study: Precision Green, observables vs FP64 within error bars, ED tests | must be validated per model class |
| I5 | **Fewer stabilisations**: adaptive Nwrap from the measured precision; cheaper stabilisation (QR without pivoting + scale tracking, CholeskyQR2 was rejected earlier: re-check on H100) | stabilisation is 25-55 % of the time | E2, precision study | must keep Precision Green < 1e-8 |
| I6 | **Kernel fusion / CUDA graphs** for the per-slice chain (wrap, update, flush) | single-chain latency (the 65 % sync time) | E2 | bitwise |
| I7 | **Delayed/submatrix block size on GPU** retuned for H100 (N_delay larger, flush as one GEMM) | flush GEMMs are the largest kernel class at L = 32 | E1, E2 | roundoff |
| I8 | **Measurements**: equal-time observables as GEMMs over all sites (batched), skip unused outputs | 5-17 % | E2 | summation order |
| I9 | **Largest sizes**: memory per walker (G, UDV stack: Nstm * 3 N^2), recompute instead of store, FP32 storage of the stack | lets L = 64+ fit with many chains | memory model | stabilisation precision |
| I10 | **Throughput after warm-up**: split a thermalised chain into many independent chains (I1) | sweeps/s per GPU | — | statistical independence must be checked (autocorrelation) |

Model-agnostic rule: every idea is a property of ALF's matrices (real/complex, dense/checkerboard), never of one Hamiltonian.

## 4. Validation protocol (every idea)

1. Bitwise or roundoff-level agreement with the reference on the regression scenarios (`test_branch_parameters.json`) where the
   idea claims it; else
2. observables within error bars on long runs (energy, a correlation function), and `test_vs_ed` for all models;
3. Precision Green / Phase must stay below the reference's order of magnitude;
4. the timing claim measured on the target GPU, with the run directory kept.

## 5. Branches

R&D branches `rd/<idea>` (free history) on top of the GPU code (`gpu-port` today, `Code-modularity-GPU` after P2.4).
Anything that survives goes through dev/canvas (docs/ai/BRANCHING.md). Results go into this file (table in section 2/3)
and `docs/ai/phase2/CG_ALF_IDEAS.md` (one line per idea: status, measured gain).

## 6. Budget

GPU hours from the Helma test/porting contingent (3000 GPU-h total, shared with production): CG-ALF experiments at most
20 GPU-h per night unless the owner says otherwise; every job logged in `docs/ai/phase2/CG_ALF_LEDGER.md`.
