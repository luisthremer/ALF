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
| I5 | (owner idea 2026-11-26, alf-ops IDEAS.md; see task T1) **Fewer stabilisations**: adaptive Nwrap from the measured precision; cheaper stabilisation (QR without pivoting + scale tracking, CholeskyQR2 was rejected earlier: re-check on H100) | stabilisation is 25-55 % of the time | E2, precision study | must keep Precision Green < 1e-8 |
| I6 | **Kernel fusion / CUDA graphs** for the per-slice chain (wrap, update, flush) | single-chain latency (the 65 % sync time) | E2 | bitwise |
| I7 | **Delayed/submatrix block size on GPU** retuned for H100 (N_delay larger, flush as one GEMM) | flush GEMMs are the largest kernel class at L = 32 | E1, E2 | roundoff |
| I8 | **Measurements**: equal-time observables as GEMMs over all sites (batched), skip unused outputs | 5-17 % | E2 | summation order |
| I9 | **Largest sizes**: memory per walker (G, UDV stack: Nstm * 3 N^2), recompute instead of store, FP32 storage of the stack | lets L = 64+ fit with many chains | memory model | stabilisation precision |
| I10 | **Throughput after warm-up**: split a thermalised chain into many independent chains (I1) | sweeps/s per GPU | — | statistical independence must be checked (autocorrelation) |
| I11 | (owner idea 2026-11-26, inbox) **Cholesky-based stabilisation again** (CholeskyQR / CholeskyQR2 instead of Householder QR for the UDV decomposition; rejected earlier in the Intel sessions, gpu-port `docs/gpu/SPEEDUPS.md`, re-open on the H100). Owner: Nwrap probably has to change with it | stabilisation (QR panels 39 % of GPU kernel time, E2) becomes GEMM-like and cheaper | precision scan Precision Green vs Nwrap for QR and CholeskyQR2 at L = 24/32, then E-timing | CholeskyQR squares the condition number of the block, so the scale separation (D) must come before it; must keep Precision Green < 1e-8 |
| I12 | (owner idea 2026-11-26, inbox) **Time-first sweep**: run the local updates through imaginary time first instead of through space first; moving between time slices is then a wrap, i.e. matrix-matrix products (GEMM), which may suit the GPU better (not the CPU) | fewer small kernels, more GEMM work per update block | prototype on `rd/`: FLOPs and time per sweep on the H100 against the space-first order (E0/E2 tools) | any fixed sequential order keeps detailed balance; compare autocorrelation times, not only time per sweep |
| I13 | (owner idea 2026-11-26, inbox) **Optimised random path in imaginary space-time**: visit the (site, time) vertices in a random or optimised order instead of the fixed sweep | shorter autocorrelation per wall-clock time; may combine with I12 | autocorrelation time of energy and a correlation function per wall-clock second, against the fixed order | random sequential order keeps detailed balance; the stabilisation points (Nwrap) must still be visited in order |

Model-agnostic rule: every idea is a property of ALF's matrices (real/complex, dense/checkerboard), never of one Hamiltonian.

## 3a. Tasks from the owner's idea inbox (`manager_tools/IDEAS.md`, formerly alf-ops)

| # | Task | What decides it | Links |
|---|---|---|---|
| T1 | **Stabilisation literature survey** (owner, 2026-11-26): newer stabilisation algorithms for DQMC (e.g. scale-separated QR variants, block/communication-avoiding QR, mixed-precision stabilisation, log-scale schemes, fewer or adaptive stabilisations), then measure the candidates | a survey table (method, cost, accuracy, GPU fit, reference), then E-measurements on the H100: cost per stabilisation and Precision Green for each candidate vs the current QR with pivoting | idea I5; E2 shows QR panels = 39 % of GPU kernel time |
| T2 | **Real GPU usage** (owner, 2026-11-26): optimise FLOPs, not "busy" time | the efficiency E per kernel class (section 1) moving towards the E1 ceilings | section 1, E0-E4, ideas I1, I6, I7, I8 |
| T3 | **Fast `Op_V_is_symmetric`** (owner, 2026-11-26): revive branch `ai/feat/op-v-symm-fast` (reference patch) | start-up time at large N with `Symm = .T.` (KNOWN_ISSUES P1: hours today) and an unchanged result of the check on every test model | KNOWN_ISSUES P1; branch `ai/feat/op-v-symm-fast` |
| T4 | **Production guards from the Intel GPU code** (owner, 2026-11-26): review the production branches for features worth taking over, e.g. guards for production runs on HPC | a list of their features (guards, recovery, restart, monitoring) with a take / leave decision per feature, each taken one as its own dev/canvas feature | tags `stp/production-mpi-intel-large-prompt`, `stp/production-mpi-intel-large-prompt-team-combined`, `stp/production-mpi-intel-unstable-guard` |

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

## 5a. Tools for R&D sessions (owner, 2026-11-26)

- `research-desk` (citation check, paper notes and comparison, referee-style review): off by default in the ALF
  development sessions, maybe on for the R&D sessions (e.g. the stabilisation survey T1).
- `research-superpowers` (full literature workflow): disabled for now; maybe enable it later for R&D.
Both are account plugins (claude.ai settings); only the owner switches them.

## 6. Budget

GPU hours from the Helma test/porting contingent (3000 GPU-h total, shared with production): CG-ALF experiments at most
20 GPU-h per night unless the owner says otherwise; every job logged in `docs/ai/phase2/CG_ALF_LEDGER.md`.

## 7. Results

### E1 — library ceilings, H100 (Helma job 1452169, 2026-11-25), TFLOP/s in real flops

| N | batch | ZGEMM | DGEMM | CGEMM | SGEMM |
|---|---|---|---|---|---|
| 256 | 1 | 20.8 | 5.2 | 13.1 | 3.8 |
| 576 | 1 | 44.5 | 30.3 | 38.0 | 22.1 |
| 1024 | 1 | 57.0 | 52.8 | 52.0 | 38.1 |
| 2304 | 1 | 61.3 | 62.3 | 54.5 | 40.6 |
| 4096 | 1 | 60.9 | 63.4 | 54.1 | 51.4 |
| 576 | 8 | 58.0 | 47.5 | 52.8 | 38.2 |
| 1024 | 32 | 62.3 | 61.1 | 54.6 | 50.4 |

Rank-k flush ZGEMM (N x N x k): N = 576: k = 32 10.3 TFLOP/s (8.2 us per call), k = 128 26.7; N = 1024: k = 32
16.7, k = 128 35.4; N = 2304: k = 32 12.4, k = 128 36.5.

Consequences: I2 (real arithmetic) = up to 4x on the GEMM share (DGEMM as fast per flop as ZGEMM from N ~ 1000).
I4 (FP32) does not pay for GEMMs on the H100 (FP64 runs on tensor cores; CGEMM/SGEMM are slower per flop) — only for
memory-bound kernels and memory size; priority lowered. I1 (batching) pays at production sizes (N = 576: 44.5 -> 58).
I7: the flush at k = 32 runs at 10-17 TFLOP/s; k = 128 doubles it.

### E0 — FLOP count (branch `rd/flopcount`, H100, production parameters, 2 sweeps incl. setup, job 1452222)

| L (N) | GFLOP total | GEMM | QR (geqrf + ormqr + qrpp) | LU + TRSM | wall s |
|---|---|---|---|---|---|
| 16 (256) | 710 | 605 | 91 | 15 | 11 |
| 24 (576) | 12 098 | 10 303 | 1 548 | 247 | 30 |
| 32 (1024) | 35 557 | 22 125 | 11 585 | 1 847 | 284 |

### E2 — GPU timeline, L = 24, 2 sweeps + setup (nsys)

GPU kernel time 2.72 s of about 29 s wall (220 160 launches): QR panel / stabilisation 39 %, GEMM 26 %, LU + phase
(`getrf_pivot`, `k_cgr_phase` 0.74 ms per call) 17 %, update kernels 7 %, wraps 5 %. Host: 168 490 stream
synchronisations (2.8 s), 21 506 device-to-host copies.

### Reading

- Theoretical limit L = 24: ~5.5 TFLOP per sweep / 61 TFLOP/s = **0.09 s per sweep**; measured 3.5 s per sweep:
  a factor of about 40. The GPU is busy ~40 % of a sweep; of that busy time the panel factorisations (QR, LU) take
  more than half, and they run far below the GEMM ceiling.
- Start-up dominates short runs (L = 32: ~260 s of 284 s).
- First targets: (a) host round trips and syncs in the update (decisions on the host), (b) stabilisation without
  the LU (phase from the Householder scalars, done on Intel, not on CUDA) and a GPU-friendlier QR, (c) start-up cost,
  (d) batching chains (I1). Next experiment E3: setup vs sweep time and host hot spots.

### E3 — setup vs sweep, host profile (job 1452309; same inputs as E0, i.e. ham_h0 = 1)

| L | NSweep 2 | NSweep 6 (L = 32: 4) | s per sweep | setup |
|---|---|---|---|---|
| 24 | 29.7 s | 82.8 s | 13.3 | ~3 s |
| 32 | 283.8 s | 568.4 s | 142 | ~0 |

Correction to E0/E2: with `ham_h0 /= 0` the device observables of `gpu-port` are off (pinning field only on
`gpu-port-unified`), so the equal-time measurements ran on the host. Host profile (nsys CPU sampling, L = 24): `zslhemm`,
`zdslsymm` (host checkerboard products for the symmetrised G), `predefined_obs_eq_{den,spinmz,green}_measure`, complex
division. The GPU kernels take ~1.4 s of the 13.3 s per sweep. So the "setup" seen in E0 was not setup but host
measurements, and the docs' 3.5 s per sweep (device observables) is the right single-chain reference.
Lessons: (1) every measurement must stay on the device — a host fallback costs an order of magnitude, silently
(the info line says only that the GPU path is active); (2) E4 onwards uses ham_h0 = 0 (device observables).

### E4 — LU-free determinant phase on CUDA (`rd/cgalf-s1`, job 1453558), ham_h0 = 0 (device observables)

Time per sweep from NSweep 6 - NSweep 2 (H100, single chain, production parameters without the defect):

| L | before (`rd/flopcount`) | with S1 | change |
|---|---|---|---|
| 24 | 2.65 s | 2.63 s | ~0 (noise) |
| 32 | 8.57 s | 7.49 s | **-12.6 %** |

Correctness: confout byte-identical, acceptance identical, energy equal to 1e-16 relative; with ALF_GPU_CGR_LU=2 the
tracked phase agrees with the LU phase to 1.1e-14 over 194 CGR calls. With device observables the plain build runs at
2.65 s per sweep at L = 24 (E3 had 13.3 s with host observables). S1 is safe and pays from L = 32 on; it also removes
the LU product that caused the GPU determinant phase underflow. Candidate for the GPU line (owner decision).
