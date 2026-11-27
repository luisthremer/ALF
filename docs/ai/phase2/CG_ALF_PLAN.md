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
| I1 | **Walker batching** (base: task T12, finish the walker): many Markov chains in one process, every kernel batched over walkers (one launch for all chains) instead of MPS processes | removes launch latency and small-GEMM inefficiency; target: single-process throughput >= best MPS layout, then above | E1 (batched vs single GEMM), E2 | bitwise per chain (same operations) |
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
| T3 | **Fast `Op_V_is_symmetric`** (owner, 2026-11-26): revive branch `ai/feat/op-v-symm-fast` (reference patch) | start-up time at large N with `Symm = .T.` (KNOWN_ISSUES P1: hours today) and an unchanged result of the check on every test model | KNOWN_ISSUES P1; branch `ai/feat/op-v-symm-fast`; done: `cm/feat/opv-symm-fast` (LEVERAGES row 1) |
| T4 | **Production guards from the Intel GPU code** (owner, 2026-11-26): review the production branches for features worth taking over, e.g. guards for production runs on HPC | a list of their features (guards, recovery, restart, monitoring) with a take / leave decision per feature, each taken one as its own dev/canvas feature | tags `stp/production-mpi-intel-large-prompt`, `stp/production-mpi-intel-large-prompt-team-combined`, `stp/production-mpi-intel-unstable-guard` |
| T5 | **Is the tracked determinant phase right in general?** (owner, 2026-11-26, with the yes to S1) Check all stabilisation variants, L > 32, checkpoint/restart, tempering and permutations keep the invariant (U = product of the tracked reflectors), models with a sign problem; then the **LU-free phase on the CPU line** too ("ich will das für cg alf") | check mode (LU vs tracked) on every case (on the GPU line: `ALF_GPU_DEBUG=1 ALF_GPU_CGR_LU=2`, step 1f; the invariant is written in the header of `gpu_stab_mod`); CPU-line timing gain | P2_4_GPU_PLAN step 1f; E4 |
| T6 | **`Lowrank_cross_slice`** (owner, 2026-11-26): research branch now, in parallel within the WIP limit. Is it mathematically sound (detailed balance, exactness of the cross-slice update)? Is it the fastest scheme on CPU and GPU? | derivation note + timing against delayed/submatrix on CPU and H100 at L = 24/32, same observables within error bars | gpu-port `Lowrank_cross_slice`; ROADMAP P2.6 (all models first) |
| T7 | **Intel ideas on CUDA/AMD** (owner, 2026-11-26): fewer submissions / block fusion, next-block prefetch, ordered observable reads (`ALFSH_*` of `intel-batched`) | time per sweep on Helma H100 (and MI210) with and without | gpu-port `docs/gpu/platforms/INTEL_BATCHED.md` |
| T8 | **Compile-time development timing flag** (owner): `#ifdef` timers, production build otherwise unchanged; revive gpu-port `Timers_mod` / `ALF_TEAM_TIMERS` as base | per-phase times match nsys within a few %; zero cost when off (byte-identical build without the flag) | ROADMAP CPU-line queue item 4 |
| T9 | **Fixed benchmark list + tuning script** (owner): models x L x hardware, tuning script with a correctness gate, results table | the table, reproducible from the script | T8 |
| T10 | **Start-up auto-selection of tuning settings** (owner: "if it is worth it") | gain of auto-selected vs fixed defaults on the T9 list larger than the selection cost | T9 |
| T11 | **New hardware** (owner: "we have to think about a way"): mechanism so tuned values (e.g. `ALFSH_*`, tuned on Max 1550 only) are re-checked on new GPUs: warn/fail when the device is not in the tuned table, or parameters-file entries — to be designed, not decided | a design note for the owner | T9, T10 |
| T12 | (owner, 2026-11-26) **Finish the walker**: Fakher's `qmc_walker_t` (Phase 6) holds G and the UDV stacks, but `main` runs one walker and per-chain state still lives in module variables (e.g. pending low-rank state, control counters, GPU device state). Move all per-chain state into the walker so one process runs N walkers; then batch the GPU kernels over walkers (idea I1) to raise throughput | N walkers in one process give the same chains as N separate runs (bitwise per chain); then sweeps/s per GPU against the MPS layout (H100: 1.18 sweeps/s at L = 24 with 8 processes) | `Prog/qmc_walker_mod.F90`, `docs/ai/codemap/engine.md`; I1, E1 (batched GEMM N = 576: 44.5 -> 58 TFLOP/s) |

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

### T6 — `Lowrank_cross_slice` (2026-11-26, session rd-cgalf; code on `gpu-port`, `Prog/Lowrank_update_mod.F90`)

**What it does.** Delayed/submatrix updates keep G in factored form $G = G_0 - X Y^T$ (submatrix: $G = G_0 - G_0(:,S)\,M\,(E_S^T - G_0(S,:))$)
and normally flush ($G_0 \leftarrow G_0 - XY^T$, one rank-$k$ GEMM) at every slice end. With `Lowrank_cross_slice = .T.`
the flush at the slice end is skipped; the pending factors are carried through the hopping wrap.

**Derivation (exactness).** The wrap between slices is a similarity transform, $G \mapsto B\,G\,B^{-1}$ with $B = e^{-\Delta\tau T}$
(up) or $B^{-1}$ (down). It is linear in G, so
$$B\,(G_0 - XY^T)\,B^{-1} = B G_0 B^{-1} - (BX)\,(Y^T B^{-1}).$$
Hence the update $G_0 \leftarrow BG_0B^{-1}$ (ALF's wrap, unchanged), $X \leftarrow BX$, $Y^T \leftarrow Y^T B^{-1}$
(`Hop_mod_mmthr`, `Hop_mod_mmthl_m1`) represents exactly the same G. A pending submatrix block is first converted to
delayed form, $X = G_0(:,S)M$, $Y_{:,m} = e_{S_m} - G_0(S_m,:)^T$, an exact identity. The vertex-wrap hook inside a slice
is the same argument with $B = e^{V}$ already used by the non-cross-slice code. So in exact arithmetic G is the same
matrix at every point as with the slice-end flush (or the standard rank-1 scheme).

**Detailed balance.** The Markov chain is defined only by (i) the order of proposals, (ii) the acceptance ratio
$R = \det[\,1 + \Delta (1 - G_{PP})\,]$ (times the bosonic weight ratio) and (iii) the random numbers. Cross-slice changes none of
these: the order is the same fixed sequential sweep, and $G_{PP}$ is read from the exact low-rank representation
(`Lowrank_propose` adds the $XY^T$ correction on the $P$ rows/columns). So the transition kernel is identical to
the standard scheme's, which satisfies detailed balance (Metropolis/heat bath per vertex, fixed sequential order:
stationarity holds per step, i.e. "balance"; ALF's usual argument). Every other reader of G (CGR / stabilisation,
measurements, global/Langevin/HMC/tempering moves) is preceded by `Lowrank_read_point` → flush (13 read points in
`main.F90`); `Global_tau` is refused. Therefore the method is **sound**; the only difference is roundoff.

**Roundoff.** The pending factors now pass through up to Nwrap hopping wraps before the flush, and $\|B\|\,\|B^{-1}\|$ grows
like $e^{2\Delta\tau\|T\|}$ per slice; the error of the delayed form grows with the number of slices carried, not
with the matrix size. Measured earlier (2026-10-31, CPU): trajectories identical to the key-off run on start/tv/kondo/
tv_flux/pv_adiab (20/20), Kondo Precision Green ratio on/off 1.47; debug shadow error ≤ 1.6e-12. Acceptable, but it
is a loss, not a gain, in precision.

**Cost model (why it gains little).** Per slice it saves one flush (ZGEMM $N\times N\times k$, $8N^2k$ flops) only when the
slice ends with $k < N_\text{delay}$ pending columns; the capacity flushes inside the slice are unchanged. It adds two
checkerboard hops on $N\times k$ and $k\times N$ blocks per slice ($O(N_\text{bond} k)$, cheap) plus the conversion
for submatrix ($O(Nk^2)$). For Hubbard the number of accepts per slice is about acceptance × N ≈ 0.4 N, i.e. 26 / 58 / 102
for L = 8 / 12 / 16 with N_delay = 32: the slice-end flush is at most one of 1-4 flushes per slice, and its saving is
partly the same GEMM done one slice later with larger k. Expected gain: a few % at most at small L, → 0 at large L.
It could matter only where few accepts per slice meet a large N_delay (low acceptance, very large N_delay on the GPU,
I7), and there it trades against carrying error over more wraps.

**CPU timing (Mac, gfortran -O3, 1 thread per run, `gpu-port` 078d806b built on `rd/cgalf-t6`).** Hubbard Mz, checkerboard,
U = 4, Beta = 4, Dtau = 0.1 (40 slices), N_delay = 32, same seed (all schemes run the identical trajectory: acceptance
equal to all digits), 2 repetitions, s per sweep (mean; range in brackets). Raw: `docs/ai/phase2/cgalf_t6_raw.txt`
(columns as `testsuite/lowrank/bench.sh`). L = 8/12 ran 3 jobs in parallel at load ~12 (noisy, ranges up to 12 %); L = 16
ran 2 in parallel at load ~5.

| L (N) | standard | delayed off | delayed on | on/off | submatrix off | submatrix on | on/off |
|---|---|---|---|---|---|---|---|
| 8 (64) | 0.052 | 0.044 (±5 %) | 0.047 (±4 %) | 1.07 | 0.056 (±4 %) | 0.053 (±3 %) | 0.94 |
| 12 (144) | 0.424 | 0.264 (±5 %) | 0.289 (±6 %) | 1.10 | 0.301 (±4 %) | 0.300 (±2 %) | 1.00 |
| 16 (256) | 2.22 | 1.455 (±1.4 %) | 1.477 (±1.2 %) | 1.015 | 1.446 (±0.1 %) | 1.468 (±0.4 %) | 1.015 |

Earlier scan (2026-10-31, Beta = 2, L = 8/12, 3 reps, quiet machine): on/off 1.074 / 0.967 / 1.016 / 1.017.

**T6 verdict.** Sound (exact up to roundoff, same Markov chain; derivation above) but **not faster on the CPU**: at L = 16,
the only low-noise point, cross-slice is 1.5 % *slower* for both schemes; at L = 8/12 the differences are inside the
noise. Delayed/submatrix themselves give 1.5x over standard at L = 16. This matches the cost model (the slice-end flush
is 1 of several flushes per slice and its saving shrinks with N). On the GPU the only plausible win is fewer small flush
kernels per slice when N_delay is much larger than the accepts per slice (I7); that is better reached by tuning N_delay.
Recommendation: do not port it to the GPU (PARKING #39 stays parked); keep it opt-in on the CPU or remove it (owner's open
question 00 of 2026-10-31). GPU timing not done (needs Helma and the owner's go; not worth it given the CPU result).
Repeat the L = 16 row at low load before quoting it outside this file.

### T1 — stabilisation literature survey (2026-11-26, desk survey; nothing measured yet)

References are from memory, not checked against the papers in this session (research-desk was off): verify before
citing. "Cost" is per stabilisation of one N×N block (N = Ndim), compared with ALF's current QR with column pivoting (QRP,
`ZGEQP3`, $\tfrac{4}{3}N^3$ real-flop equivalents ×4 for complex, BLAS-2 bound in the pivoting).

| Method | Idea | Cost vs QRP | Accuracy | GPU fit | Reference |
|---|---|---|---|---|---|
| UDV with QRP (ALF today, STAB1/2/3) | $B\cdots B = UDV$, D separated by pivoted QR every Nwrap slices | 1 | reference; stable for any β | poor: pivoting is a column-norm search per step (BLAS-2, host syncs); 39 % of H100 kernel time (E2) | Loh et al., PRB 1989 / Int. J. Mod. Phys. B 1992; Assaad & Evertz, LNP 739 (2008); ALF 2.0, SciPost Phys. Codebases (2022) |
| Pre-pivoted QR | pivot once by the column norms of $B\,U D$ (sort, permute), then *unpivoted* blocked QR (`geqrf`) | ~0.5-0.7 (all BLAS-3) | same as QRP in practice (D is already graded, so norm sorting ≈ QRP pivots) | good: `geqrf` is the fastest QR in cuBLAS/cuSOLVER/oneMKL | Bai, Lee, Li, Xu, Linear Algebra Appl. 435, 659 (2011); used by QUEST |
| Scale-separated / log-scale D | store D as log or (mantissa, exponent); build G from $D_{\max}, D_{\min}$ split (STAB3, `STABLOG`) | ≈1 (+O(N)) | removes over/underflow for large β·U; needed at L ≥ 24 production | neutral (O(N) work) | Bauer, SciPost Phys. Core 2, 011 (2020) (`StableDQMC.jl`); ALF STAB3 / gpu-port `STABLOG` |
| Unpivoted QR + scale tracking | QR without any pivoting, D from $|R_{ii}|$, rely on the graded structure | ~0.5 | can lose digits when the column scales are not ordered; ALF test 108 (gpu-port `docs/gpu/reports/108-stablog-pivot.md`) | very good | Bauer 2020 (comparison of QR variants); ALF report 108 |
| CholeskyQR2 / shifted CholeskyQR3 | $A^\dagger A = R^\dagger R$, $Q = AR^{-1}$, repeated | ~0.3-0.5, all GEMM/SYRK | needs κ(block) ≲ 1e8 (CQR2) or ≲ 1e16 (sCQR3): only *after* D is factored out (apply to $B\,U$ columns scaled by D) | excellent (GEMM ceiling, E1: 60 TFLOP/s) | Fukaya, Nakatsukasa, Yanagisawa, Yamamoto, ScalA 2014; Fukaya et al., SIAM J. Sci. Comput. 42, A477 (2020) — idea I11 |
| TSQR / CAQR (communication-avoiding) | tree of small QRs on row blocks; reduces panel latency | ~1 flops, fewer syncs | as Householder QR (stable) | good for tall-skinny panels; for square N×N mainly fewer kernel launches; batched over walkers (I1) | Demmel, Grigori, Hoemmen, Langou, SIAM J. Sci. Comput. 34, A206 (2012) |
| Randomized QRCP (HQRRP) | choose pivots of a whole block from a small sketch $\Omega A$, then blocked QR | ~0.5-0.6, BLAS-3 | pivots as good as QRP for rank revealing | good | Martinsson, Quintana-Ortí, Heavner, van de Geijn, SIAM J. Sci. Comput. 39, C96 (2017) |
| Mixed precision | wraps/updates in FP32 (or TF32), stabilisation and G in FP64; or FP32 panel + FP64 refinement | 0.5 on the GEMM share (FP32 vector) | must keep Precision Green < 1e-8; Nwrap shrinks | H100: no gain for GEMM (E1: FP64 tensor cores faster than FP32), gain on memory-bound parts only; A30/MI210/Max 1550 differ | Yang, Fox, Sanders, SIAM J. Sci. Comput. 43, A1723 (2021) (mixed-precision block Householder QR); idea I4 |
| Adaptive Nwrap | choose Nwrap at run time from the measured Precision Green (stabilisation only as often as needed) | 1/(Nwrap gain) | controlled by construction | neutral; fewer syncs | no dedicated paper found from memory; ALF already reports Precision Green; idea I5 |
| SVD-based (`gesvj`, Jacobi) | UDV via SVD instead of QR | 3-10 | most accurate | poor | Bauer 2020 (shows QR suffices) |

**Reading.** The best GPU candidates are (1) **pre-pivoted unpivoted QR** (same accuracy as QRP in the literature, all BLAS-3,
no code change outside the stabilisation kernel), (2) **CholeskyQR2 after scale separation** (GEMM speed; owner's I11), with
(3) **adaptive Nwrap** on top of either. Mixed precision is low priority on the H100 (E1). Next step: E-measurement on the
H100 (Helma, owner's go needed): cost per stabilisation and Precision Green vs Nwrap for QRP, pre-pivoted QR and CholeskyQR2
at L = 24/32 on the production model.
