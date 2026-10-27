# What makes the ALF GPU port hard

Ranked by expected impact on the port, as of the end of Phase 0 (2026-10-27). Each item gives
the evidence (code with file:line, or a measurement) and the proposed strategy. Measurements:
[PROFILE.md](PROFILE.md) (CPU time per kernel), [BASELINE.md](BASELINE.md) (CPU references,
validation lessons), [MACHINES.md](MACHINES.md). Numbers marked *[spec]* are vendor data sheet
values; *[est]* marks back-of-envelope estimates that Phase 1 must measure.

Notation: N = number of sites (orbitals per flavour), L = time slices, F = flavours (2 here),
Nw = Nwrap, a = acceptance (≈0.43 in our runs).

## The bar a GPU has to clear

ALF's CPU parallelism is embarrassingly parallel: one independent Markov chain per MPI rank,
one rank per core. A GPU therefore does not compete with one core but with **all the cores it
replaces**: 16 Zen 3 cores on nividia (one A30), and 112/8 = 14 Sapphire Rapids cores per PVC stack
on SuperMUC-NG Phase 2. To win, one GPU chain must finish a sweep at least ~15× faster than a
single-core chain, or the GPU must run several chains at once. All the items below should be read
against that bar. Single-core s/sweep (PROFILE.md): 0.6 s at N=144, 8 s at N=400, 193 s at N=1024
(checkerboard, L=40).

## 1. Sequential rank-1 local updates (highest impact)

**What:** every site and time slice proposes a field flip. Its acceptance ratio needs the current
G, and an accepted flip changes all of G at once (`ZGERU`, `Prog/upgrade_mod.F90:269`; accept
test against `ranf_wrap()` at `:222`). The chain is strictly sequential: 2·L·N decisions per sweep,
each depending on the previous update.

**Evidence:**
- Rank-1 update share of CPU time: 17 → 24 → 33 % (checkerboard) and 10 → 15 → 20 % (dense)
  for N = 144 → 400 → 1024. It is the largest kernel at large N and its share keeps growing,
  because it is memory-bound BLAS-2 while the rest is mostly BLAS-3 (PROFILE.md).
- Bandwidth floor on the A30 *[est]*:
  - at N=1024 one ZGERU reads and writes 2 × 16.8 MB of G, about 36 µs at 933 GB/s *[spec]*;
  - 2·L·N·F·a ≈ 70 000 accepted flips per sweep, so about 2.5 s/sweep;
  - the 16-core node does the same work in 64 s / 16 = 4 s per chain-sweep.
  So even a perfect rank-1 kernel only roughly matches the node.
- Latency floor *[est]*: each decision is a device→host→device round trip if the accept test stays
  on the host, as the RNG stream requires (item 7). At ~10 µs per round trip:
  - N=1024: 81 920 decisions ≈ 0.8 s/sweep;
  - N=144: 11 520 decisions ≈ 0.12 s/sweep, i.e. 3× *slower* than the node's 0.038 s per
    chain-sweep, from latency alone.
  R3 §5 expects higher launch latency on PVC than on CUDA.

**Strategy:**
- Phase 2 keeps rank-1 updates, as the brief requires, so results can be checked against the CPU
  code flip by flip. The ratio and the update run on the device; per site, only a few scalars
  (ratio inputs) go to the host, and the accept/reject decision stays on the host.
- The update code is structured around a "flush" interface so that delayed or submatrix updates
  can replace it later.
- Delayed or submatrix updates (Sun & Xu) are the real fix. They turn 2·L·N BLAS-2 updates into
  ~2·L·N/n_d GEMMs, and with submatrix updates the host only needs a small k×k matrix (Γ) per block.
  That split is natural for a GPU: scalar Markov-chain logic on the host on small matrices, large
  GEMMs on the device. This becomes the main performance lever after Phase 2. **Phase 1 measures
  rank-1 vs rank-k on the device to size this.**

## 2. Pivoted QR in every stabilization, including STABLOG

**What:** the UDV decomposition, CGR (equal-time G) and CGR2_2 (2N×2N, for τ) all call
`QDRP_decompose`, which is `ZGEQP3`, i.e. QR with column pivoting (`Prog/QDRP_decompose_mod.F90:78-82`).
- In `STABLOG`, only the UDV decomposition pre-pivots by column norm and disables LAPACK
  pivoting (`IPVT = 1`, `Prog/udv_state_mod.F90:503-542`).
- **CGR (`Prog/cgr1_mod.F90:305-312`) and CGR2_2 (`Prog/cgr2_2_mod.F90:346, 405, 420`) set
  `IPVT = 0` in every stabilization variant, so they do full QRCP under STABLOG too.** The
  research reports and the brief assumed STABLOG already avoids pivoted QR everywhere; it doesn't
  (logged in DISCREPANCIES.md).

**Evidence:**
- Stabilization takes 42–55 % (checkerboard) and 26–32 % (dense) of the CPU time. CGR2_2 alone is
  14–29 %, the largest stabilization kernel. CGR is 7–15 % (PROFILE.md).
- No GPU QRCP exists in cuSOLVER or in oneMKL (OpenMP dispatch or SYCL) (R2 §2, R3 §3, to be
  confirmed against the library docs in Phase 1). MAGMA has `zgeqp3_gpu` for CUDA/HIP only, and
  it is BLAS-2-bound.
- A host fallback is expensive *[est]*: one CGR2_2 at N=1024 takes ≈5.7 s on one core (45 s/sweep
  over 8 calls). Multithreaded host LAPACK helps, but 4N² complex values must cross PCIe each way
  per call.
- oneMKL's Fortran `!$omp dispatch` runs `geqrf`, `unmqr` and `getrf` partly on the host, and
  `ungqr` entirely on the host (R3 §3), so the Intel path needs a SYCL shim or our own kernels.

**Strategy:**
- On the GPU, replace QRCP with **pre-pivoting by column norm followed by unpivoted `geqrf`**
  (+ `ormqr`/`ungqr`) in all three places, as STABLOG's decompose already does.
- This changes the numerics of CGR/CGR2_2 relative to upstream. It must be validated **on the CPU
  first** (Phase 1 ticket): `Control_PrecisionG` / `Precision tau` and the observables against QRCP
  at increasing β (up to a strongly ill-conditioned regime), for both the default and LOG scales.
- If pre-pivoting is not accurate enough, the fallbacks in order are:
  1. MAGMA on NVIDIA;
  2. a hybrid that computes the column norms and pivot order on the device, then runs blocked
     QR (this is pre-pivoting again);
  3. host QRCP.
- `GPU=none` keeps the upstream QRCP unchanged.

## 3. Small matrices: the GPU is under-filled and launch-bound

**What:** N is 100–1000 in typical runs (ALF's own docs say the matrices are small:
`Documentation/conclusion.tex`). An N=144 complex matrix is 330 KB and a ZGEMM on it is ~24 MFLOP
(≈5 µs at the A30's 5.2 TFLOP/s FP64 *[spec]*, comparable to a launch), far too little to fill
a GPU (≈1800 FP64 units on the A30) or a PVC stack. Every routine loops over
flavours on the host (e.g. `Prog/wrapgr_mod.F90:110-117`), and nothing is batched.

**Evidence:**
- PROFILE.md: at N=144, 31 % of the time is in ALF's own Fortran loops and 35 % in BLAS-1/2;
  only 23 % is BLAS-3.
- Per sweep: 46 080 `Op_Wrapup/Wrapdo` calls, 20 M `ZDSLSYMM` calls (2×2 checkerboard bonds),
  138 k small allocations (`fields_make`) (12×12 call counts).

**Strategy:**
- Fuse instead of calling per element: one kernel per checkerboard family (all bonds of the
  family in parallel), one kernel per diagonal `Op_V` wrap (row and column scaling of G for a
  site), and both flavours per launch.
- Batch flavours, and later several chains per GPU (strided-batched GEMM, one chain per CUDA
  stream / SYCL queue).
- Set honest expectations: BENCHMARKS.md will report the crossover N where one GPU beats the
  16-core node. Below it, CPU MPI stays the better choice, and that is an acceptable outcome.

## 4. Object-oriented Fortran on the hot path

**What:**
- Polymorphic `class(ContainerElementBase)` hopping exponentials in a `DynamicMatrixArray`, with
  type-bound `lmult/rmult` dispatch inside `Hop_mod` loops.
- `class(ham_base)` with Hamiltonian submodules (`ham%Obser`, `ham%S0`) called inside the slice
  loop.
- `Type(Operator)` with **pointer** components that may be MPI shared-memory windows.
- `CLASS(UDV_State)` with a defined assignment.
- Allocations inside inner loops (Upgrade2, `nsigma_single%make`, `ZDSLSYMM` work arrays,
  `QDRP_decompose`, CGR/CGR2_2 buffers). R1 §4, §10.

**Evidence:**
- Device support for polymorphism and deep copy is weak or absent in both compilers: an ifx
  derived-type reallocation bug (Aurora #155), and `declare mapper` only arrived in ifx 2025.2
  (R3 §2).
- On the CPU, allocation is cheap (≤ 1.8 %, PROFILE.md). On a GPU, device allocation inside a
  loop would serialise.

**Strategy:**
- A **flat device mirror**, built once at setup from the host objects:
  - G(N,N,F);
  - per checkerboard family: bond index lists plus 2×2 exponentials for ±1 and ±½;
  - dense exp(±T), exp(±T/2);
  - per Op_V: site index plus the diagonal exponential table per field value;
  - nsigma(N_op, L);
  - UDV workspaces and `udvst`;
  - preallocated scratch.
- All device code works on plain contiguous arrays. Classes and dispatch stay on the host. A
  setup-time check refuses operator types the GPU path does not handle (non-diagonal Op_V,
  continuous fields, time-dependent couplings in the MVP), so unsupported cases fall back to the
  CPU code instead of silently computing wrong results.

## 5. Two hopping representations with opposite GPU profiles

**What:**
- Checkerboard: ~2N 2×2 bond updates per B_T, applied in sequence, hand-coded in `ZDSLSYMM`.
- Dense: one real N×N exp(T) applied to complex G via `ZLARCM`/`ZLACRM`, which internally runs
  two `DGEMM`s on split real and imaginary parts.

**Evidence:**
- Checkerboard B-products cost 2–6 % on the CPU; dense ones cost ≈42 % at every N.
- Dense is 1.5–1.7× slower per sweep at N ≥ 400 (PROFILE.md).
- The two variants give different Trotter errors (E/N differs by 0.006–0.007,
  BASELINE.md), so each needs its own reference.

**Strategy:**
- Dense becomes GEMMs on the device. Right-multiplication by a real matrix is one DGEMM on the
  interleaved real view of the complex matrix (2N×N). Left-multiplication needs either a complex
  copy of exp(T) (ZGEMM, 4× the flops) or a split/transposed layout; Phase 1 benchmarks both.
- Checkerboard needs custom kernels: the bonds within a family are disjoint, so one kernel per
  family updates N/2 bond pairs in parallel. The family order (product order) must match the CPU
  exactly, for the 1e-10 agreement.

## 6. Measurements are user code that reads G 3L times per sweep

**What:**
- `ham%Obser` runs 2L−1 times per sweep and `ham%ObserT` L+1 times (12×12 call counts).
- Each is a Wick-contraction loop with scatter-adds into lattice bins through `latt%imj`
  indirection (`Predefined_Obs_mod`). The loops are written per Hamiltonian, in Fortran.

**Evidence:**
- Measurements take 17 % at N=144 and 8 % at N=1024 (checkerboard).
- Copying G to the host instead costs *[est]*: 3L × F × 16 N² bytes ≈ 4 GB per sweep at N=1024,
  about 0.2 s at ~20 GB/s PCIe.

**Strategy:**
- MVP: copy G (and the τ matrices) to the host per measurement and keep the CPU observable code.
  This is simple and exact, and it gives the Hamiltonian-independent fallback for every model.
- Then port the predefined observables used by Hubbard (Green, SpinZ/XY, Den, kinetic and
  potential energy) to the device, using per-thread partial sums and reductions instead of
  atomics so results are reproducible.
- Custom Hamiltonians keep the host fallback.

## 7. Validation: the RNG and chaotic trajectories

**What:**
- `ranf_wrap` is the intrinsic `RANDOM_NUMBER`, seeded with `RANDOM_SEED`
  (`Libraries/Modules/random_wrap_mod.F90:48-93`), so the stream depends on the compiler.
- Device BLAS sums in a different order, so round-off differs from the CPU.

**Evidence:**
- Same compiler, same seeds: GNU LOG vs GNU default gave identical trajectories (1e-13) over
  600 sweeps × 16 ranks.
- Different streams (GNU vs ifx): scalars agree, but the provisional full-array χ² fails on
  correlation functions because of a slow AF mode with underestimated errors (BASELINE.md).

**Strategy:** three test levels.
1. **Kernel tests with fixed fields** (no Monte Carlo): G after a wrap, after stabilization, after
   a τ step, CPU vs GPU in the same build, relative difference ≤ ~1e-10.
2. **Same-compiler trajectory test:** nvfortran CPU vs nvfortran GPU (and ifx vs ifx), same seeds.
   Observables should coincide until a round-off-driven accept/reject flip; we report the sweep of
   first divergence. This requires the accept decision and RNG calls to stay on the host, in the
   CPU order (item 1).
3. **Statistical tests:** test_vs_ed on the GPU, plus a redesigned compare (curated observables,
   autocorrelation-aware errors, parameters with short autocorrelation). This is a Phase 1 CPU
   ticket; `compare.py` stays provisional until then.

## 8. Host-side numerical checks and phase bookkeeping

**What:** at every stabilization the code runs:
- `Control_PrecisionG` / `Control_Precision_tau` (full-matrix comparison of wrapped vs recomputed
  G, `control_mod`);
- `testscale` NaN and scale checks;
- `Op_phase`, which is O(N_op·L) on the host;
- `DET_C` (ZGETRF) for the phase in CGR (`Prog/cgr1_mod.F90:303`).

**Strategy:** device-side reductions return a few scalars (max and mean deviation, NaN flags,
phase factors). The precision numbers must stay at CPU levels and be reported (correctness rules),
so these checks are kept, not dropped.

## 9. Intel PVC can only be tested at one remove

**What:**
- Agents cannot reach SuperMUC-NG.
- nividia has ifx 2025.3 but no `ocloc`, so no AOT `spir64_gen` build for PVC, only JIT `spir64`
  plus host fallback (MACHINES.md).
- oneMKL offload coverage is incomplete for exactly our LAPACK routines (item 2).
- ifx offload still has open ICE and wrong-answer bugs (R3 §2).
- PVC launch latency is likely higher than on CUDA (R3 §5).
- The ifx CPU profile uses `-static -xHost`, which an offload build must drop.
- The Max 1550 is end-of-life.

**Strategy:**
- Keep the device code to the conservative OpenMP subset both compilers handle.
- Compile every change with ifx (`spir64` JIT) and run it in host fallback on nividia (CI-like).
- Write the SYCL shim as small, separately testable C++.
- Give the owner a validation recipe with expected values and tolerances, and job scripts with
  `OMP_TARGET_OFFLOAD=MANDATORY` so a host fallback fails loudly.
- Pin compiler versions (SuperMUC versions: open question for the owner).

## 10. Keeping `GPU=none` identical and not slower

**What:** every GPU hook touches hot CPU routines (upgrade, wraps, stabilization, measurements)
that also carry many preprocessor variants (STAB1/2/3/LOG/QRREF, MPI, TEMPERING, HDF5; R1 §10).

**Strategy:**
- The GPU path lives in new modules, selected at configure time (`-DALF_GPU`).
- Upstream routines get at most a guarded early branch (`#ifdef ALF_GPU`), so the `GPU=none`
  object code is unchanged.
- The linear-algebra layer's HOST backend calls exactly the same LAPACK/BLAS routines with the
  same arguments.
- Gate per ticket:
  - unit tests;
  - GNU LOG/default trajectory identity vs the stored references;
  - `test_branch_parameters.json` via `compare_dirs.py`;
  - CPU s/sweep within noise of BASELINE.md, for GNU and ifx.

## Lower-impact items, noted

- **Memory:** `udvst` is NSTM·F·(2·16·N²) bytes ≈ 1.3 GB at N=1024, L=200, Nw=10, which fits on
  the A30 (24 GB) and a PVC stack (64 GB). Several chains per GPU will be limited by this.
- **MPI, tempering, global moves, Langevin/HMC, projector:** out of scope for now (FEATURES.md).
  Tempering would need GPU-aware MPI or host staging of G and `udvst`.
- **Mixed precision:** not planned. DQMC stabilization needs FP64, and FP32 for delayed-update
  intermediates is untested in the literature (R2 §2).
