# ALF GPU port – plan

Status 2026-10-29: Phase 2 in progress – 201–206 done, 207 (τ pass) **on hold** (owner decision and verification
standard pending). Process change proposed on branch `gpu-port-process-lanes` (two lanes, [VERIFICATION.md](VERIFICATION.md));
no lane starts before the owner's go-ahead. Phase 1 gate passed 2026-10-27 (ADR-001/002 accepted). Why each piece is here: [HARD.md](HARD.md) (items referred to
as H1–H10). Measurements: [PROFILE.md](PROFILE.md), [BASELINE.md](BASELINE.md).

## 1. Architecture (decided: [ADR-001](ADR-001-programming-model.md), [ADR-002](ADR-002-alf-la-and-data-model.md))

```
host (unchanged ALF control flow, RNG, accept/reject, I/O, MPI, classes)
  │  setup: flatten host objects → device mirror      per site: ratio inputs ↔ decision
  ▼
gpu_state module (plain contiguous device arrays, allocated once)
  G(N,N,F) · B_T tables (checkerboard families / dense exp) · Op_V diag exp tables · nsigma
  UDV workspaces + udvst · τ matrices · preallocated scratch
  │
  ├─ device kernels (OpenMP target, conservative subset): checkerboard family update,
  │   Op_V row/column scaling, rank-1 update + ratio, column norms / pre-pivot, reductions,
  │   predefined Hubbard observables
  └─ alf_la backend layer (~12 ALF-level ops on device arrays + stream/queue handle)
        HOST   : the same LAPACK/BLAS calls as upstream (GPU=none path, test oracle)
        CUDA   : cuBLAS / cuSOLVER (nvfortran modules or ISO_C_BINDING)
        ONEMKL : !$omp dispatch for BLAS; small C++ SYCL USM shim for geqrf/ormqr/ungqr/getrf/getrs (+ trsm)
                 (the module provides la_getrs/la_trsm, not la_getri; 202)
```

- **Directive language:** OpenMP 5.x target offload: `target enter/exit data`, `target update`,
  `target teams distribute parallel do` / `loop`, `use_device_addr`. No `declare mapper`, USM,
  polymorphism, device allocation or type-bound procedures in device code (R3 §7; HARD H4, H9).
  The nividia smoke tests passed for nvfortran (`-mp=gpu` + cuBLAS via `use_device_addr`) and for
  ifx (`spir64` JIT, host fallback) (MACHINES.md).
- **Configure:** `GPU=none|nvidia|intel` orthogonal to MACHINE. `none` adds nothing (H10).
  `nvidia` uses nvfortran with `-mp=gpu -gpu=cc80[,cc90] -cudalib=cublas,cusolver -DALF_GPU
  -DALF_LA_CUDA`. `intel` uses ifx with `-fiopenmp -fopenmp-targets=spir64[_gen] -qmkl -DALF_GPU
  -DALF_LA_ONEMKL` plus the SYCL shim, without `-static`.
  - New compiler profile `NVHPC` (nvfortran), also usable with `GPU=none` as the same-compiler CPU
    reference.
  - New machine profile `SUPERMUC-NG-PVC`.
- **Device residency:** G and the B tables stay on the device for a whole run. Every transfer is
  listed in the table below and in the code, next to the transfer.
- **Stabilization on the device:** pre-pivoting by column norm followed by unpivoted
  `geqrf`/`ormqr`/`ungqr` in decompose, CGR and CGR2_2 (104; approved by the owner). Exception:
  under STABLOG, CGR and CGR2_2 use host ZGEQP3 (108, ADR-001 §4).
- **Updates:** Phase 2 does rank-1 updates on the device with the decision on the host (H1, H7).
  The code is structured as `propose → ratio → decide → apply/defer → flush`, so delayed and
  submatrix updates slot in later.

### Planned host↔device transfers (MVP)

| When | Direction | Size | Why |
|---|---|---|---|
| setup | H→D | B tables, Op_V tables, nsigma, initial UDV | build the mirror |
| per site and slice (2·L·N per sweep) | D→H, H→D | a few complex scalars | ratio to host, decision (and new field) back; keeps the RNG order identical to the CPU |
| per site, per flavour (205, as built) | D→H | 1 complex element G(P(1),P(1),nf), 16 B | the only G entry the host ratio needs (N_non_zero ≤ 1); wrap: 1 launch per site and flavour for N = 1 operators, 2 otherwise |
| per accepted flip (205, as built) | H→D | 1 complex scalar nsigma_d(n,t), 16 B | wrap kernels read the field on the device; per flavour 2 more launches (rank-1 update) |
| per WRAPGRUP / WRAPGRDO call (205, **interim**, removed by 206/208) | H→D at entry, D→H at exit | F·16N² B each way | stabilization, τ pass and measurements still use the host GR; marked `INTERIM 205` in `gpu_update_mod` |
| per measurement (3L per sweep), MVP only | D→H | G (and GT0/G0T/GTT for τ), F·16N² B each | host observables; removed once Hubbard observables run on the device |
| per stabilization | D→H | a few scalars | precision, NaN checks, phase (H8) |
| per CGR / CGR2_2 call, **STABLOG only** | D→H, H→D | N² (CGR) / 4N² (CGR2_2) complex, each way | host ZGEQP3: no GPU pivoting meets the precision criteria (108, ADR-001 §4) |
| 206: per QR of the default decompose and CGR (`la_qrpp`) | D→H, H→D | 8N B (column norms), 4N B (permutation) | the pre-pivot sort runs on the host |
| 206: per STABLOG decompose | D→H, H→D | 8N B (log norms), 4N B (permutation) | the CPU's selection sort on the host |
| 206: per default decompose / per CGR | D→H | 24 B (scale test) / 16 B (phase) | NaN and range test of the scales; Phase of the determinant |
| 206 **INTERIM** (removed by 207/208): per CGR | D→H, H→D | 16N² B each way (new G, wrapped G) + 32 B | GR is host-authoritative until 208; the precision check reduces on the device (`Control_PrecisionG_vals`) |
| 206 **INTERIM** (host-update path only, `ALF_GPU_NOUPDATE=1`; skipped with device updates, 205): per WRAPUR/WRAPUL block | H→D | N_op·Nwrap·16 B | HS fields changed by the host updates (`nsigma_d`) |
| 206 **INTERIM** (removed by 207): per backward stabilization, Ltau = 1 | D→H | 32N² + 16N B per flavour | `TAU_M` (host) reads `udvst` |
| 302b (control runs only, `ALF_GPU_FAULT_STAB=eps:<x>` in the sweep): per CGR with precision check | H→D | 16N² B | perturbed `st_Gn` back to the device for the precision reduction |
| end of bin / run | D→H | nsigma (confout), and G only if needed | output |

## 2. Phases, gates and tickets

Ticket IDs: Phase p uses `p01…`. Each ticket gets a file `docs/gpu/tickets/NNN-name.md` (template
`000-template.md`) when it is opened. Owner agent in brackets. "dev" means `fortran-gpu-dev`,
"test" means `test-runner`, "rev" means `numerics-reviewer`.

### Phase 0 – baseline and understanding ✅
Fork, agents, CLAUDE.md, nividia probe, CPU baseline (GNU + ifx unit tests, test_vs_ed both
PASS, Hubbard references), CPU profile, HARD.md, this plan.

### Phase 1 – architecture decision by prototype ✅ (gate passed 2026-10-27)

| ID | Ticket | Owner | Depends |
|---|---|---|---|
| 101 | **Kernel benchmark harness** `gpu_bench/` outside `Prog/`, with HOST and CUDA backends and device-resident complex N×N (N = 64…1024, F = 2): (a) B_T propagation, checkerboard family kernel and dense GEMM (real exp(T): DGEMM on the interleaved view vs ZGEMM); (b) Op_V diagonal wrap; (c) stabilization: column-norm pre-pivot + `geqrf` + `ungqr`/`ormqr` vs host ZGEQP3; (d) rank-1 vs rank-k update (k = 1…128) incl. per-update host round-trip latency; (e) flavour batching. Correctness vs HOST at 1e-10 (or the documented bound for QR) | dev | – |
| 102 | OpenACC variant of 101 (a, d) on nvfortran, only if OpenMP is ≥1.5× slower than OpenACC on a kernel that matters. Done as a narrow launch study; OpenACC not adopted | dev | 101 |
| 103 | **ifx build of 101** (`spir64` JIT), run in host fallback on nividia; oneMKL dispatch + SYCL shim prototype for `geqrf`/`ungqr`/`ormqr` (compile and host-fallback correctness only, no PVC here) | dev | 101 |
| 104 | **CPU accuracy study: pre-pivoted unpivoted QR in CGR and CGR2_2** (CPU only, behind a test-only macro). Compare `Precision Green`/`tau` and observables vs upstream QRCP for β = 4…20+, U = 4…8, default and LOG, cb and dense. Decides H2 | dev + rev | – |
| 105 | **NVHPC profile** in `configure.sh` (nvfortran, `GPU=none`): unit tests, test_vs_ed, 12×12 references. Gives the same-compiler CPU baseline for H7 | dev, test | – |
| 106 | **Statistical comparison redesign** (`testsuite/gpu_reference/`): curated observables, autocorrelation-aware errors (binning analysis / τ_int), plus a short-autocorrelation reference set (e.g. U=2 or β=2, many bins). Rerun GNU vs ifx, which must pass | dev | – |
| 107 | **ADR-001** programming model + **ADR-002** `alf_la` interface and data model, written from 101–104 (+ 102, 108). Accepted | manager | 101–104 |
| 108 | STABLOG pivot study (variants of pre-pivoting in CGR/CGR2_2 under STABLOG). Result: host ZGEQP3 under STABLOG | dev + rev | 104 |

Gate deliverables: benchmark numbers (BENCHMARKS.md, per kernel, crossover estimates), ADRs, the
104 verdict, the revised Phase 2 ticket list.

### Phase 2 – NVIDIA MVP (gate)
Scope:
- Hubbard on the square lattice (predefined Hamiltonian), finite T;
- local rank-1 updates;
- checkerboard and dense hopping;
- Ltau = 0/1;
- default and STABLOG stabilization;
- single GPU, no MPI.

Unsupported inputs are refused at setup (H4).

| ID | Ticket | Depends |
|---|---|---|
| 201 | `configure.sh` `GPU=` switch, `-DALF_GPU`, NVHPC GPU flags (nvfortran via the existing PGI profile, 105); `GPU=none` byte-identical behaviour gate | 105, ADR |
| 202 | `alf_la` module: HOST backend (exact upstream calls) + CUDA backend, incl. `la_qrpp` and `la_geqp3_host`; unit tests per op | ADR |
| 203 | `gpu_state` device mirror + setup flattening + unsupported-feature checks | 201 |
| 204 | Device B_T (checkerboard families, dense) and Op_V wraps; kernel tests with fixed fields vs CPU ≤ 1e-10 | 202, 203 |
| 205 | WRAPGR up/down on the device: ratio, host decision, rank-1 apply; `flush` interface | 204 |
| 206 | Stabilization on the device: WRAPUR/UL + decompose, CGR (both NVAR), pre-pivoted QR per 104; under STABLOG CGR via `la_geqp3_host` (108); precision checks as reductions | 202, 204 |
| 207 | τ pass on the device: PROPR/PROPRM1, CGR2_2 (STABLOG: `la_geqp3_host`), `Control_Precision_tau` | 206 |
| 208 | Measurements: MVP host copy of G; then predefined Hubbard observables (eq + τ) on the device | 205, 207 |
| 209 | Validation suite: kernel tests, same-compiler trajectory test (NVHPC CPU vs GPU, first-divergence sweep), test_vs_ed on GPU, statistical compare vs references (106), precision numbers | 205–208 |
| 210 | BENCHMARKS.md: s/sweep CPU (16 cores, nividia) vs GPU for N = 64…1024, ±Ltau, per-kernel breakdown (incl. rank-1 share), crossover. Development numbers; the official CPU baseline is decided in Phase 3 | 209 |
| 211 | CPU regression gate: unit tests, `test_branch_parameters.json` / `compare_dirs.py`, GNU + ifx s/sweep within noise | 201–208 |
| 213 | STABLOG same-seed excess above 1e-9 at 12×12 (206 code; first bin 184, max 1.214e-9): cause, fix or a pre-registered calibrated criterion ([draft](tickets/213-log-same-seed-excess.md)) | 212 |
| 212 | **Verification hardening** (code): gate all `ALF_GPU_*` debug/fault variables so production runs refuse them; provenance header in the binary's output; "GPU debug env" and "fault armed" info lines; graded fault `ALF_GPU_FAULT_TAU=eps:<x>`; `ALF_GPU_TAUDIAG_SEED`; the cause of the default-path non-reproducibility and of 207's default paired spread, then a deterministic default path (fixed-order `k_colnorms`/`k_lognorm`/`k_prec`, explicit cuSOLVER/cuBLAS modes; decision 2a of 2026-10-29); run-twice-identical and "τ pass leaves the sweep unchanged" gates. Must land before any further ticket closes | harness approved |

Every ticket is accepted under [VERIFICATION.md](VERIFICATION.md): pre-registered acceptance file, T2 harness run by
the manager on a frozen candidate, positive controls, independent verifier and reviewer. 207 closes as described in
[tickets/207-closing-plan.md](tickets/207-closing-plan.md); order in the NVIDIA lane: 212 → 207 → 208 → 209 …

### Two lanes under one manager (proposed 2026-10-29; starts only with the owner's go-ahead)

The NVIDIA lane (Phase 2 tickets, `nividia`, `~/alf-gpu/work/t2NN`) and the Intel lane (pulled-forward Phase 3
tickets, `intel`, `~/alf-gpu/work/t3NN`) run in parallel on two machines that never compete for a GPU. Only the
manager commits; each lane works on its own sub-branch (`gpu-port-NNN-name`) and both are accepted under
[VERIFICATION.md](VERIFICATION.md) with the same gates.

**Why now:** under `GPU=intel` only the 205 updates and the checkerboard wraps run on the GPU today, and the oneMKL
backend is a guard without calls. The longer Intel lags, the more the design fits only NVIDIA (two nvfortran bugs and
two cuSOLVER quirks are already worked around); oneMKL on a real Intel GPU will have its own, and finding them early is
cheaper.

**Prerequisites of the Intel lane:** a consistent oneAPI on `intel` (matching ifx and MKL releases: since 2026-10-29
ifx 2026.2 + MKL 2026.2, see MACHINES.md), the owner's approval of the process diff, and the harness. All met
2026-10-30.

**File ownership while both lanes run**

| File / area | Owner | Rule |
|---|---|---|
| `Prog/gpu_tau_mod.F90`, `Prog/tau_m_mod.F90`, observables (208) | NVIDIA lane | – |
| `Libraries/Modules/alf_la_mod.F90`: the `ALF_LA_ONEMKL` bodies, `la_init` interop; `alf_la_sycl_shim.cpp`; `configure.sh` intel block; `testsuite/alf_la.tests` ONEMKL mode | Intel lane | HOST and CUDA code paths unchanged; `GPU=none` identity per change |
| the `ALF_LA_ONEMKL` guards in `gpu_stab_mod.F90:137–140`, `gpu_state_mod.F90` (la_init/finalize), `gpu_wrap_mod.F90:98–102` | Intel lane | one-line guard edits only |
| everything else in `gpu_stab_mod`, `gpu_state_mod`, `gpu_wrap_mod`; `PLAN.md`; `Prog/Makefile` | shared | the manager sequences edits: one lane at a time, the other rebases |
| `alf_la` public interface (`alf_la_mod.F90` public list and argument lists) | frozen | no signature change while both lanes run; a new routine needs the manager's decision |

**τ path on Intel:** removing the 302 guard would switch on the 206 stabilization **and** the 207 τ pass at once. The
τ path under `GPU=intel` gets its own switch (refused at setup until 302c), so 302b validates 204–206 alone.

### Phase 3 – Intel / SuperMUC-NG Phase 2 enablement (Intel lane; owner runs the PVC jobs)

| ID | Ticket | Depends |
|---|---|---|
| 301 | `SUPERMUC-NG-PVC` profile: ifx, `-fiopenmp -fopenmp-targets=spir64_gen -Xs "-device pvc"`, oneMKL, `mpiifx`, HDF5. Flags checked against current LRZ/Intel docs and the owner's `module avail` | owner's `module avail` |
| 302a | ONEMKL backend of `alf_la` + `alf_la.tests` routine tests and `t_order` on the B580; per routine measured GPU / host / hybrid; STABLOG keeps host ZGEQP3; interop queue question (ADR-002 §1) ([draft](tickets/302a-onemkl-backend.md)); acceptance in two parts: B580 tests on `intel`, GPU=none and CUDA identity on `nividia` | 202 (interface), consistent oneAPI (owner, decision 6), harness |
| 302b | 204–206 device paths on the B580 with the 206 gates (self-tests, same-seed vs CPU, precision); τ path still refused | 302a |
| 302c | τ path (207) on the B580 | 302b, 207 committed |
| 303 | ifx build of the full GPU path on `intel` (B580, `spir64` JIT and AOT), 209 suite on the B580 (correctness; not PVC performance) | 302c, 209 |
| 304 | Slurm scripts (`test` partition, one rank per tile by default, CCS mode as an option, `I_MPI_OFFLOAD_CELL=tile`, `OMP_TARGET_OFFLOAD=MANDATORY`) + validation recipe with expected values and tolerances | 301, 303 |
| 305 | Evaluate the owner's results; BENCHMARKS.md PVC column; fixes | 304 |
| 306 | Small dense linear algebra on PVC, host vs device per routine and size (planned 2026-10-30, owner): on the B580, oneMKL `zgeqrf` costs a fixed ~0.45 ms per call at any size, and host MKL with copies was faster at every size measured (8 to 400: 0.02–2.3 ms host vs 0.6–7.8 ms oneMKL; 302a, implementer 3). Measure the same on PVC (QR family, `la_qrpp`, LU, trsm; N = 36…1024; including the H↔D copies) and decide per routine whether a size threshold "small → host, large → device" pays off on the GPU node; the threshold is a measured, documented constant, not a guess | 304, 302a |

Validation tail of the Intel lane: the same gates as NVIDIA under the same harness (kernel self-tests with fault
controls, same-seed against the CPU, precision, compare2, `GPU=none` identity), run on `intel`; B580 numbers are
correctness evidence, not PVC performance.

### Later phases – outline only, wait for the owner's go-ahead

| ID | Work | Notes |
|---|---|---|
| 401 | **Delayed updates on the CPU** (Sun & Xu, PRB 109, 235140): complex G, op_dim > 1 vertices (Δ diagonal in the rotated basis, already ALF's structure), N_SUN via Ratio^N_SUN | Flush before anything that reads G: wraps (Op_Wrapdo/up), measurements, stabilization, global moves. `GPU=none` default stays upstream until the owner decides |
| 402 | **Submatrix updates on the CPU** (Sun & Xu, SciPost Phys. 18, 055) | Γ on the host, small k×k |
| 403 | GPU versions of 401/402 | Host Γ + device GEMM flush; the expected main GPU speedup (H1) |
| 404 | Performance work: flavour batching, several chains per GPU, fused checkerboard kernels, device measurements for more observables, CUDA graphs / `nowait` | |
| 405 | MPI + multi-GPU: one rank per GPU/stack, GPU-aware MPI for tempering | Testable only on SuperMUC |
| 406 | Idea: the Metropolis accept decision on the device, with each slice's random numbers uploaded in advance, so the RNG stream and same-seed comparisons stay unchanged | idea only |
| 407 | **Vertices with `Op_V%N_non_zero > 1` on the device** (other Hamiltonians; planned 2026-10-30, owner): today the device update refuses them (`gpu_state_mod`, "device update supports <= 1"). Rank-n update: n×n ratio determinant and Woodbury apply per vertex. Measure where the small n×n work should run (custom device kernel, batched library calls, host), on A30, B580 and PVC; the fixed per-call library cost seen in 302a/306 makes the answer rank-dependent. Together with 401–403 (delayed/submatrix updates batch exactly these updates) | owner's choice of Hamiltonians; 306 for PVC |

## 3. Validation strategy (summary, details in HARD H7 and BASELINE)

1. Kernel level: fixed fields, CPU vs GPU in one build, ≤ ~1e-10 relative (documented per test).
2. Trajectory level: same compiler, same seeds, CPU vs GPU, report the first divergent sweep.
3. Statistical level:
   - test_vs_ed on the GPU (|Δ| < 1e-3, < 3σ);
   - redesigned compare vs the CPU references, same hopping variant;
   - χ² reported;
   - `Precision Green`/`tau`/phase reported next to the CPU values.
4. CPU path: unit tests, branch comparison, GNU LOG/default identity, performance within noise.
5. Every gate has a positive control that must FAIL; verdicts are computed by `verdict.py` from stamped raw outputs
   ([VERIFICATION.md](VERIFICATION.md)).

## 4. Risks

| # | Risk | Likelihood / impact | Mitigation |
|---|---|---|---|
| R1 | Pre-pivoted QR not accurate enough in CGR/CGR2_2 at large β | default: resolved (104); STABLOG: occurred (108) | STABLOG uses host ZGEQP3 for CGR/CGR2_2; MAGMA (NVIDIA) as later option |
| R2 | GPU slower than the 16-core node below large N (rank-1 latency, small GEMMs) | high / medium | Report the crossover honestly; delayed/submatrix updates (403) and chains per GPU (404) are the levers |
| R3 | nvfortran OpenMP offload markedly slower than OpenACC | medium / medium | 102 measures; keep directives few and localized |
| R4 | ifx offload bugs (ICE, wrong answers), no local PVC, no `ocloc` | high / medium | Conservative subset, ifx JIT builds + host fallback on every change, version pinning, owner recipe |
| R5 | Same-compiler trajectory test diverges early because device round-off flips a decision | medium / low | Report the first-divergence sweep; the statistical tests remain authoritative |
| R6 | oneMKL/SYCL and OpenMP interop on one device immature | medium / medium | Prototype in 103, shim kept small and separately testable |
| R7 | Touching hot CPU routines slows `GPU=none` | low / high | New modules + guarded early branches; 211 gate per ticket |
| R8 | Disk (19 GB free, 93 % used, 2026-10-29) and node size on nividia | medium / medium | Evidence kept until the gate (≤ 1.5 GB live per ticket), archive ≤ 2 GB, CPU reference cache ≤ 2 GB, driver refuses below 5 GB free; the owner decides on the old t104–t108 directories (~2.7 GB) |
| R9 | Max 1550 is EOL; SuperMUC versions differ from docs | – / medium | `alf_la` limits lock-in; owner checks `module avail` |

## 5. Owner decisions

Decided at the Phase 1 gate (2026-10-27):
1. `ocloc` on nividia: not installed. AOT `spir64_gen` compile checks run on the `intel` test machine, which has it.
2. Pre-pivoted QR in decompose/CGR/CGR2_2 for GPU builds (default stabilization): **approved**. `GPU=none` keeps
   upstream QRCP.
3. STABLOG: host ZGEQP3 for CGR/CGR2_2; no further pivot study unless profiling shows the host QR dominates.
4. Per-site update: host decision in the MVP; OpenACC not adopted; device decision and delayed/submatrix updates
   later (401–404, need the owner's go-ahead).
5. ADR-001 and ADR-002 accepted as written.

Still open:
1. oneAPI/ifx, MKL, MPI and HDF5 module versions on SuperMUC-NG Phase 2: needed in Phase 3 (301); until then pin
   **ifx 2026.2 + oneMKL 2026.2** for `GPU=intel` on `intel` (pin change of 2026-10-30, DECISIONS.md row 13; was ifx
   2025.3). The ifx on `nividia` (2025.3, CPU-only stand-in builds) is not the Intel-lane toolchain.
2. Official CPU baseline for BENCHMARKS.md (nividia 16-core VM vs a SuperMUC full-node CPU run): decided in
   Phase 3. Phase 2 reports nividia numbers as development numbers.
