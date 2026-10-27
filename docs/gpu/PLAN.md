# ALF GPU port – plan

Status 2026-10-27, end of Phase 0. Why each piece is here: [HARD.md](HARD.md) (items referred to
as H1–H10). Measurements: [PROFILE.md](PROFILE.md), [BASELINE.md](BASELINE.md).

## 1. Architecture (hypothesis; Phase 1 ADR confirms or overturns it)

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
        ONEMKL : !$omp dispatch for BLAS; small C++ SYCL USM shim for geqrf/ormqr/ungqr/getrf/getri
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
  `geqrf`/`ormqr`/`ungqr`, everywhere QRCP is used today (decompose, CGR, CGR2_2). This is subject
  to the Phase 1 CPU accuracy study (H2).
- **Updates:** Phase 2 does rank-1 updates on the device with the decision on the host (H1, H7).
  The code is structured as `propose → ratio → decide → apply/defer → flush`, so delayed and
  submatrix updates slot in later.

### Planned host↔device transfers (MVP)

| When | Direction | Size | Why |
|---|---|---|---|
| setup | H→D | B tables, Op_V tables, nsigma, initial UDV | build the mirror |
| per site and slice (2·L·N per sweep) | D→H, H→D | a few complex scalars | ratio to host, decision (and new field) back; keeps the RNG order identical to the CPU |
| per measurement (3L per sweep), MVP only | D→H | G (and GT0/G0T/GTT for τ), F·16N² B each | host observables; removed once Hubbard observables run on the device |
| per stabilization | D→H | a few scalars | precision, NaN checks, phase (H8) |
| per CGR / CGR2_2 call, **STABLOG only** | D→H, H→D | N² (CGR) / 4N² (CGR2_2) complex, each way | host ZGEQP3: no GPU pivoting meets the precision criteria (108, ADR-001 §4) |
| end of bin / run | D→H | nsigma (confout), and G only if needed | output |

## 2. Phases, gates and tickets

Ticket IDs: Phase p uses `p01…`. Each ticket gets a file `docs/gpu/tickets/NNN-name.md` (template
`000-template.md`) when it is opened. Owner agent in brackets. "dev" means `fortran-gpu-dev`,
"test" means `test-runner`, "rev" means `numerics-reviewer`.

### Phase 0 – baseline and understanding ✅ (this gate)
Fork, agents, CLAUDE.md, nividia probe, CPU baseline (GNU + ifx unit tests, test_vs_ed both
PASS, Hubbard references), CPU profile, HARD.md, this plan.

### Phase 1 – architecture decision by prototype (gate)

| ID | Ticket | Owner | Depends |
|---|---|---|---|
| 101 | **Kernel benchmark harness** `gpu_bench/` outside `Prog/`, with HOST and CUDA backends and device-resident complex N×N (N = 64…1024, F = 2): (a) B_T propagation, checkerboard family kernel and dense GEMM (real exp(T): DGEMM on the interleaved view vs ZGEMM); (b) Op_V diagonal wrap; (c) stabilization: column-norm pre-pivot + `geqrf` + `ungqr`/`ormqr` vs host ZGEQP3; (d) rank-1 vs rank-k update (k = 1…128) incl. per-update host round-trip latency; (e) flavour batching. Correctness vs HOST at 1e-10 (or the documented bound for QR) | dev | – |
| 102 | OpenACC variant of 101 (a, d) on nvfortran, only if OpenMP is ≥1.5× slower than OpenACC on a kernel that matters | dev | 101 |
| 103 | **ifx build of 101** (`spir64` JIT), run in host fallback on nividia; oneMKL dispatch + SYCL shim prototype for `geqrf`/`ungqr`/`ormqr` (compile and host-fallback correctness only, no PVC here) | dev | 101 |
| 104 | **CPU accuracy study: pre-pivoted unpivoted QR in CGR and CGR2_2** (CPU only, behind a test-only macro). Compare `Precision Green`/`tau` and observables vs upstream QRCP for β = 4…20+, U = 4…8, default and LOG, cb and dense. Decides H2 | dev + rev | – |
| 105 | **NVHPC profile** in `configure.sh` (nvfortran, `GPU=none`): unit tests, test_vs_ed, 12×12 references. Gives the same-compiler CPU baseline for H7 | dev, test | – |
| 106 | **Statistical comparison redesign** (`testsuite/gpu_reference/`): curated observables, autocorrelation-aware errors (binning analysis / τ_int), plus a short-autocorrelation reference set (e.g. U=2 or β=2, many bins). Rerun GNU vs ifx, which must pass | dev | – |
| 107 | **ADR-001** programming model + **ADR-002** `alf_la` interface and data model, written from 101–104 | manager | 101–104 |

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
| 201 | `configure.sh` `GPU=` switch, `-DALF_GPU`, NVHPC GPU flags; `GPU=none` byte-identical behaviour gate | 105, ADR |
| 202 | `alf_la` module: HOST backend (exact upstream calls) + CUDA backend; unit tests per op | ADR |
| 203 | `gpu_state` device mirror + setup flattening + unsupported-feature checks | 201 |
| 204 | Device B_T (checkerboard families, dense) and Op_V wraps; kernel tests with fixed fields vs CPU ≤ 1e-10 | 202, 203 |
| 205 | WRAPGR up/down on the device: ratio, host decision, rank-1 apply; `flush` interface | 204 |
| 206 | Stabilization on the device: WRAPUR/UL + decompose, CGR (both NVAR), pre-pivoted QR per 104; precision checks as reductions | 202, 204, 104 |
| 207 | τ pass on the device: PROPR/PROPRM1, CGR2_2, `Control_Precision_tau` | 206 |
| 208 | Measurements: MVP host copy of G; then predefined Hubbard observables (eq + τ) on the device | 205, 207 |
| 209 | Validation suite: kernel tests, same-compiler trajectory test (NVHPC CPU vs GPU, first-divergence sweep), test_vs_ed on GPU, statistical compare vs references (106), precision numbers | 205–208 |
| 210 | BENCHMARKS.md: s/sweep CPU (16 cores, nividia) vs GPU for N = 64…1024, ±Ltau, per-kernel breakdown (incl. rank-1 share), crossover | 209 |
| 211 | CPU regression gate: unit tests, `test_branch_parameters.json` / `compare_dirs.py`, GNU + ifx s/sweep within noise | 201–208 |

Every ticket gets review (rev) and acceptance re-run by the manager.

### Phase 3 – Intel / SuperMUC-NG Phase 2 enablement (gate; owner runs the jobs)

| ID | Ticket |
|---|---|
| 301 | `SUPERMUC-NG-PVC` profile: ifx, `-fiopenmp -fopenmp-targets=spir64_gen -Xs "-device pvc"`, oneMKL, `mpiifx`, HDF5. Flags checked against current LRZ/Intel docs and the owner's `module avail` |
| 302 | ONEMKL backend of `alf_la`: `!$omp dispatch` BLAS + SYCL USM shim for the LAPACK pieces that dispatch would run on the host; documented per routine |
| 303 | ifx build of the full GPU path on nividia (`spir64` JIT), host-fallback run of the 209 suite |
| 304 | Slurm scripts (`test` partition, 8 ranks/node = 1 per stack, FLAT mode, `I_MPI_OFFLOAD_CELL=tile`, `OMP_TARGET_OFFLOAD=MANDATORY`) + validation recipe with expected values and tolerances |
| 305 | Evaluate the owner's results; BENCHMARKS.md PVC column; fixes |

### Later phases – outline only, wait for the owner's go-ahead

| ID | Work | Notes |
|---|---|---|
| 401 | **Delayed updates on the CPU** (Sun & Xu, PRB 109, 235140): complex G, op_dim > 1 vertices (Δ diagonal in the rotated basis, already ALF's structure), N_SUN via Ratio^N_SUN | Flush before anything that reads G: wraps (Op_Wrapdo/up), measurements, stabilization, global moves. `GPU=none` default stays upstream until the owner decides |
| 402 | **Submatrix updates on the CPU** (Sun & Xu, SciPost Phys. 18, 055) | Γ on the host, small k×k |
| 403 | GPU versions of 401/402 | Host Γ + device GEMM flush; the expected main GPU speedup (H1) |
| 404 | Performance work: flavour batching, several chains per GPU, fused checkerboard kernels, device measurements for more observables, CUDA graphs / `nowait` | |
| 405 | MPI + multi-GPU: one rank per GPU/stack, GPU-aware MPI for tempering | Testable only on SuperMUC |

## 3. Validation strategy (summary, details in HARD H7 and BASELINE)

1. Kernel level: fixed fields, CPU vs GPU in one build, ≤ ~1e-10 relative (documented per test).
2. Trajectory level: same compiler, same seeds, CPU vs GPU, report the first divergent sweep.
3. Statistical level:
   - test_vs_ed on the GPU (|Δ| < 1e-3, < 3σ);
   - redesigned compare vs the CPU references, same hopping variant;
   - χ² reported;
   - `Precision Green`/`tau`/phase reported next to the CPU values.
4. CPU path: unit tests, branch comparison, GNU LOG/default identity, performance within noise.

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
| R8 | Disk (25 GB free) and node size on nividia | medium / low | Clean up after each run; keep references compact |
| R9 | Max 1550 is EOL; SuperMUC versions differ from docs | – / medium | `alf_la` limits lock-in; owner checks `module avail` |

## 5. Open decisions for the owner
1. Install `ocloc` system-wide on nividia (enables AOT `spir64_gen` PVC builds)?
2. oneAPI/ifx, MKL, MPI and HDF5 module versions on SuperMUC-NG Phase 2 (`module avail intel-toolkit`)?
3. Is a numerically different (pre-pivoted) QR in CGR/CGR2_2 acceptable for GPU builds if 104 shows
   equal precision? (`GPU=none` keeps upstream QRCP either way.)
4. Is nividia's 16-core VM the official CPU comparison for BENCHMARKS.md, or should the owner run a
   full-node CPU number (e.g. SuperMUC Phase 2's 112 cores) later?
