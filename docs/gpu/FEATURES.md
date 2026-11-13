# Features: CPU vs NVIDIA vs Intel vs AMD

Date 2026-11-13. Generated from the code at commit b7bf0bb3 (branch `gpu-port-unified`). The earlier text of this file (2026-11-06) is replaced.
Line numbers are anchors and drift; search for the name.

## How to read this page

Entries in the four backend columns:
- **device**: the work runs on the GPU.
- **host-assisted (what stays on the host)**: the GPU path is active, but the named part runs on the host.
- **refused -> runs on CPU path**: `gpu_supported` returns a reason, `info` says `GPU path (gpu_state): INACTIVE, reason: ...`, and the whole run uses
  the CPU code. Never a silent wrong answer. Reason strings and conditions: [NOT_IMPLEMENTED.md](NOT_IMPLEMENTED.md).
- **stops the run**: the Hamiltonian itself refuses the setting (all builds, CPU included). There is no CPU fallback.
- **not built**: no code for this backend.
- **not verified**: the code or the documents do not settle it. Nothing was guessed.

Three facts decide most entries.
1. The three GPU backends share the same Fortran modules (`gpu_state/wrap/lowrank/update/stab/tau/obs_mod`). The checks of `gpu_supported`
   (`Prog/gpu_state_mod.F90:499-617`) do not depend on the backend. A feature accepted there is compiled for NVIDIA, Intel and AMD alike.
   The backends differ in the `alf_la` layer and in a few kernels (section "Backend differences").
2. "device" therefore says what the code does. Whether a backend has **run** the feature is a separate question, answered in "Run evidence".
3. `gpu_supported` tests the structure of the model (`Op_V`, `Op_T`), not its name. Any user Hamiltonian that passes these tests takes the GPU path.
   Device observables exist only for the predefined Hubbard model (`gpu_obs_setup` is called only from `Hamiltonian_Hubbard_smod.F90:788`).

Backends: NVIDIA = nvfortran, OpenMP offload + cuBLAS/cuSOLVER (`ALF_LA_CUDA`); Intel = ifx, OpenMP offload + oneMKL + SYCL shim (`ALF_LA_ONEMKL`);
AMD = gfortran amdgcn offload + hipBLAS/hipSOLVER + HIP kernels (`ALF_LA_HIP`). CPU = `GPU=none`.

## Table

| Feature | CPU | NVIDIA | Intel | AMD | Evidence |
|---|---|---|---|---|---|
| Finite-T DQMC, predefined Hubbard (`N_FL` = 1 or 2, field types 1, 2, rank-1 diagonal vertices) | yes | device | device | device | `gpu_state_mod.F90:552-580`, `main.F90:709-714` |
| Projector (T = 0) | yes | refused -> runs on CPU path | refused -> runs on CPU path | refused -> runs on CPU path | `gpu_state_mod.F90:526-528` |
| Hamiltonians other than Hubbard | yes | accepted only if every `Op_V` is diagonal, field type 1 or 2, `N_non_zero` <= 1, no `g_t`, no weight reconstruction; then host-assisted (observables on the host, one D->H of G per measured slice) | same | same | `gpu_state_mod.F90:545-583`, `main.F90:1033,1109`. Stock Kondo, tV, Spin_Peierls, Z2_Matter, LR_Jz_Test, LRC, Plain_Vanilla: not run in this audit, refusal reason not verified per model |
| HS field types 1, 2 | yes | device | device | device | `gpu_state_mod.F90:559-564` |
| HS field types 3, 4 | yes | refused -> runs on CPU path | refused -> runs on CPU path | refused -> runs on CPU path | same lines |
| Non-diagonal vertex, or `N_non_zero` > 1 (rank > 1) | yes | refused -> runs on CPU path | refused -> runs on CPU path | refused -> runs on CPU path | `gpu_state_mod.F90:555-558,570-574` |
| Local updates, standard (rank-1) | yes | device (decision on host, same random-number order) | device | device | `gpu_update_mod.F90:43-50`, `gpu_lowrank_mod.F90:225-229` |
| Delayed and submatrix updates (`Sequential_scheme`, `N_delay`) | yes (CPU: rank up to 8) | device (submatrix: block matrix `M` is kept on the host) | device; async block path opt-in (`ALFSH_ASYNCBLK`) | device (fused `hk_blk_end` kernel) | `gpu_lowrank_mod.F90:33-55,128,1334`; `CONSOLIDATION.md` merge 6 |
| Cross-slice delay (`Lowrank_cross_slice`) | yes | refused -> runs on CPU path | refused -> runs on CPU path | refused -> runs on CPU path | `gpu_state_mod.F90:541-543` |
| Global moves, global-in-tau moves (`N_Global_tau` > 0) | yes | refused -> runs on CPU path | refused -> runs on CPU path | refused -> runs on CPU path | `gpu_state_mod.F90:532-534` |
| Langevin, HMC | yes | refused -> runs on CPU path | refused -> runs on CPU path | refused -> runs on CPU path | `gpu_state_mod.F90:529-531` |
| MALA (sequential, global, global-tau) | yes | refused -> runs on CPU path | refused -> runs on CPU path | refused -> runs on CPU path | `gpu_state_mod.F90:538-540` |
| `Propose_S0` | yes | refused -> runs on CPU path | refused -> runs on CPU path | refused -> runs on CPU path | `gpu_state_mod.F90:535-537` |
| Parallel tempering with state exchange (`TEMPERING` without `PARALLEL_PARAMS`) | yes (MPI) | refused -> runs on CPU path (compile-time test, whole build) | same | same | `gpu_state_mod.F90:516-522` |
| Parallel parameters (`PARALLEL_PARAMS`, "PP": independent groups, no exchange) | yes (MPI) | device (own mirror per rank, `Exchange_Step` skipped) | device (same code; PP build documented) | device (same code) | `gpu_state_mod.F90:516-520`, `PP_GPU.md`, `CONSOLIDATION.md` merge 1 |
| MPI independent walkers (one chain per rank) | yes | device (rank `r` uses OpenMP device `r mod n_devices`; MPS for several ranks per GPU) | device (rank `r` on tile `r mod 8`; no MPS, one rank per tile is best) | not verified | `gpu_state_mod.F90:552`, `MPI_GPU.md`, `platforms/intel_pvc_supermuc_ng2.md`; AMD documents give no MPI run |
| Checkerboard hopping | yes | device (group kernels; dense matrix for `Ndim` <= 800) | device | device | `gpu_state_mod.F90:86,958-978` (dense part: `HOPD_NMAX = 800`, `gpu_wrap_mod.F90:104`) |
| Dense hopping | yes | device (GEMM) | device (GEMM) | device (GEMM) | `gpu_state_mod.F90:978` |
| Symmetric Trotter (`Symm`) | yes | device | device | device | `gpu_obs_mod.F90:268-290`, `gpu_wrap_mod.F90` (`gpu_symm`) |
| Time-dependent hopping `Op_T` with `g_t` (space-time point defect) | yes | device | device | device | `gpu_state_mod.F90:586-611`, `Hop_mod.F90:364-391`, `TIMEDEP_HOPPING.md` |
| Time-dependent coupling of `Op_V` (`g_t`) | yes | refused -> runs on CPU path | refused -> runs on CPU path | refused -> runs on CPU path | `gpu_state_mod.F90:576-580` |
| Flavour symmetry reconstruction (`Calc_Fl`, `weight_reconstruction`) | yes | refused -> runs on CPU path | refused -> runs on CPU path | refused -> runs on CPU path | `gpu_state_mod.F90:544-546` |
| Stabilization, default (pre-pivoted QR) | yes | device | device (oneMKL QR; opt-in own device QR `ALFSH_DEVQR`) | device (own HIP QR for m <= 2304, else hipSOLVER) | `gpu_stab_mod.F90:148`, `alf_la_hip_mod.F90:518-523` |
| Stabilization LOG | yes | host-assisted (pivoted QR `la_geqp3_host` on the host) | host-assisted (same) | host-assisted (same) | `gpu_stab_mod.F90:43-46,857-858`, `alf_la_mod.F90:2078-2095` |
| STAB1, STAB2, STAB3 | yes | refused -> runs on CPU path (compile-time, whole build) | same | same | `gpu_state_mod.F90:523-525` |
| `QRREF` stabilization | yes | not verified (not tested by `gpu_supported`) | not verified | not verified | `gpu_state_mod.F90` has no test for it |
| Equal-time observables, Hubbard | yes | device | device | device (`hk_*` kernels) | `gpu_obs_mod.F90:76-182,231` |
| Equal-time observables, other models | yes | host-assisted (host `Obser`, D->H of G per slice) | same | same | `main.F90:1033,1109` |
| Defect local observables `SpinZ_local`, `SpinZ_localtau` | yes | device (accumulators, D->H once per bin) | device | device | `gpu_obs_mod.F90:104-109,372`, `TIMEDEP_HOPPING.md` section 3 |
| Time-displaced observables (`Ltau = 1`), Hubbard | yes | device (tau pass and tau observables) | device (code); not run on the Max 1550 in the documents read | device (code); not run in the documents read | `gpu_tau_mod.F90:243-248`, `gpu_obs_mod.F90:179-183,592` |
| Defect (`ham_h0` /= 0) with `Ltau = 1` | pinning field switched off | pinning field switched off, then as above | same | same | `Hamiltonian_Hubbard_smod.F90:249-258` (all builds) |
| Defect with Projector, Langevin/HMC, or `LOBS_ST`/`LOBS_EN` /= 0 | stops the run | stops the run | stops the run | stops the run | `Hamiltonian_Hubbard_smod.F90:265-269` (all builds) |
| Custom (user) observables | yes | host-assisted (only for a model that passes `gpu_supported`; the host measures) | same | same | `main.F90:1033,1109` |
| Entanglement (Renyi) | yes | CPU code only; the predefined Hubbard model does not call it | same | same | `Predefined_Obs_mod.F90` (grep: not in `Hamiltonian_Hubbard_smod.F90`) |
| HDF5 output | yes | host side; used in production (H100 `data.h5`, A30 L = 16) | host side; production builds use HDF5 | not verified | `platforms/helma_h100.md` section 3, `platforms/intel_pvc_supermuc_ng2.md:15-24`, `README.md:23` |
| Rank-to-GPU binding | n/a | `CUDA_VISIBLE_DEVICES` per rank (`rank_wrap.sh`) | `ZE_FLAT_DEVICE_HIERARCHY=FLAT`, `ZE_AFFINITY_MASK` | not verified | `RUNTIME_SWITCHES.md` section 7 |

## Backend differences (code level)

| Item | NVIDIA | Intel | AMD | Evidence |
|---|---|---|---|---|
| `alf_la` routines without a device implementation (`LA_GUARD` stops the run) | none (not under `ALF_LA_OM`) | none | none | see below |
| Device QR | cuSOLVER (deterministic mode) | oneMKL; opt-in `ALFSH_DEVQR=1/2`, `ALFSH_DEVUNMQR=1` | own HIP QR and Q-apply by default; `ALF_HIP_QR=0` gives hipSOLVER | `alf_la_sycl_shim.cpp:119-148`, `alf_la_hip_mod.F90:518` |
| Async submatrix block path | no | opt-in `ALFSH_ASYNCBLK=1` | own `hk_blk_*` kernels, no switch | `gpu_lowrank_mod.F90:128,1334` |
| Small kernels (copies, gathers, observables) | OpenMP target | OpenMP target | HIP kernels (`alf_hip_kernels.hip.cpp`) | `gpu_obs_mod.F90:212-225` (`k_cpy`) |
| Deferred synchronisation | stream shared with OpenMP regions | queue of the OpenMP interop | `ALF_HIP_ASYNC`, default on only with the patched libgomp plugin | `RUNTIME_SWITCHES.md` section 5 |
| Test hooks `ALF_LA_TEST_*` | CUDA hooks | oneMKL hooks | cannot be armed | `RUNTIME_SWITCHES.md` sections 4, 5 |

`LA_GUARD` (`alf_la_mod.F90:117-123`) is defined only for `ALF_LA_OM` (Intel and AMD), and it heads only the host branch of 16 routines
(`alf_la_mod.F90`: lines 1523, 1553, 1586, 1614, 1647, 1676, 1710, 1738, 1773, 1804, 1841, 1869, 1902, 1934, 2013, 2095). Under `ALF_LA_OM` that branch
is preprocessed away: each of these routines calls its `om_*` twin, and `alf_la_om_mod` (oneMKL) and `alf_la_hip_mod` (HIP) both define every `om_*`
routine these 16 call. So in this tree no public `la_*` routine is without a device implementation. Not provided on any backend: `la_gemm_rc`,
`la_getri` (no ALF call site, `alf_la_mod.F90:100-102`).

## Run evidence (what has actually been run)

Source: `CONSOLIDATION.md` verification table (2026-11-13), platform pages.

| Backend | Machine, toolchain | Evidence | Scope |
|---|---|---|---|
| AMD | MI210, gfortran 14 GCN + HIP | GPU vs same-compiler CPU 1.15e-14 (L = 6), 5.0e-14 (L = 12); unified build bitwise equal to the pre-merge tip at L = 6, 12, 16, 24; also with `ALF_HIP_ASYNC=0`, `ALF_HIP_QR=0`. 2026-11-13 | Hubbard defect model, L = 6..24 |
| Intel | Max 1550 (SuperMUC-NG Phase 2), ifx 2025.3 / 2025.1 | GPU vs same-compiler CPU 2.2e-15 (L6F), 1.5e-14 (L12F), 5.5e-11 (L6T); V7 switches 1.36e-10 (L6T); unified build equal to the p5 build apart from 3.3e-16 in `Den_eq`, `SpinT_eq`, `SpinZ_eq`. 2026-11-13. Production chains and 8..32 rank throughput runs earlier (`platforms/intel_pvc_supermuc_ng2.md`) | Hubbard defect model, L = 6..24 |
| NVIDIA A30 | nvfortran 26.9 | Production runs and acceptance earlier (`BENCHMARKS.md`, `TIMEDEP_HOPPING.md` section 5, `MPI_GPU.md`, `PP_GPU.md`). The unified-commit comparison was **not run** (card busy, `CONSOLIDATION.md`) | Hubbard, defect, MPI, PP |
| NVIDIA H100 | nvhpc 25.7 (Helma) | L = 8 correctness 5.2e-13, restart, benchmarks to L = 32 (`platforms/helma_h100.md`, 2026-11-11). The row of the unified commit in `CONSOLIDATION.md` is still a placeholder | Hubbard defect model, PP, MPS |

Not run on any backend in the documents read: any model other than Hubbard on the device; `Ltau = 1` on the Max 1550 and the MI210; MPI on AMD;
the `Symm = F` and `Mz = T` combinations on AMD beyond what the L = 6..24 defect runs used (the defect runs use `Symm = T`, `N_FL = 2`).
Intel L6T (`Mz = T`, tolerance 1e-10) fails with the V7 switches (1.36e-10) and is noise-limited; see `ISSUE_INTEL_VS_A30.md` section 7.

See `NOT_IMPLEMENTED.md` for the refusal reasons and open gaps, `RUNTIME_SWITCHES.md` for switches, `PARKING.md` (repository root) for the old task list, and
`Known_issues.md` (repository root) for known problems.
