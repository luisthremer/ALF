# Not implemented, refused, and open gaps of the GPU path

Date 2026-11-13. Generated from the code at commit b7bf0bb3 (branch `gpu-port-unified`). The earlier text (2026-11-06, priority numbers of `PARKING.md`) is replaced.
What does run is in [FEATURES.md](FEATURES.md). Line numbers are anchors and drift; search for the reason string.

## 1. What a user sees

Two kinds of refusal exist. They behave differently.

**A. `gpu_supported` refuses (GPU builds only).** The run does not stop. The whole run uses the CPU code path. Setup prints, to stdout and `info`:
`GPU path (gpu_state): INACTIVE, reason: <string>` (`Prog/gpu_state_mod.F90:641`). If `ALF_GPU_SELFTEST` is also set, the run stops with that reason
(`gpu_state_mod.F90:642-646`). A run on a GPU build is therefore only on the GPU if `info` says `GPU path (gpu_state): ACTIVE`. Check this line.

**B. The Hubbard Hamiltonian refuses (every build, CPU included).** `Ham_Set` writes a message to stderr and stops with `ERROR_HAMILTONIAN`.
There is no fallback.

Order of the checks in `gpu_supported`: first match wins, so only one reason is shown even if several apply.

## 2. Refusals of `gpu_supported` (`Prog/gpu_state_mod.F90:499-617`)

| # | Reason string (`info`) | Condition | Line | Needed to lift it |
|---|---|---|---|---|
| 1 | `disabled by ALF_GPU_DISABLE=1` | environment `ALF_GPU_DISABLE=1` (allowed in production) | 513-515 | Not a gap. Unset the variable. |
| 2 | `tempering is not supported on the GPU path` | build has `TEMPERING` and not `PARALLEL_PARAMS` (compile-time test, so the whole build is inactive) | 520-522 | State exchange between replicas: host exchange plus device mirror of the exchanged `nsigma` and the Green functions. `PARALLEL_PARAMS` is supported (`PP_GPU.md`). |
| 3 | `STAB1/STAB2/STAB3 stabilization is not supported on the GPU path (default and LOG only, ticket 206)` | build has `STAB1`, `STAB2` or `STAB3` (compile-time) | 523-525 | Device versions of those stabilization variants. Use the default or `LOG` build. |
| 4 | `projector code is not supported on the GPU path` | `Projector = .T.` | 526-528 | Projector wraps, `Tau_p`, trial wave function propagation on the device (`tau_p`, CGRP). |
| 5 | `Langevin/HMC is not supported on the GPU path` | `get_Langevin()` or `get_HMC()` | 529-531 | Force evaluation and the leap-frog update on the device. |
| 6 | `global moves are not supported on the GPU path` | `get_Global_moves()`, `get_Global_tau_moves()` or `N_Global_tau > 0` | 532-534 | Global and global-in-tau proposals with a device `Green` update. Hubbard itself has no global moves. |
| 7 | `Propose_S0 is not supported on the GPU path` | `get_Propose_S0()` | 535-537 | Device (or host-assisted) `S0` proposal. |
| 8 | `MALA updates are not supported on the GPU path` | sequential, global or global-tau MALA, or `N_Global_tau_MALA > 0` | 538-540 | Force terms on the device. |
| 9 | `Lowrank_cross_slice (cross-slice delay) is not supported on the GPU path` | `get_Lowrank_cross_slice()` | 541-543 | Delayed state that survives a time-slice change on the device. |
| 10 | `Calc_Fl with weight reconstruction (ham%weight_reconstruction) is not supported on the GPU path` | `reconstruction_needed` | 544-546 | Flavour reconstruction of the weight and of G (`Calc_Fl`); would open Plain_Vanilla and similar models. |
| 11 | `no OpenMP offload device found` | `omp_get_num_devices() < 1` | 547-549 | Not a gap. Run on a node with a device. |
| 12 | `Hamiltonian <name>: Op_V(n,nf) is not diagonal` | any `Op_V` with `diag = .F.` | 555-559 | Non-diagonal vertices (Kondo, tV, Spin_Peierls, Z2_Matter, LR_Jz_Test). |
| 13 | `Hamiltonian <name>: Op_V(n,nf) has field type k (only 1 and 2 supported)` | field type not 1 or 2 | 560-564 | HS field types 3 and 4 (continuous fields, LRC). |
| 14 | `Hamiltonian <name>: field type of Op_V(n,nf) differs from flavour 1` | `type` differs between flavours | 565-568 | Per-flavour field type in the device tables (`V_type` is one per operator). |
| 15 | `Hamiltonian <name>: Op_V(n,nf) has N_non_zero = k (device update supports <= 1)` | `N_non_zero > 1` (rank > 1) | 569-574 | Rank-k updates on the device (ticket 205 implements rank 1 only). |
| 16 | `Hamiltonian <name>: Op_V(n,nf) has a time-dependent coupling g_t` | `Op_V%get_g_t_alloc()` | 575-580 | Time-dependent `Op_V` couplings in the update and wrap kernels. Time-dependent `Op_T` is supported. |
| 17 | `Hamiltonian <name>: hopping element (n,nf) is of unknown type` | `Op_T` element is neither real, complex nor `OpT_time_dependent` | 588-594 and 605-608 | Another `Op_T` container type in `Hop_mod_gpu_element`. |
| 18 | `time-dependent hopping: size(g_t) /= Ltrot` | `g_t` of a time-dependent `Op_T` has a length other than `Ltrot` | 596-598 | Not a gap; the model must give one `g_t` per slice. |
| 19 | `time-dependent hopping element too large (matrices per slice exceed 2 GB)` | `16*8*nh^2*Ltrot > 2e9` after support reduction (`td_reduce`) | 599-603 | Generate the per-slice matrices on the device instead of storing them. The point defect reduces to `nh` = 1. |
| 20 | `Hamiltonian <name>: static hopping element (n,nf) follows a time-dependent one (GPU path requires time-dependent elements last)` | a static element comes after a time-dependent one in `Op_T(:,nf)` | 609-612 | Order of the elements in the model, or mixed-order support in `gpu_wrap_mod`. |

Notes:
- Rows 12 to 16 are tested for every `Op_V(n,nf)` and are the only model-dependent refusals. The test is on structure, not on the Hamiltonian name.
  Which refusal a stock model other than Hubbard gets was not run in this audit (not verified); the old table (2026-11-06) said Kondo, tV, Spin_Peierls,
  Z2_Matter, LR_Jz_Test, Plain_Vanilla and LRC are all refused.
- Not tested by `gpu_supported` and therefore not refused: `QRREF`, `LOBS_ST`/`LOBS_EN` windows, `Ltau`, `Symm`, `N_SUN`, HDF5, the number of MPI ranks.
  Whether they work on each backend is in FEATURES.md ("not verified" where nothing documents it).
- After the GPU path is active, a few conditions stop the run instead of refusing: `Sequential_scheme` other than standard, delayed or submatrix
  (`gpu_lowrank_mod.F90:235`), pending sites above `N_delay` inside a block (`:652`), `N_Global_tau > 0` in the device update (`gpu_update_mod.F90:486,568`,
  unreachable after row 6), NaN scale with the default stabilization (`gpu_stab_mod.F90:589`, message "Switch to LOG is required").

## 3. Silent downgrades (GPU path active, part of it on the host)

These are not refusals. `info` shows `GPU path ... ACTIVE`, and the measurement stays on the host. Nothing is printed except the missing
line `GPU observables: device equal-time measurement ACTIVE` (`gpu_obs_mod.F90:241`).

| Condition | Where | Effect | Needed |
|---|---|---|---|
| Model is not the predefined Hubbard model | `gpu_obs_setup` is called only at `Hamiltonian_Hubbard_smod.F90:788` | host `Obser`, one D->H of G per measured slice (`main.F90:1033,1109`) | Device accumulators for other models; custom observables stay on the host by design. |
| `Projector_` | `gpu_obs_mod.F90:103` | host observables (the projector is refused earlier in any case) | none |
| Defect (`ham_h0` /= 0) with `N_FL` /= 2, or `Obs_local` missing or not of size 2 | `gpu_obs_mod.F90:104-109` | host observables | Device defect observables for `N_FL` = 1. |
| `N_FL` not 1 or 2, or `Obs_eq` of unexpected size, or `Ndim` /= `N*Norb` | `gpu_obs_mod.F90:110-113` | host observables | Generalised observable tables. |
| Lattice displacement table `imj` not a bijection | `gpu_obs_mod.F90:125` | host observables | Generic gather for such lattices. |
| `Ltau = 1` with the defect | `gpu_obs_mod.F90:179` (`.not. o_loc`) | tau observables on the host. Moot since 2026-11-13: with `Ltau = 1` the Hubbard model switches the defect off (`Hamiltonian_Hubbard_smod.F90:249-258`) | none |
| `ALF_GPU_NOOBS=1`, `ALF_GPU_NOUPDATE=1` | `gpu_obs_mod.F90:231`, `gpu_update_mod.F90:99-105` | host observables / host wraps | Debug switches (need `ALF_GPU_DEBUG=1`). |
| Stabilization `LOG` | `gpu_stab_mod.F90:857-858` | pivoted QR (`ZGEQP3`) on the host, N x N matrix goes D->H and the factors H->D | Device pivoted QR. |
| Submatrix scheme | `gpu_lowrank_mod.F90:33-55` | the k x k matrix `M` and its O(k^2) proposal ratio stay on the host (no device work on accept) | By design; the async block path (Intel) and `hk_blk_*` (AMD) shorten the cost. |
| Random numbers, accept/reject, field bookkeeping, bin I/O, MPI reduction | `main.F90`, `gpu_update_mod.F90` | host | By design. |

## 4. Refusals of the Hubbard Hamiltonian (every build, `Hamiltonian_Hubbard_smod.F90`)

| Message (stderr) | Condition | Line | What happens | Needed |
|---|---|---|---|---|
| `Ham_Set: the pinning field (ham_h0 /= 0) needs Projector = .F., sequential updates (no Langevin, no HMC) and the default measurement window (LOBS_ST = 0, LOBS_EN = 0).` | `ham_h0` /= 0 together with `Projector`, Langevin, HMC, `LOBS_ST` /= 0 or `LOBS_EN` /= 0 | 265-270 | Run stops, `ERROR_HAMILTONIAN` (rules of 2026-11-13, `CONSOLIDATION.md`) | A normalisation of `SpinZ_local`, `SpinZ_localtau` that counts every measurement and weights only visited end points; needs the same change in the host code, `gpu_obs_mod.F90` and the HIP kernel. |
| `Ham_Set: Set N_fl=2 if you want to use the pinning field` | `ham_h0` /= 0 and `N_FL` /= 2 | 295-298 | Run stops | none (model rule) |
| `Ham_Set: Ltau = 1, the pinning field ham_h0 is set to 0. Value in the parameter file: ...` | `Ltau = 1` and `ham_h0` /= 0 | 249-258 | Not a refusal: the field is set to 0, `info` gets `Pinning field : 0 (switched off because Ltau = 1)` (line 365) | Time-displaced observables with the defect (host `ObserT` needs different left and right slices; `gpu_symm(M,..,t1,t2)` exists, the device tau observables do not use it). |

## 5. Backend-specific gaps

### 5.1 `alf_la` (`LA_GUARD`)
`LA_GUARD` stops the run when a routine without device code is called under `ALF_LA_OM` (Intel, AMD) (`Libraries/Modules/alf_la_mod.F90:117-123`).
In this tree it never fires: the 16 sites (lines 1523, 1553, 1586, 1614, 1647, 1676, 1710, 1738, 1773, 1804, 1841, 1869, 1902, 1934, 2013, 2095) all sit in
the `#else` branch of `#if defined(ALF_LA_OM)`, and each routine has an `om_*` call in the `OM` branch. `alf_la_om_mod` (oneMKL) and `alf_la_hip_mod`
(HIP) define all of them. Routines not provided on any backend: `la_gemm_rc`, `la_getri` (no ALF call site, `alf_la_mod.F90:100-102`).
A new `la_*` routine without an `om_*` twin would stop an Intel or AMD run at its first call, not at setup.

### 5.2 Intel (`GPU=intel`)
- `libalfdevqr.so` and `libalfsb.so` are loaded by `dlopen` only if a switch asks for them; a missing library or symbol stops the run (`alf_la_sycl_shim.cpp:119-148`).
- Device QR falls back to oneMKL when the device limits are not met (large M).
- oneAPI 2025.3 compiles, but the OpenMP-to-SYCL interop crashes at run time with the Level Zero V2 adapter; 2025.1 works
  (`platforms/intel_pvc_supermuc_ng2.md`). Build without `omp_get_interop_*` in `omp_lib` needs the `bind(C)` declarations in `alf_la_mod.F90`.
- No MPS: ranks sharing a tile time-slice; one rank per tile.
- Documents disagree on the toolchain of the 2026-11-13 verification: `CONSOLIDATION.md` says ifx 2025.3, the platform page says only 2025.1 works at run time. Not resolved here (not verified).

### 5.3 AMD (`GPU=amd`)
- Needs GCC >= 14 with amdgcn offload; `-static` and `-ffast-math` are removed (`configure.sh:622-632`).
- The own HIP QR handles m <= 2304 only; larger matrices use hipSOLVER (`alf_la_hip_mod.F90:518-523`).
- `ALF_HIP_ASYNC` is on by default only if the patched libgomp GCN plugin registered its hook; without the patch every HIP call synchronises.
- `ALF_LA_TEST_*` hooks cannot be armed on HIP (declared, never called), so the negative controls of the `alf_la` order and info tests are not available there.
- MPI, PP and HDF5 on AMD: no run in the documents (not verified).
- The chunked `corr` kernel of the amd branch is not merged (`CONSOLIDATION.md`, `SPEEDUPS.md`).

### 5.4 NVIDIA (`GPU=nvidia`)
- nvhpc 24.11 builds but the L = 8 run aborts with `la_qrpp: NaN or Inf in the matrix`; use 25.7 (`platforms/helma_h100.md` section 2).
- cuSOLVER needs `cusolverDnSetDeterministicMode` (CUDA 12.6 or newer libraries, `cuda12.6` or `cuda12.9` with nvhpc 25.7).
- Several ranks on one GPU without the MPS daemon are 2.6 to 15 times slower per chain (`CLAUDE.md`, `platforms/helma_h100.md` section 4).
- `ALF_GPU_CC` exported at run time stops a production run (guard, `RUNTIME_SWITCHES.md` section 0).
- The unified-commit comparison on the A30 was not run; the H100 row of `CONSOLIDATION.md` is still a placeholder.

## 6. Opt-in switches not cleared for production (`RUNTIME_SWITCHES.md`)

`ALFSH_*`, `ALF_HIP_*` and `ALFDEVQR_*` are outside the production guard, so nothing stops a run that sets them. Check them by hand.

| Switch | Backend | Status | Source |
|---|---|---|---|
| `ALFSH_ASYNCBLK=1` | Intel | Bitwise identical to off in tests (V2). Release step "V2 next": same-seed check of the production binary first. | `ISSUE_INTEL_VS_A30.md` section 9 |
| `ALFSH_DEVQR=1`, `=2` | Intel | Rounding-level differences from oneMKL. Waits for the numerics review. `=2` also waits for the L6T rule and the tau-slot review. | same, sections 7 and 9 |
| `ALFSH_DEVUNMQR=1` | Intel | Same status as `DEVQR=2`. | same |
| `LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576` | Intel | "V1", recommended for production, bitwise identical in tests; confirm with one same-seed run. | same, section 9 |
| `ALFSH_QRNB`, `ALFSH_DEVQR_LIB`, `ALFDEVQR_*` | Intel | Tuning and test variables, no clearance. | `RUNTIME_SWITCHES.md` section 4 |
| `ALF_HIP_QR` (default on), `ALF_HIP_ASYNC` | AMD | Default on. Bitwise equal to hipSOLVER: not determined. Async output was bit-identical to sync mode. Acceptance run with both off and on at L = 6, 12 (`CONSOLIDATION.md`). | `RUNTIME_SWITCHES.md` section 5 |
| `ALF_GPU_TL`, `ALF_GPU_TR`, `ALF_GPU_VB` | all | Tuned on the A30. Bitwise independence of the value: not determined. Debug runs only. | `RUNTIME_SWITCHES.md` section 3 |
| `ALF_GPU_OBS_SKIP` | all | Skipped observables stay zero in the output. Never use for results. | same |
| `ALF_UNSTABLE_THRESHOLD` | all | Allowed in production (SuperMUC uses 1e-4). | same |

Open numerical item on Intel: L6T (`Mz = T`, same seeds as the CPU) is 1.36e-10 with the V7 switches against a tolerance of 1e-10. The criterion is noise-limited
(oneMKL itself fails 2 of 5 seed sets). A new pre-registered rule needs the owner's approval (`ISSUE_INTEL_VS_A30.md` section 7).

## 7. Performance and portability items still open

Taken from the documents, not from the code; their status was not re-checked here.
- Several chains per GPU (`MPI_GPU.md`: one rank per A30 recommended; aggregate throughput rises with ranks, per-chain time does not).
- Intel: about 76 000 small kernel launches per sweep, about 60 % of host time waiting for the GPU; further gains need fewer, larger kernels
  (`platforms/intel_pvc_supermuc_ng2.md`). Device `ZGETRF` and the observables phase are listed as not done (`ISSUE_INTEL_VS_A30.md` section 7).
- Memory footprint, about 1650 N^2 bytes (old text, not re-measured).
- Batched multi-walker rebuild: assessed, no-go recommended (`ISSUE_INTEL_VS_A30.md` section 8).
- Faulty GPU tile found in production on SuperMUC (`ISSUE_INTEL_VS_A30.md` section 10).
