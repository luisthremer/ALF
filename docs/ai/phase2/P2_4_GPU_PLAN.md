# P2.4 — GPU port onto Code-modularity-GPU: plan

Written 2026-11-26 by the manager session from three read-only surveys of `fork/gpu-port` (078d806b) and of
`Code-modularity` (CM). Numbers are `git diff --numstat` counts; nothing below was built unless marked "measured".
The roadmap section P2.4 in [`ROADMAP.md`](ROADMAP.md) stays the frame; this file says how to do it.

## 1. What there is to port

`gpu-port` is built on upstream before the refactor (merge base 521ddfca). It has three layers:

| Layer | Files (lines) | Vendor |
|---|---|---|
| Linear-algebra API `alf_la` (`la_gemm`, `la_geqrf`, `la_ormqr`, `la_getrf`, `la_trsm`, `la_qrpp`, ...; backend by macro `ALF_LA_HOST` / `_CUDA` / `_ONEMKL` / `_HIP`) | `Libraries/Modules/alf_la_mod.F90` (2115) | mixed: HOST common; CUDA inline (`bind(C)` to cuBLAS/cuSOLVER); ONEMKL via `!$omp dispatch` |
| Vendor shims | HIP: `alf_la_hip_mod` 834, `alf_hipk_mod` 559, `alf_la_hip_shim.cpp` 135, `alf_hip_kernels.hip.cpp` 603, `alf_hip_qr.hip.cpp` 349. SYCL: `alf_la_sycl_shim.cpp` 304, `_blk.cpp` 102, `_devqr.cpp/.hpp` 444, `_qr2.hpp` 456 | AMD / Intel only |
| Device code (OpenMP 5 target offload, one conservative subset for nvfortran, ifx, gfortran-amdgcn; ADR-001) | `Prog/gpu_state_mod` 1902 (device mirror, `gpu_supported`, env guard), `gpu_wrap_mod` 2689, `gpu_update_mod` 633, `gpu_lowrank_mod` 1727, `gpu_stab_mod` 2122, `gpu_obs_mod` 729, `gpu_tau_mod` 2634 = 12 436 | common, with small `#if ALF_LA_HIP` / `ALF_LA_ONEMKL` branches (HIP kernels replace some OpenMP loops on AMD) |
| Hooks into upstream routines | `main.F90` +377, `Wrapgr_mod` +310/-44, `tau_m_mod` +10, Hubbard `Obser` +16, `Hop_mod` +93 (read accessors), `Operator_mod` (`Op_get_E_exp`), `control_mod` (`Control_Precision*_vals`), `configure.sh` +108, Makefiles | common, behind `ALF_GPU` |

Already on gpu-port, so not separate work: MPI on GPU, `PARALLEL_PARAMS`, delayed/submatrix on the device.
Not on the GPU (`gpu_supported` refuses and the run falls back to the CPU path): tempering without PP, STAB1/2/3,
projector, Langevin/HMC/MALA, global moves, `Propose_S0`, non-diagonal or rank > 1 `Op_V`, field types other than 1/2,
time-dependent `g_t` (20 refusals, `docs/gpu/NOT_IMPLEMENTED.md`). These are P2.5/P2.6.

Branches with work that is not on gpu-port (three-dot diffs): `fix/det-phase-underflow` (GPU determinant phase,
class A, host part common, one HIP line), `speedup/nvidia` (S1 phase without LU, looks vendor-neutral but measured on
CUDA only; S2 CUDA-only), `intel-batched` and tag `stp/production-mpi-intel-large-prompt-team-combined` (SYCL batched
update, device QR, Intel only), the amd production branch (`alf_exit`, HIP det-phase line), `gpu-port-unified` (the
point-defect model: model work, not vendor work, P2.6).

## 2. Where it attaches in Code-modularity

The refactor moved the code that gpu-port hooks into. The hooks must be re-homed, not merged:

| gpu-port hook | CM place |
|---|---|
| `main.F90` setup: `gpu_state_setup`, `gpu_stab_wrapul/store_l/cgr` | `qmc_setup_mod` / `qmc_walker_mod` (`Walker_alloc` :125, `Walker_init_stack`, `Walker_init_GR`, `Recompute_GR_phase` :229) |
| `main.F90` sweeps: `gpu_stab_wrapur/wrapul/exchange/cgr`, `gpu_obs_*`, `gpu_g_host` | `qmc_sweep_mod` (`Sweep_up` :75, `Sweep_down` :135) |
| `Wrapgr_mod` early `gpu_wrapgrup/do` | same file (`WRAPGRUP` :87, `WRAPGRDO` :253) |
| `tau_m_mod` early `gpu_tau_pass` | same file (`TAU_M` :84) |
| Hubbard `Obser` / `Alloc_obs` | same file (`Hamiltonians/Hamiltonian_Hubbard_smod.F90`) |
| device data from the walker | `qmc_walker_t` fields `GR`, `udvl/udvr`, `udvst`, `Stab_nt` (qmc_walker_mod :74-87) |
| `gpu_lowrank_mod` on gpu-port's `Lowrank_update_mod` (2506 lines) | CM's own `lowrank_update_mod` (1835 lines, different interfaces; no `Lowrank_cross_slice`) |

CM's `linalg_mod` (16 `la_z*` host wrappers, callers in cgr1, cgr2_2, udv_state, upgrade, lowrank, ...) is not enough
as the only seam: the element-wise steps between the BLAS calls (pivoting, D scaling, phase loops) and the wraps,
updates and measurements are bespoke kernels in the `gpu_*` modules, and `alf_la` has no device zhemm, ztrmm, zgemv,
zscal, zcopy, zlacpy. See decision D1.

CI rules the gpu-port code breaks today and that the port must fix: bare `USE` in all seven `gpu_*` files (CM needs
`ONLY:`), direct LAPACK calls in `Prog/` (`gpu_tau_mod` 20 lines: route through `linalg_mod`), implicit `SAVE` not checked.
Every CI job builds with GPU=none, so each `#ifdef ALF_GPU` block must compile out cleanly in all STAB/MPI variants.

## 3. Steps (each a dev/canvas feature into `Code-modularity-GPU`, merged one after another)

The common part is most of the code, so the "scaffold" of the roadmap is split into five small features. Vendor
code comes only in step 2.

| Step | Feature | Content | Done when |
|---|---|---|---|
| 1a | `gpu-build` | `configure.sh GPU=none\|host\|nvidia\|intel\|amd` (only none and host active; the others stop with "not yet on this line"); Makefile switch with stamp file. No new linear-algebra module (D1) | GPU=none bitwise identical (16 scenarios + preprocessor identity of every changed file); CI green |
| 1b | `gpu-state` | Hubbard only (D3): `gpu_state_mod` on `qmc_walker_t` (device mirror, `gpu_supported`, `ALF_GPU_DISABLE`, env guard), hooks in setup/walker, read accessors in Hop_mod/Operator_mod | as 1a; GPU=host run: the mirror round-trips G bitwise |
| 1c | `gpu-wrap-update` | `gpu_wrap_mod`, `gpu_update_mod`, `gpu_lowrank_mod` adapted to CM's `lowrank_update_mod`; Wrapgr hooks; their `alf_la` calls become `linalg_mod` calls (missing wrappers added there with a case in test 13-linalg) | as 1a; GPU=host vs CPU trajectories agree (configurations byte-identical where the algorithm is the same) |
| 1d | `gpu-stab` | `gpu_stab_mod` with the GPU determinant-phase fix of `fix/det-phase-underflow` (class A, must come with it) | as 1c, plus the N >~ 2950 phase test |
| 1e | `gpu-obs-tau` | `gpu_obs_mod`, `gpu_tau_mod` (LAPACK through `linalg_mod`), hooks in `tau_m` and Hubbard `Obser` (Hubbard only, D3) | as 1c for observables and Ltau = 1 |
| 1f | `gpu-s1` | Speedup S1 (owner yes, 2026-11-26): determinant phase from the tracked det(U) of the Householder reflectors instead of an LU (`rd/cgalf-s1`), no LU; documented (why no LU; invariant: U stays the exact product of the tracked reflectors; check mode); the same approach for every vendor where possible (Q11). Own feature after 1d (scheduled after 1e, which was already running) | as 1d; phase identical to the LU path where both are defined; check mode passes |
| 2 | `gpu-nvidia`, `gpu-intel`, `gpu-amd` (parallel dev branches from the scaffold) | device backend inside `linalg_mod` (D1): CUDA (cuBLAS/cuSOLVER), ONEMKL + SYCL shims, HIP + kernels; GPU=none keeps the plain calls. Per-vendor extras later (S2, intel-batched). `gpu-intel`: switches `ALFSH_OBSCORR=1`, `ALFSH_BLKFUSE=5`, `ALFSH_BLKPF=1` on by default after their own check on Max 1550, active values written to the `info` file (owner, Q4b); investigate the `Mz = T` deviation (L6 1.36e-10, L12 1.21e-9) with claude-bughunter (Q4c) | builds and the trajectory checks of `docs/gpu/VERIFICATION.md` on A30/H100, B580/Max 1550, MI210 |
| 3-4 | integrate, consolidate | as in ROADMAP P2.4 | as in ROADMAP |

Not in P2.4 (each its own feature later, with its own check): the CPU-path changes gpu-port made without a guard
(D4: queue in ROADMAP.md "CPU-line queue"; `QR_PREPIVOT` variant 1 only as a test-only flag in 1d/1e), the test
Hamiltonian `LR_Jz_Test`, the defect model, tempering on the GPU.

## 4. Decisions (owner, 2026-11-26)

- **D1 — decided: the backend dispatch lives in Fakher's `linalg_mod`** (as ROADMAP P2.4 says). No separate
  `alf_la_mod`. Consequences: the `gpu_*` code calls the `la_z*` wrappers; wrappers gpu-port needs and CM lacks
  (column-norm pre-pivoted QR `la_qrpp`, batched GEMM, real x complex GEMM, LWORK queries) are added to `linalg_mod`
  with a test case in `testsuite/matmod.tests/13-linalg` when the module that needs them is ported; the device
  backend (device addresses, workspaces, asynchronous error checks) goes into `linalg_mod` in step 2, behind the
  vendor macros, while GPU=none keeps the plain calls bitwise. Fakher's 6G.E1-E7 (merged into the CPU line
  2026-11-26) already route cgr1, cgr2_2, udv_state, upgrade, wrapul/wrapur and mymats through `linalg_mod`.
- **D2 — decided: yes**, a `GPU=host` test build (step 1a) runs the GPU code path on the CPU for tests on the Mac
  and in CI.
- **D3 — decided: yes**, step 1 (1b-1e) covers the Hubbard model only (marked in the table above); other models
  keep the CPU path through `gpu_supported`, as on gpu-port. Other models: P2.6.
- **D4 — decided: option (b)** (interview 2026-11-26): gpu-port's 8 unguarded CPU-path changes come as separate
  CPU-line features later (ROADMAP.md "CPU-line queue"), so GPU=none stays byte-identical to `Code-modularity-CPU`.
  Exceptions: `QR_PREPIVOT` as a test-only flag (default off) in 1d/1e; `ALF_UNSTABLE_THRESHOLD` first in the queue.
- **S1 — decided: yes** (Q3, 2026-11-26): step 1f above. CG-ALF checks that tracking is right in general (T5).

## 5. Probes

- **GPU=host build of gpu-port (measured 2026-11-26, Mac arm64, gfortran, noMPI, `-DALF_GPU`, HOST backend, no
  offload target).** Everything compiles. Three things block a run, all small: (1) the link needs a LAPACK with
  `zgeqrt`/`zgemqrt` (Homebrew `lapack`; the system one lacks them; used by `gpu_tau_mod`); (2) the env guard
  `gpu_env_read` stops when `/proc/self/environ` is missing (Linux only); (3) `gpu_supported` refuses when
  `omp_get_num_devices() < 1`, and `omp_set_default_device(Mod(rank, 0))` would divide by zero. With (2) and (3)
  patched in a scratch copy the GPU path ran ACTIVE on the host (mirror, wraps, updates, stabilisation, tau pass,
  device observables): Hubbard 4x4, beta 5, 5 sweeps x 2 bins, default parameters. Result vs the CPU path of the same
  binary (same seeds): `confout_0` byte-identical; `Ener_scal`, `Kin_scal`, `Pot_scal` agree to about 1e-11 relative
  (pre-pivoted QR instead of QRP); Precision Green max 1.8e-8. Run time 34 s against 0.13 s (OpenMP host fallback),
  so only small sizes for tests. D2 is feasible: step 1a adds `GPU=host`, a host-device mode in `gpu_supported`
  (`omp_get_num_devices() == 0` allowed only in that build) and a macOS-safe env reader.
