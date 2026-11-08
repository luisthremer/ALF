# Features: CPU vs NVIDIA GPU vs Intel GPU

Status 2026-11-06 (branch `gpu-port-delayed-submatrix`). Values:
- **supported**: implemented and validated (what "validated" means is given per row);
- **untested**: implemented or expected to work, not validated;
- **refused**: the GPU build refuses it at setup (`gpu_supported`, reason in `info`) and runs the CPU code path instead, never a silent wrong
  answer;
- **not implemented**: no GPU code, and the GPU path does not check for it (CPU-only feature that the model or run does not use).

GPU columns: NVIDIA = A30 (nvfortran, OpenMP offload + cuBLAS/cuSOLVER); Intel = Arc B580 (ifx, oneMKL). Ponte Vecchio, Helma and SuperMUC
have not been run. The GPU path is taken only if every check passes; otherwise the whole run uses the CPU path.

| Feature | CPU | NVIDIA GPU | Intel GPU | Notes |
|---|---|---|---|---|
| Finite-T DQMC | supported | supported (Hubbard) | supported on the B580 (Hubbard) | validated against ED (4-site chain, within 1σ) and CPU statistics |
| Projector (T=0) | supported | refused | refused | `tau_p`, CGRP |
| Hubbard (predefined) | supported | supported | supported on the B580 | Mz and `N_FL = 1`, square, honeycomb, bilayer, N-leg ladder tried |
| Other Hamiltonians | supported | refused | refused | Plain_Vanilla (flavour reconstruction); Kondo, tV, Spin_Peierls, Z2_Matter, `LR_Jz_Test` (non-diagonal vertices); LRC (field type 3) |
| HS fields type 1, 2 | supported | supported | supported | diagonal, rank-1 vertices only |
| HS fields type 3, 4 | supported | refused | refused | |
| Non-diagonal or rank > 1 vertices | supported | refused | refused | |
| Local sequential updates (rank-1) | supported | supported | supported | host decision, same random-number order |
| Delayed and submatrix updates | supported (CPU, rank up to 8, non-diagonal) | supported (rank-1 diagonal; test branch) | supported on the B580 | `Sequential_scheme`, `N_delay`; schemes agree bin for bin with standard; submatrix uses a host block of `G0(B,B)` |
| Cross-slice delay (`Lowrank_cross_slice`) | supported | refused | refused | |
| Global-in-τ and global space-time moves | supported | refused | refused | Hubbard itself has none |
| Langevin / HMC | supported | refused | refused | |
| MALA (sequential, global, global-τ) | supported | refused | refused | merged from master; MALA vertices use the standard path on the CPU |
| Parallel tempering | supported (MPI) | refused | refused | |
| MPI independent walkers | supported | supported (NVIDIA, tested 2026-11-08; rank mod number of devices) | refused (MPI builds) | one rank per GPU keeps the time per chain; more ranks per GPU share it, see [MPI_GPU.md](MPI_GPU.md) |
| Checkerboard hopping | supported | supported | supported | custom group kernels (and a dense matrix for Ndim <= 800, see `Known_issues.md`) |
| Dense hopping | supported | supported | supported | GEMM |
| Symmetric Trotter (`Symm`) | supported | supported | supported | `Hop_mod_Symm` on the device (dense A·M·B, see `Known_issues.md`) |
| Time-dependent hoppings (`Op_T` with `g_t`, `OpT_time_dependent`) | supported | supported on branch `gpu-delayed-submatrix-defect-fullgpu` (equal time, `Ltau = 0/1`, all update schemes) | supported (same) | the time-dependent elements must follow all static ones in `Op_T(:,nf)`; see `TIMEDEP_HOPPING.md` |
| Time-dependent couplings of `Op_V` (`g_t`) | supported | refused | refused | |
| Flavour symmetry reconstruction | supported | refused | refused | |
| Stabilization default | supported | supported | supported | pre-pivoted QR on the device |
| Stabilization LOG | supported | supported | supported | host QR with pivoting |
| STAB1/2/3, QRREF | supported | refused | refused | |
| `Propose_S0` | supported | refused | refused | |
| Equal-time observables (predefined) | supported | supported (Hubbard, on the device) | supported (Hubbard) | |
| Time-displaced observables (`Ltau = 1`) | supported | supported (Hubbard, on the device) | supported (Hubbard) | τ pass and τ observables on the device |
| Custom (user) observables | supported | not implemented | not implemented | only matters for models the GPU path refuses |
| Entanglement (Rényi) | supported | not implemented | not implemented | |
| HDF5 output | supported | untested | untested | host side, unchanged |
| OpenMP threading | BLAS threads only | – | – | |

## Not yet validated
- Ponte Vecchio (SuperMUC-NG Phase 2), Helma (H100), SuperMUC builds.
- Any model other than Hubbard on the device; anything beyond 62×62 at β = 5 (memory limit of a 24 GB card).
- Larger `Nwrap` (precision degrades quickly), HDF5 with GPU.

See `PARKING.md` for the open task list and `Known_issues.md` for the known problems.
