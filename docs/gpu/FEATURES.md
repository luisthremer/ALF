# Features: CPU vs NVIDIA GPU vs Intel GPU

Rows from R1 §9, checked against `Prog/` where noted. Values:
- **supported**: implemented and validated;
- **untested**: should work, not validated;
- **planned (Pn)**: in the plan for phase n;
- **not supported**: GPU builds refuse it at setup or fall back to the CPU code path, never a
  silent wrong answer;
- **out of scope**.

Status 2026-10-27 (end of Phase 0): no GPU code exists yet. CPU column = upstream; validated on
nividia means unit tests + test_vs_ed + Hubbard references (BASELINE.md).

| Feature | CPU | NVIDIA GPU | Intel GPU | Notes |
|---|---|---|---|---|
| Finite-T DQMC | supported (validated on nividia, GNU + ifx) | planned (P2) | planned (P3) | |
| Projector (T=0) | supported | out of scope | out of scope | `tau_p`, CGRP |
| Hubbard (predefined) | supported (validated) | planned (P2) | planned (P3) | MVP model |
| Other Hamiltonians (Plain_Vanilla, Kondo, tV, LRC, Z2_Matter, Spin_Peierls) | supported | not supported | not supported | need generic Op_V/observable support |
| HS fields: type 1 Ising, type 2 discrete 4-value | supported | planned (P2: type 2 as used by Hubbard Mz; type 1 untested) | planned (P3) | diagonal Op_V only in MVP |
| HS fields: type 3 continuous, type 4 complex | supported | not supported | not supported | |
| Non-diagonal Op_V (op_dim > 1) | supported | not supported (MVP) | not supported | needed for tV etc.; fits delayed/submatrix scheme |
| Local sequential updates | supported | planned (P2, rank-1) | planned (P3) | delayed/submatrix later (401–403) |
| Delayed / submatrix updates (Sun & Xu) | not in upstream | later (403) | later (403) | CPU first (401, 402) |
| Global-in-τ moves | supported | not supported | not supported | |
| Global space-time moves | supported | out of scope | out of scope | |
| Langevin / HMC | supported | out of scope | out of scope | |
| Parallel tempering / parallel params | supported (MPI) | out of scope (future: GPU-aware MPI) | out of scope | |
| MPI independent walkers | supported | future (405) | future (405) | 1 rank per GPU / PVC stack |
| Checkerboard hopping | supported (validated) | planned (P2) | planned (P3) | custom family kernel |
| Dense hopping | supported (validated) | planned (P2) | planned (P3) | GEMM |
| Symmetric Trotter (`Symm`) | supported | untested → planned after P2 | same | extra Hop_mod_Symm products; test_vs_ed uses Symm=T, so needed for GPU test_vs_ed |
| Time-dependent couplings | supported | not supported | not supported | |
| Flavour symmetry reconstruction | supported | untested | untested | |
| Stabilization default | supported (validated) | planned (P2) | planned (P3) | QRCP replaced by pre-pivoted QR on GPU, pending ticket 104 |
| Stabilization STABLOG | supported (validated) | planned (P2) | planned (P3) | |
| STAB1/2/3, QRREF | supported | not supported | not supported | |
| Equal-time observables (predefined) | supported (validated) | planned (P2: host copy first, then device) | planned (P3) | |
| Time-displaced observables (Ltau=1) | supported (validated) | planned (P2) | planned (P3) | |
| Custom (user) observables | supported | planned via host copy of G | same | |
| Entanglement (Rényi) | supported | not supported | not supported | |
| HDF5 output | supported (used for pyALF) | planned (host side, unchanged) | planned | |
| OpenMP threading | BLAS threads only | – | – | |

## Future work
- **MPI and multi-GPU:** one rank per GPU (NVIDIA) or per PVC stack (8 per SuperMUC node).
  GPU-aware MPI for tempering, which exchanges G, `udvst`, udvl/udvr and nsigma. Not testable on
  nividia (one GPU); needs SuperMUC.
- **Global updates, Langevin/HMC, projector:** out of scope for now. Langevin/HMC do not use local
  updates, so they need a separate port (forces, CG-type solvers).
- **More Hamiltonians:** generic non-diagonal Op_V and device observables.
- **Delayed and submatrix updates:** PLAN.md 401–403.
