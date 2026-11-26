# CG-ALF T5 part 1: is the tracked determinant phase right in general? (2026-11-26)

Question (owner, with the yes to S1): does the CGR phase from the tracked det(U) (P2.4 step 1f, merge a3e1220f on
`Code-modularity-GPU`) agree with the LU phase in every case the GPU path accepts?

Answer: **yes, in every accepted case tested.** The largest difference |phase(tracked) − phase(LU)| over all runs is
**8.1e-15 (L = 56, N = 3136)** (rounding level; the check mode aborts above 1e-8). Real-phase models (Mz-Hubbard without flux, attractive
SU(2)) agree **exactly** (0). No case broke the invariant; no class-A candidate. The negative control (`ALF_GPU_FAULT_STAB=det`)
is stopped at the first comparison in both stabilisation builds.

## Setup
- Code: `Code-modularity-GPU` 889f4633 (contains the 1f merge a3e1220f). Host: melmac, GNU, `GPU=host` (OpenMP host
  fallback: the device code path runs on the host, so this checks the algorithm and the invariant, not vendor kernels).
- Builds (`. ./configure.sh GNU <MODE> [LOG] GPU=host`): `def` noMPI, `log` noMPI LOG, `mpi` MPI, `pp` PARALLEL_PARAMS,
  `temp` Tempering — reference BLAS/LAPACK; `defmkl`, `logmkl` the same as def/log linked against MKL (gf_lp64, gnu_thread)
  for L >= 24 (the reference BLAS is too slow there).
- Every run: `ALF_GPU_DEBUG=1 ALF_GPU_CGR_LU=2` (tracked phase used, LU phase computed and compared at every CGR; the info
  file reports the number of comparisons and the max difference). Parameters: `Scripts_and_Parameters_files/Start/parameters`
  with the edits in the table (default L1 = L2 = 4, NSweep = 6, NBin = 2, Beta = 5, Nwrap = 10, Ltau = 1, Mz = .T. N_FL = 2;
  "su2" = Mz = .F., N_FL = 1).
- Driver: `docs/ai/phase2/t5/t5run.py <build> <case> KEY=VAL ... [--np N] [--dirs Temp_0,Temp_1] [--restart] [--env K=V]
  [--omp N]` (`--restart`: second run from `confin_*` = the first run's `confout_*`); `t5_scen.py <build>` runs all
  serial scenarios of `testsuite/test_branch_parameters.json` in check mode.

## Results (CGR comparisons / max |tracked − LU|)

### Stabilisation variants, schemes, Nwrap, N_FL, restart (L = 4 unless noted)
| case | edits | def | log (STABLOG) |
|---|---|---|---|
| mz-std | Mz, N_FL 2 | 240 / 0 | 240 / 0 |
| su2-std | N_FL 1 | 120 / 1.4e-15 | 120 / 1.3e-15 |
| mz-delayed-N4 | delayed, N_delay 4 | 240 / 0 | 240 / 0 |
| mz-submatrix-N4 | submatrix, N_delay 4 | 240 / 0 | 240 / 0 |
| su2-delayed-N4 | | 120 / 1.4e-15 | 120 / 1.3e-15 |
| su2-submatrix-N4 | | 120 / 1.4e-15 | 120 / 1.3e-15 |
| mz-ltau0 | Ltau 0 | 240 / 0 | 240 / 0 |
| su2-nwrap2 / 25 / 50 | Nwrap 2, 25, 50 | 600 / 1.6e-15; 48 / 1.4e-15; 24 / 9.5e-16 | 600 / 1.8e-15; 48 / 1.4e-15; 24 / 9.5e-16 |
| mz-nwrap50 | Nwrap 50 | aborted "unstable" (control_mod) — the CPU path aborts the same way (`ALF_GPU_DISABLE=1`), not a phase issue | same |
| mz-h0 | ham_h0 1 | 240 / 0 | 240 / 0 |
| mz-bilayer-hc | Bilayer_honeycomb 2x2, U2 = 2 | 240 / 0 | 240 / 0 |
| su2-symmF-cbF | Symm, Checkerboard .F. | 120 / 1.3e-15 | 120 / 1.2e-15 |
| mz-beta20 | Beta 20 | 960 / 0 | 960 / 0 |
| su2-beta20-U8 | Beta 20, U 8 | 480 / 2.3e-15 | 480 / 1.5e-15 |
| L8-mz | L 8 | 240 / 0 | 240 / 0 |
| L12-su2-delayed | L 12, N_delay 16, 2 sweeps | 20 / 3.4e-15 | 20 / 4.4e-15 |
| L16-mz-delayed | L 16, N_delay 16, 1 sweep | 20 / 0 | 20 / 0 |
| restart-mz (2nd run) | from confin_0 | 240 / 0 | 240 / 0 |
| restart-su2-delayed (2nd run) | from confin_0 | 120 / 1.5e-15 | 120 / 1.4e-15 |
| restart-sign-nphi1 (2nd run) | from confin_0 | 240 / 1.2e-15 | 240 / 1.2e-15 |
| scen-Hubbard_N_Leg_Ladder | test_branch scenario | 800 / 1.3e-15 | 800 / 1.4e-15 |
| scen-Hubbard_Bilayer | test_branch scenario | 400 / 2.4e-15 | 400 / 2.4e-15 |

### Sign / phase problem (all GPU-accepted)
| case | edits | def | log |
|---|---|---|---|
| sign-mz-chem | ham_chem 0.5 (real sign problem) | 240 / 0 | 240 / 0 |
| sign-mz-chem-beta10 | ham_chem 0.5, Beta 10, L 6 | 480 / 0 | 480 / 0 |
| sign-mz-hc-chem | Honeycomb 3x3, ham_chem 0.4 | 240 / 0 | 240 / 0 |
| sign-mz-nphi1 | N_Phi 1, ham_chem 0.3 (complex phase, N_FL 2) | 240 / 1.4e-15 | 240 / 1.5e-15 |
| sign-mz-twist | Phi_X 0.3, Phi_Y 0.17, Bulk .F., chem 0.2 | 240 / 1.3e-15 | 240 / 1.4e-15 |
| sign-su2-nphi1 | N_FL 1, N_Phi 1 | 120 / 1.2e-15 | 120 / 1.5e-15 |
| sign-su2-twist-chem | Phi_X 0.3, chem 0.4 | 120 / 1.3e-15 | 120 / 1.2e-15 |
| attr-su2-Uneg | U = −4 | 120 / 0 | 120 / 0 |

### Large N (MKL builds; Symm .F., Beta 2, delayed N_delay 32, N_Phi 1, ham_chem 0.3: complex phases)
| case | N = Ndim | sweeps | defmkl | logmkl |
|---|---|---|---|---|
| L24 Mz | 576 | 4 | 32 / 4.5e-15 | 32 / 3.1e-15 |
| L32 SU(2) | 1024 | 3 | 12 / 7.8e-15 | |
| L40 Mz / SU(2) | 1600 | 2 | Mz: 16 / 6.0e-15 | SU(2): 8 / 6.2e-15 |
| L56 Mz (N > ~2950 class) | 3136 | 1 | 8 / 6.7e-15 | 8 / 6.7e-15 |
| L56 SU(2) | 3136 | 1 | 4 / 8.1e-15 | |

### MPI, parallel parameters, tempering
| case | build | result |
|---|---|---|
| mpi-mz-std, np 2 | mpi | 240 / 0 |
| mpi-su2-delayed, np 2 | mpi | 120 / 1.5e-15 |
| mpi-sign-nphi1, np 2, + restart | mpi | 240 / 1.3e-15; restart 240 / 1.7e-15 |
| pp-mz-std, np 4, 2 groups (Temp_0/1), + restart | pp | each group 240 / 0, restart 240 / 0 |
| pp-sign-su2-nphi1, np 2, 2 groups | pp | 120 / 1.5e-15, 120 / 1.2e-15 |
| temp-mz-std, np 2 | temp | INACTIVE (refused: "tempering is not supported on the GPU path") — no device exchange to check |

### Negative control
| case | def | log |
|---|---|---|
| fault-det (Mz), fault-det-su2 (`ALF_GPU_FAULT_STAB=det`, sign of the tracked phase flipped) | stopped at comparison 1, diff 2.0 | same |

## Refused by `gpu_supported` (no device phase to check; the run keeps the CPU path)
STAB1/STAB2/STAB3 (compile time), Tempering (exchange; PP is accepted), projector (Hubbard_PAM*, Plain_Vanilla scenarios),
Langevin/HMC, global and global-tau moves (LRC, Z2_Matter, Hubbard_Global_moves scenarios), Propose_S0, MALA,
`Sequential_scheme = "auto"`, N_delay = −1, weight reconstruction (Hubbard_Plain_Vanilla), non-diagonal Op_V (Kondo, tV,
Spin_Peierls), field type 3 (Hubbard `Continuous = .T.`), time-dependent couplings. So the sign-problem models reachable on the
GPU path are Hubbard variants (chemical potential, flux, twist, honeycomb); all were checked above.

## Notes
- The invariant held through every U writer exercised here: wraps + decompose (default and STABLOG), resets, the
  initial upload (`gpu_stab_put_udv`, first CGR of every run and every restart from confin), the copies of the backward
  sweep, and MPI/PP (one device mirror per rank).
- Max difference grows slowly with N (1e-15 at N = 16, 8e-15 at N = 1024, 8.1e-15 at N = 3136), 6–7 orders below the abort limit.
- Not covered here: vendor backends (P2.4 step 2: each vendor runs `ALF_GPU_CGR_LU=2` in its verification, design note §7);
  long production runs (rounding-level phase differences may flip a Metropolis decision after many sweeps; not a correctness issue).
- Raw outputs: melmac `~/work/ALF_P2/t5/` (`runs/<build>/<case>/`, `*results*.txt`).
