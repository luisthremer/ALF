# CG-ALF S1 on the CPU line: CGR determinant phase without LU (cpu-s1, 2026-11-26)

Same method as P2.4 step 1f on the GPU line (`gpu_stab_mod`), for the CPU `CGR` (CG-ALF T5 part 2).

## What changed
- `udv_state_mod`: every `UDV_State` keeps `det_U` = det(U)/|det U| and `det_ok`. `reset` sets 1; `decompose`
  takes it from the Householder reflectors of its QR (det H_i = −tau_i/conj(tau_i), 1 for tau_i = 0;
  `Householder_det_phase`) times the phase it puts into column 1 of U (pivot sign, R-diagonal phases); assignment
  and `MPI_Sendrecv` (tempering exchange) copy it. Modulus 1, no raw product, nothing accumulates over slices.
- `cgr1_mod::CGR` (default, STAB3, STABLOG): the phase of det(U_R^dagger U_L) is det_U(R) conj(det_U(L)); the
  ZGETRF of U_R^dagger U_L (`DET_C_PHASE`) is no longer called by default. STAB1/STAB2 keep their own LU (they need it
  for G); the projector (`CGRP`) is unchanged.
- Environment `ALF_CGR_LU` (read once by `CGR_work_alloc`): 0 tracked (default), 1 LU, 2 both — the tracked phase is
  used and compared with the LU at every CGR; a difference above 1e-8 stops the run; `CGR_work_dealloc` prints
  `CGR phase check (ALF_CGR_LU=2): CGR calls compared N, max |phase(tracked) - phase(LU)| X` on stdout (every rank).
  `ALF_CGR_FAULT=det` flips the tracked phase (negative control). A CGR with an untracked U stops the run
  ("untracked U in CGR"), never a silent LU fallback. The GPU path reads its own `ALF_GPU_CGR_LU` (`gpu_stab_mod`);
  each path reads only its own variable.
- Revert: the default is one commit ("cgr1_mod: determinant phase from tracked det(U), no LU"); `git revert` of it
  makes the LU the default again with the check mode and the tests intact (checked with `git merge-tree`).

## Verification
- Mode 2 on all 16 serial scenarios of `test_branch_parameters.json` in default, STABLOG and STAB3 builds (Hubbard
  variants, Kondo, tV, LRC, Z2_Matter, Spin_Peierls, global moves, Langevin/HMC): max |tracked − LU| 3.0e-15; projector
  scenarios and STAB1/STAB2 builds compare nothing (by design). MPI (2 ranks) and tempering (state exchange): 2.1e-15.
- Tests: `49-udv-det-phase` (reflector formula vs LU of Q with a tau = 0 column and a complex last reflector, negative
  control; decompose on both sides N = 4/50/200 vs `DET_C_PHASE`; N = 3200 diagonal vs the analytic phase) and
  `cgr-phase-e2e` (SU(2) with flux, Mz with chemical potential and delayed updates, Kondo: modes 0/1/2, mode 2
  byte-identical to mode 0, tracked vs LU within 1e-10, fault det and an invalid mode stop the run; in tempering builds
  the states move between ranks, which checks the det_U exchange — removing it makes the test fail).
- Commits 1–3 bitwise: 16 scenarios byte-identical to `Code-modularity-CPU`; commit 4 (default change): confout
  byte-identical in all 16 (same Markov chain), observables differ by at most 1.2e-13 (relative) in 6 of them.

## Timing (melmac, GNU -O3, 1 thread per run, pinned; same binary with `ALF_CGR_LU=1` vs `0`; Hubbard Mz, Beta 5,
Nwrap 10, Ltau 0, N_delay 16; median of 3; wall time per sweep including start-up)

Reference BLAS (the CI build):

| L | N | scheme | LU s/sweep | tracked s/sweep | change |
|---|---|---|---|---|---|
| 8 | 64 | standard | 0.0998 | 0.0956 | −4.2 % |
| 8 | 64 | delayed | 0.2066 | 0.2015 | −2.5 % |
| 12 | 144 | standard | 1.7930 | 1.7479 | −2.5 % |
| 12 | 144 | delayed | 0.9122 | 0.8894 | −2.5 % |
| 16 | 256 | standard | 8.9438 | 8.7787 | −1.8 % |
| 16 | 256 | delayed | 11.5383 | 11.0838 | −3.9 % |

MKL (gf_lp64, 1 thread):

| L | N | scheme | LU s/sweep | tracked s/sweep | change |
|---|---|---|---|---|---|
| 8 | 64 | standard | 0.0690 | 0.0654 | −5.2 % |
| 8 | 64 | delayed | 0.1249 | 0.1229 | −1.6 % |
| 12 | 144 | standard | 0.9400 | 0.9306 | −1.0 % |
| 12 | 144 | delayed | 0.3987 | 0.3924 | −1.6 % |
| 16 | 256 | standard | 5.9590 | 5.9162 | −0.7 % |
| 16 | 256 | delayed | 6.4498 | 6.3783 | −1.1 % |

On the CPU the LU is one N^3/3 step next to the CGR's own QR and GEMMs and the updates, so the gain is 1–5 % (on the
H100 it was −12.6 % at L = 32, where the LU was latency bound). The other gain is robustness: the default CGR no longer
forms an LU product, so the determinant-phase underflow class (N >~ 2950) cannot occur on this path.
Script: `docs/ai/phase2/cpu_s1_timing.py`.

## Follow-up
- Move the mode-2 summary from stdout to the info file once b-batch-2 (which edits `control_mod`) is merged.
