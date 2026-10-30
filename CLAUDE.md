# CLAUDE.md — ALF delayed/submatrix local updates (branch `feature/delayed-submatrix`)

Research branch: delayed and submatrix variants of ALF's sequential local update.
Priorities: Markov-chain correctness > default behaviour unchanged > speed > elegance.
Never push or open merge requests. Do not install software; ask the user.

## Layout
- Repo: `ALF/` (cloned from `git@github-work:ALF-QMC/ALF.git`; the GitLab URL needs a login).
  Baseline commit: `f81b5cbc` (master, 2026-08-20).
- Papers (PDFs, not in git): `../papers/`.
- Notes: `docs/algorithm-notes.md`; design/plan: `docs/superpowers/specs|plans/`.
- Local run dirs `Run/`, `runs/` are excluded via `.git/info/exclude`.

## Environment
- Apple M1 Pro (10 cores, 32 GB), macOS 27.2, GNU Fortran (Homebrew GCC) 16.2.0, CMake 4.4.3.
- BLAS/LAPACK: Apple Accelerate (vecLib), linked through `-llapack -lblas`. OpenBLAS is present at
  `/opt/homebrew/opt/openblas` (unused so far). ALF uses no complex-valued BLAS functions
  (only DZNRM2 → real), so Accelerate's f2c return convention is not a problem.
- Default flags (configure.sh GNU): `-std=f2008 -cpp -O3 -ffree-line-length-none -ffast-math -fopenmp`.
- python3 + numpy 2.5.3 available; **pyALF and scipy are not installed** (ALF's `test_branch`/`test_vs_ed`
  harnesses need pyALF). Use the Fortran `Analysis/ana.out` for analysis.
- Shell gotcha: `rm` is aliased to `trash`; use `/bin/rm` in scripts. zsh globs: quote `--include='*.F90'`.

## Commands
Environment variables do not persist between shell calls, so always source configure.sh in the same command:
```bash
cd ALF && source configure.sh GNU noMPI NO-INTERACTIVE >/dev/null 2>&1 && make all -j8      # ~20 s
cd ALF && source configure.sh GNU noMPI NO-INTERACTIVE >/dev/null 2>&1 && make program -j8  # ALF.out only
```
Run the Hubbard example (6x6, U=4, β=5, 20 sweeps × 5 bins, ~5 s):
```bash
cd ALF && cp -r Scripts_and_Parameters_files/Start Run && cd Run && ../Prog/ALF.out && ../Analysis/ana.out Ener_scal
```
ALF test suite (CTest; 39 tests, all pass at baseline):
```bash
cd ALF && source configure.sh GNU noMPI NO-INTERACTIVE >/dev/null 2>&1 && mkdir -p testsuite/tests && cd testsuite/tests && cmake .. && make -j8 && ctest -j8
```

## Baseline (commit f81b5cbc, Start example)
- Precision Green Mean/Max 3.24e-11 / 1.35e-7; Phase 4.55e-15; tau 7.71e-12 / 7.02e-8;
  acceptance 0.4296; CPU 5.39 s; Ener_scalJ: −29.822 ± 0.130.

## Code map: sequential local update
- `Prog/main.F90`: the sweep. Up: `DO NTAU=0,LTROT-1 → WRAPGRUP`; at `Stab_nt` → `WRAPUR` + `CGR`
  (+ `Control_PrecisionG/P`, `Op_phase`); measurements via `ham%Obser` (after `Hop_mod_Symm` if Symm).
  Down: `WRAPGRDO`, measurements, `WRAPUL` + `CGR` at stabilization; `Tau_m`/`Tau_p` for time-displaced.
- `Prog/Wrapgr_mod.F90`: `WRAPGRUP` / `WRAPGRDO` contain the whole per-slice sequential update.
  Per vertex `n`: `Op_Wrapup/Op_Wrapdo` (type 1 = rotate + propagate, type 2 = rotate back),
  `nsigma%flip`, `ham%S0`, the `T0_proposal` ranf, `Upgrade2`. `Wrapgr_Random_update`/`PlaceGR` = Global_tau moves.
- `Prog/upgrade_mod.F90`: `Upgrade2`: ratio (Sylvester, all flavours, `^N_SUN`, `e^{Z1 α}`, γ ratio),
  Metropolis `Weight = S0·T0·|Re(Phase·R)/Re(Phase)|` vs `ranf_wrap()`, phase update, rank-N_P G update.
- `Prog/Operator_mod.F90`: `Operator` type (`N, N_non_zero, diag, U, E, P, g, g_t, alpha, type`), `Op_Wrapup/do`, `Op_phase`.
- `Prog/Hop_mod.F90`: hopping propagation (`HOP_MOD_mmthr`, `_mmthl_m1`, …).
- `Prog/QMC_runtime_var_mod.F90`: `&VAR_QMC` namelist, getters/setters, MPI broadcast.
- Hubbard: `Prog/Hamiltonians/Hamiltonian_Hubbard_smod.F90`; `Mz=.T.` → `Predefined_Int_U_MZ`
  (N_FL=2, N=1, diag, type 2, g=±sqrt(Δτ U/2), α=−1/2); `Mz=.F.` → `Predefined_Int_U_SUN` (N_FL=1).
- G convention: `G = (1+B)^{-1}`, `G_xy = <c_x c_y^†>`; `B_τ = e^{V_M}⋯e^{V_1} e^{-ΔτT}`.

## Papers (all verified on arXiv/OpenAlex 2026-10-30)
- [SX-S] F. Sun, X. Y. Xu, "Boosting determinant quantum Monte Carlo with submatrix updates: Unveiling
  the phase diagram of the 3D Hubbard model", SciPost Phys. 18, 055 (2025), doi:10.21468/SciPostPhys.18.2.055,
  arXiv:2404.09989. **This is "Xu's submatrix update" for DQMC.** (OpenAlex lists a wrong second
  author, "Yongbiao Xiao"; arXiv and the PDF say Xiao Yan Xu.) PDF: `../papers/arXiv_2404.09989.pdf`.
- [SX-D] F. Sun, X. Y. Xu, "Delay update in determinant quantum Monte Carlo", PRB 109, 235140 (2024),
  doi:10.1103/PhysRevB.109.235140, arXiv:2308.12005. The companion "Xu" paper (delay update; it also
  covers extended interactions and T=0). PDF: `../papers/arXiv_2308.12005.pdf`.
- [G11] E. Gull, P. Staar, S. Fuchs, P. Nukala, M. S. Summers, T. Pruschke, T. C. Schulthess, T. Maier,
  "Submatrix updates for the continuous-time auxiliary-field algorithm", PRB 83, 075122 (2011),
  doi:10.1103/PhysRevB.83.075122, arXiv:1010.3690. PDF: `../papers/arXiv_1010.3690.pdf`.
- [N09] P. K. V. V. Nukala, T. A. Maier, M. S. Summers, G. Alvarez, T. C. Schulthess, "Fast update
  algorithm for the quantum Monte Carlo simulation of the Hubbard model", PRB 80, 195111 (2009),
  doi:10.1103/PhysRevB.80.195111. Submatrix-type update for Hirsch–Fye. No arXiv version found;
  not downloaded (paywalled).
- [A08] G. Alvarez, M. S. Summers, D. E. Maxwell, M. Eisenbach, J. S. Meredith, J. M. Larkin, J. Levesque,
  T. A. Maier, P. R. C. Kent, E. F. D'Azevedo, T. C. Schulthess, "New algorithm to enable 400+ TFlop/s
  sustained performance in simulations of disorder effects in high-Tc superconductors", Proc. SC'08,
  doi:10.1109/SC.2008.5218119. Delayed updates (Hirsch–Fye). No arXiv version found; not downloaded.
- Related, not used: McDaniel et al., J. Chem. Phys. 147, 174107 (2017), arXiv:1708.00744 (delayed Slater-determinant updates, ab-initio QMC).

## Decisions
- (2026-10-30) Rank-1 eligibility, flush points, frame handling: see `docs/algorithm-notes.md` and the design spec.

## Tolerances
- To be set in Phase 2 (target ~1e-10 relative for the unit level).

## Open questions
- See the design spec, section "Open questions for Checkpoint 1".
