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
Bit-identity regression of the branch's default path vs baseline `f81b5cbc` (6 short sets, builds
`../ALF_baseline` worktree on first use, prints `<set> IDENTICAL|DIFFERENT`, exit 0 iff all identical;
details in `testsuite/lowrank/README.md`):
```bash
cd ALF && testsuite/lowrank/regress_vs_baseline.sh      # BRANCH_EXTRA="K=V ..." adds branch-only keys
```
Trajectory comparison of two `Sequential_scheme` values (same seeds, same Metropolis decisions required;
builds/reuses the trace binary `builds/ALF_trace.out`, `-DALF_TRACE_UPDATES`; exit 0 = pass, 1 = BUG/mismatch;
a ROUNDOFF divergence prints `DIVERGED-ROUNDOFF at proposal k (margin ...)` and exits 0; `--rebuild` forces a
rebuild; details in `testsuite/lowrank/README.md`):
```bash
cd ALF && testsuite/lowrank/compare_trajectory.sh start standard standard [N_delay=8 SEED=42 ...]
```
After adding source files under `Prog/`, rerun `cmake ..` in `testsuite/tests` (test 36 globs `Prog/*.o`).
CTest unit test 36 links whatever `Prog/*.o` the last `make program` produced: a trace build (`ALF_FLAGS_EXT=-DALF_TRACE_UPDATES`) or any
other flagged build leaves incompatible objects. Always rebuild the normal binary (no `ALF_FLAGS_EXT`) and re-run `cmake ..` before `ctest`.
ALF test suite (CTest; 39 tests at baseline, 40 with `36-lowrank-unit`, all pass):
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
- (2026-10-30, Checkpoint 1, approved as proposed) Eligible = diagonal rank-1 vertices in every flavour; other
  vertices → flush + standard Upgrade2. Parameters `Sequential_scheme` / `N_delay` in `&VAR_QMC`. A trace write
  in `upgrade_mod.F90` only under `#ifdef ALF_TRACE_UPDATES`. Trajectory tolerance for bins: 1e-8 relative.
  Submatrix flushes before re-proposing a site already in the block.
- Plan: `docs/superpowers/plans/2026-10-30-delayed-submatrix.md` (execution: subagent-driven-development).
- (2026-10-30, Task 5, pre-registered validation protocol for delayed/submatrix) The scalar rule |z| < 3 is decisive: any
  scalar failure is a bug. If a candidate fails ONLY an array rule, run ONE replicate with fresh seeds (candidate seed+10
  vs standard seed 2); the scheme passes iff the replicate passes everything. An array failure that replicates is a bug.
  Measured standard-vs-standard false-fail rate of the array rule (4 seeds, 18 comparisons): 3/18 literal, 2/18 with the
  default particle-hole merge; density `Green_eqJK` is mildly above calibration (pooled fraction 7.4 %, rms z 1.09).

## Statistical harness (Task 5)
- Sets in `testsuite/lowrank/sets.sh` (`STAT_SETS`); `stat_runs.sh`/`stat_compare.py`, see `testsuite/lowrank/README.md`.
  Reference chains: `runs/stat/<set>/standard_s1` (seed 1), compared with `standard_s2` for the harness validation.
- Doped set: **mu = -1.0, beta = 5** (6x6, U=4, Mz, `NBin=80 n_skip=10`, ~105 s per job); average sign 0.58(2) / 0.59(2)
  (seeds 1/2; the requested window is [0.3, 0.8]). Pilot (mu=-1.0, NSweep=50 NBin=20, 26-36 s, seeds 1/2):
  beta=4: 0.926(11); beta=5: 0.577(21) / 0.592(17); beta=6: 0.224(19) / 0.210(24); beta=7: 0.083(17) / 0.088(21).
  So beta=5 is the smallest beta in {5,6,7} in the window. At beta=4 the sign is ~0.92 for every mu: seed 1,
  NSweep=50 NBin=20: mu=-0.6: 0.941(9), -0.8: 0.920(10), -1.0: 0.926(11), -1.2: 0.927(11), -1.4: 0.954(4),
  -1.6: 0.970(5), -1.8: 0.978(4); short pilots (NSweep=20 NBin=5): -0.3: 0.978(9), -0.5: 0.945(22), -0.8: 0.975(11).
- Statistical array rule (Controller rulings, Task 5 fix rounds 1-2): binomial 3-sigma bound on the fraction of |z|>2 over
  symmetry-unique entries (copies, sign flips and, by default, particle-hole partners G(k)+G(k+Q)=2 merged only when
  they agree in BOTH runs; `--no-ph` disables); see the README.

## Tolerances
Unit level (`testsuite/Prog.tests/36-lowrank-unit.F90`, ctest `36-lowrank-unit`; measured 2026-10-30):
- `TOL_UNIT = 1e-10` (relative ratio error and max|G_cand − G_ref|/max|G_ref|) for the cases mz, density,
  continuous, repeat (cond(1+A) ≲ 1e7). `TOL_ILL = 1e-8` for `mz_illcond` (cond(1+A) ~ 1.9e8).
- Measured drift of the standard path (real `Upgrade2`) against an exact recomputation `(1+A)^{-1}`
  (A tracked explicitly), max over 3 slice ends; cond = cond_2(1+A) at the last slice; ε = 2.2e-16:

  | case | dir | cond(1+A) | std drift | drift/(cond·ε) |
  |---|---|---|---|---|
  | mz | up / down | 4.6e6 / 6.3e6 | 1.6e-13 / 5.5e-14 | 1.6e-4 / 3.9e-5 |
  | density | up / down | 2.7e2 / 2.0e2 | 1.4e-14 / 4.2e-15 | 0.22 / 0.09 |
  | continuous | up / down | 1.5e3 / 1.2e3 | 7.2e-15 / 5.2e-15 | 0.021 / 0.019 |
  | repeat | up / down | 2.6e6 / 7.7e6 | 6.4e-14 / 1.5e-12 | 1.1e-4 / 9.0e-4 |
  | mz_illcond | up / down | 1.5e7 / 1.9e8 | 3.5e-12 / 3.5e-11 | 1.1e-3 / 8.5e-4 |

- Justification: candidate and reference each lie within about the drift of exact, so |cand − ref| ≲ 2× drift.
  TOL_UNIT = 1e-10 vs max well-conditioned drift 1.5e-12 → 65× margin; TOL_ILL = 1e-8 vs 3.5e-11 →
  ~280× margin. The margin covers N_delay-dependent summation order in the low-rank schemes.
  The drift grows with cond(1+A) (mz/repeat reach ~1e6–1e7 because diagonal wraps and flips make A
  strongly non-normal; mz_illcond uses an ill-conditioned inter-slice similarity), but it stays
  well below cond·ε (ratio 4e-5 … 0.22), i.e. at the double-precision floor.
- Asserted in the test: errors ≤ TOL (per case); drift ≤ TOL_UNIT/2 (not for mz_illcond, where the
  double-precision "exact" inverse is itself unreliable); drift ≤ 10·cond(1+A)·ε for every case. The
  last bound is loose for the high-cond cases (actual ratio ≤ 1e-3 there); it only catches gross loss.
- Inter-slice similarity: S = I + 0.1·R/√Ndim (cond(S) ≈ 1.3); mz_illcond: S = I + 0.1·R (cond(S) ≈ 20).

## Phase 3 results: delayed (verified 2026-10-30, commit 3de8a5f9)
- Unit (ctest 36): delayed vs real Upgrade2, all cases × up/down × N_delay {1,7,32,1000}: max ratio err 2.0e-12,
  max G err 1.4e-12 (mz_illcond); ≤ 8.6e-13 / 7.4e-14 for the well-conditioned cases. Kernel equivalence
  (Upgrade2 Final vs Upgrade_lowrank, 200 proposals): 0 toggle mismatches, |ΔPhase| ≤ 1e-15, identical next ranf.
- Trajectory (trace build, NSweep=3 NBin=2): IDENTICAL decisions vs standard for start, doped, density, continuous,
  kondo (mixed eligible/ineligible), pv_adiab (projector + adiabatic g_t + reconstruction), z2gt (Global_tau, N_delay=5),
  start N_delay=3.
- Statistical (delayed s3 vs standard s1): half, doped (sign 0.603(11) vs 0.590(10)), density: PASS; all scalar |z| < 3.
  Precision Green max ratio delayed/standard: 1.11, 0.25, 0.61 (stability no worse).
- Regression: 9 sets bit-identical to baseline f81b5cbc. Debug build (-DALF_LOWRANK_DEBUG): shadow-G max rel. err
  2e-14 over ~1e5 checks; self-tests fire.
- Implementation notes: lazy init at the first WRAPGRUP/WRAPGRDO (main.F90 calls Wrapgr_alloc only for Global_tau);
  G is flushed at every slice end, so the effective delay is min(N_delay, accepts per slice); submatrix rejected at
  init until Phase 4.

## Phase 4 results: submatrix (verified 2026-10-31, commit ab81f3fd)
- Unit (ctest 36), all cases × up/down × N_delay {1,7,32,1000}, max ratio / G err:
  - submatrix ≤ 3.7e-13 / 7.5e-14 in the well-conditioned cases (mz, density, continuous, repeat, mixed, asym, repeat_skip);
  - mz_illcond 2.6e-10 / 8.8e-11 (delayed 2.0e-12 / 1.4e-12), within TOL_ILL. This is the forward error of the explicit
    Γ̃⁻¹ (cond(Γ̃) up to ~1e6 in that case), not an ordering bug (two independent reviews).
  - Kernel equivalence: 0 toggle mismatches, identical next ranf. Mutations (no Ginv rescale, no repeated-site flush)
    fail 15–30 rows.
- Trajectory: IDENTICAL decisions standard vs submatrix for start, doped, density, continuous, kondo, pv_adiab, z2gt
  (N_delay=5), tv (fallback, byte-identical, message once), start with N_delay=3, Propose_S0=.T. and N_delay=100000
  (peak RSS 8 MiB). delayed vs submatrix also IDENTICAL (start, doped). The max relative ratio difference of ~5e-8 occurs
  only at |R| ~ 1e-5 (cancellation in d), the same for delayed; median 3e-15.
- Statistical (submatrix s4 vs standard s1): half, doped (sign 0.601(10) vs 0.590(10)), density: PASS, all scalar |z| ≤ 2.
  Same seed gives the same trajectory over the whole stat run: confout_0 of submatrix_s4 == standard_s4 and
  delayed_s3 == standard_s3 for all three sets.
- Stability, paired (same trajectory) Precision Green max sub/std: half 0.99, doped 1.03, density 1.92
  (2.2e-11 vs 1.1e-11); Phase max density 1.9; means 1.03–1.12. Delayed paired: 1.00 / 0.94 / 0.89.
- Regression: 9 sets bit-identical to baseline; ctest 40/40. Debug shadow max 3e-14; submatrix self-test fires
  (down-sweep only for Hubbard: the up-sweep hook never touches a pending site).
- Implementation notes: capacity capped at min(N_delay, vertices per slice[, Ndim for submatrix]); flush before
  proposing at a site already pending (any flavour); O(1) slot-invariant check in accept.

## Phase 5 results: benchmarks (verified 2026-10-31, benchmark binary 26b354c7, M1 Pro, 1 thread, Mz Hubbard Beta=2)
- Speedup vs standard at N_delay=32 (`docs/benchmarks/size_scan.md`): L=12 1.5x, L=18 2.0/2.1x, L=20 2.2/2.3x (delayed/submatrix);
  L=16 (N=256) is an outlier for the low-rank schemes (speedup dips to 1.5/1.6x; L=18 is not slower than L=16).
- N_delay scan (`docs/benchmarks/ndelay_scan.md`, generated by `bench_report.py --ndelay`): best at N_delay 16 (L=12) and 32 (L=20) for both schemes;
  within the 2.3 % noise of the best: 16-32 at L=12, 32 (delayed) / 32-64 (submatrix) at L=20; 64 is +1 to +11 % slower than 32; N_delay=256 is +33 % /
  +44 % slower than 32 at L=20, N_delay=4 +48 % / +60 %. **Default N_delay = 32** (unchanged), chosen as the single value for N = 144-400, untested
  outside. Green precision max within 0.5 % of standard for all runs; the submatrix Green mean rises mildly with N_delay (+2.6 % at L=12, +4.2 % at
  L=20 from 4 to 256); at 32 the mean differs from standard by +1.2 % (L=12) / 0.0 % (L=20) for submatrix, +0.1 % / -0.7 % for delayed.
- Profile at L=20 (`docs/benchmarks/profile.md`): standard 63 % rank-1 ZGERU, 27 % stabilisation; submatrix 60 % stabilisation, 19 % flush ZGEMM, 18 % wraps.
  THREADS=8 gives no gain. The ~24 s setup at L=20 is `validate_ham_variables` -> `Op_V_is_symmetric` (gfortran matmul).

## Open questions
- See the design spec, section "Open questions for Checkpoint 1".
