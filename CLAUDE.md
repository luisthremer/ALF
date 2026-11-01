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
ALF test suite (CTest; 39 tests at baseline, 40 with `36-lowrank-unit`; 43 since Phase 6 with `37-Op-V-symm-fast` and the
negative entries `36-lowrank-unit-pending` / `-capacity` (Lowrank_set_N_delay must terminate; matched on its error message); all pass):
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

## Final report
- `docs/benchmarks/report.md`: method, validation summary (all five levels, with numbers), benchmarks, recommendation,
  caveats, reproduction commands, commit list.
- Final verification (2026-10-31, code commit 70176feb): ctest 40/40 in GNU noMPI, MPI and TEMPERING builds; regress 9/9 IDENTICAL;
  6 trajectory comparisons IDENTICAL (report §3.6).

## Open questions (project 1)
- See the design spec, section "Open questions for Checkpoint 1".

# Project 2: full implementation (started 2026-10-31; prompt `../PROMPTS/02-manager-prompt.md`)

Phases 0-6 run locally and autonomously; the first stop is the cluster request after Phase 6, or on **2026-12-01**
at the latest (deadline guard: SuperMUC-NG Phase 1 ends ~31.12.2026). No cluster access before explicit permission.
Design/plan approvals are replaced by independent reviewer subagents; every decision is recorded below with its reason.

## Decisions (project 2)
- (2026-10-31, user) Features first, locally; clusters afterwards as a separately approved phase.
- (2026-10-31, user) Scope: automatic N_delay; rank-N vertices (diagonal and non-diagonal) for both schemes; real-run
  validation of MPI, tempering, PARALLEL_PARAMS locally; HDF5 validated on the clusters; projector validation;
  cross-slice delay as a research prototype.
- (2026-10-31, user) Large-N setup: fast path in baseline `Op_V_is_symmetric` (all `Op_V` diagonal -> true; otherwise
  only operators with overlapping supports). Identical result; regression tier byte-identical.
- (2026-10-31, user) Cluster jobs: built and submitted by me via ssh after permission, inside the workspace dirs only;
  each benchmark campaign needs approval with a core-hour estimate; builds/validation/smoke on test/short partitions
  (<= 1 node, <= 30 min) need no further approval.
- (2026-10-31, user) No local HDF5 build (it would download and compile HDF5).
- Tolerances of project 1 stay (unit 1e-10 / ill-conditioned 1e-8, pre-registered statistical rules). New tier 6
  (parallel): observables independent of MPI rank and thread count within errors; tempering/PP runs agree with
  independent single runs.

## Phase 0: starting state (verified 2026-10-31, HEAD de4f2a4c, code 70176feb)
- `make all` ok; ctest 100 % of 40; `regress_vs_baseline.sh` 9/9 IDENTICAL; trajectories start standard vs
  submatrix and standard vs delayed: IDENTICAL decisions over all 21600 proposals.
- Audit blocking items re-checked locally: B2 (`Op_V_is_symmetric`, `Operator_mod.F90:1025`, two dense N x N matmul per
  vertex and flavour; called only if `Symm .and. .not. override_symmetric_decomposition_check`; on failure it only
  warns) confirmed by reading; I4 (scheme validated lazily in `Lowrank_init` at the first WRAPGRUP/WRAPGRDO) confirmed.
  B3/B4 (Intel toolchain, cluster stacks) cannot be checked locally: `ifx`, `ifort`, `mpiifort` are not installed.

## Phase 1: large-N enablement (verified 2026-10-31, commits 0cca1d3a, 482bda03)
- `Op_V_is_symmetric` fast path (`Prog/Operator_mod.F90`): all `O` diagonal (read from `O`, not `Op%diag`) -> `.true.`
  without products; otherwise block-local column updates `M(:,Q) <- M(:,Q)(I + g O_Q)` (repeated `P` indices accumulate
  as in the dense embedding), two Ndim^2 heap matrices instead of four, same tolerance/scale.
  Test 37 (`37-Op-V-symm-fast`): 270 random sets vs a verbatim copy of the dense version (0 mismatches); a last-write-wins
  mutation of the accumulation fails 30 cases; Ndim = 1600 diagonal < 1 ms, Ndim = 600 palindromic bonds 0.018 s
  (the dense version did not finish within 120 s: RED run).
- Setup time (Start template, Mz Hubbard, Symm=.T., time until the first WRAPUL frame, sampled; Mac M1 Pro, 1 thread):
  L=20 (N=400): baseline 24.5 s -> 0.2 s; **L=60 (N=3600): 1.9 s**. The initial stabilisation (ZGEQP3 in `WRAPUL`) is not
  included; at N=3600 it is the dominant start-up cost (> 40 s, baseline code).
- `validate_sequential_scheme` (`QMC_runtime_var_mod.F90`) runs right after the namelist read/broadcast on every rank:
  `"delayd"` -> "unknown Sequential_scheme", `Delayed`/`N_delay=0` and `N_delay=-3` -> "must be >= 1", each after 0.01 s at L=40;
  **superseded by Phase 6 (P6-1):** `N_delay = 0` (new default) is automatic (static rule), `-1` automatic with warm-up timing,
  `>= 1` fixed; only `N_delay < -1` is rejected (by `validate_sequential_scheme` and `Lowrank_init`);
  `standard` with any N_delay unchanged. Lowrank_init keeps its own check.
- regress 9/9 IDENTICAL; ctest 41/41. Symmetry warning identical to baseline for all 9 regression sets (none warns).
- Independent review (Opus subagent): no critical; important "repeated P never on the new path in test 37" fixed (case 9);
  stack temporaries fixed; NaN/Inf note added. Not fixed (cosmetic): under MPI every rank prints the validation error.

## Phase 2: design (2026-10-31; commits 04b7fa8b, d491282a and the review fold-in)
- Spec `docs/superpowers/specs/2026-10-31-rankN-design.md` (D1–D14), formulas `docs/algorithm-notes.md` Part II, plan
  `docs/superpowers/plans/2026-10-31-rankN.md` (10 tasks, subagent-driven), ledger `.superpowers/sdd/2026-10-31-rankN/`.
- Numerical prototype `docs/scripts/check_rankN.py` (7 models incl. tV/Kondo/Jz/SaSb-zz/repeated sites, up/down, capacity,
  cross-slice with complex hopping): `RESULT: PASS`, max error 9.4e-13; mutations fail (6.2, assertion, 1.2, 29).
- Independent re-derivation (Opus subagent, no access to the notes; `docs/independent-derivation-rankN.md`, its check
  `check_rankN_independent.py`: ALL CHECKS PASSED, ≤ 1.4e-14): agrees on ratio, rank-r update, L/R table, delayed hooks
  (R^T, not R^†), cross-slice; **its submatrix M-form was adopted** (D3) because my first draft (Γ̃^{-1} + Q, two accept
  paths) lost accuracy (1.1e-7 here, 9.5e-5 in the reviewer's L=12 bond test) with many repeated-site updates.
- Independent spec review (Opus subagent): approve with changes. Resolved: unconditional hook calls after both wraps (D8);
  padding moved into the post-wrap hook, no pre-wrap hook (D8); cross-slice refused with Global_tau, flush/assert at every
  G read point, LOBS constraint (D10); padding-cascade decision rule (D13); Δ = 0 columns (D14); complex hopping in the
  prototype; extra unit cases (E = 0 diagonal entries, padding past cap) and cross-slice trajectory sets.
  Not adopted: E_Exp in the hooks (private in the Operator type; roundoff only); keeping project 1's rank-1 submatrix path
  selectable (git history keeps it; tiers 2–5 are re-run on the M-form).
- Decisions without the user (recorded here per the prompt): submatrix-T (arXiv:2602.03656) is not implemented (D9); SaSb-zz
  is covered at the unit tier only, Jz at all tiers via the test model `LR_Jz_Test` (D11); N_delay now bounds the pending
  rank (columns / pending sites), unchanged for rank-1 (D5).

## Phase 3 results: rank-N delayed (verified 2026-10-31, code 5f493935 + 2fed83f6; tiers run at bf1e55a2..)
Command: `testsuite/lowrank/tiers.sh delayed 3` (summaries `runs/tiers_delayed_summary*.txt`, logs `runs/tiers/delayed/`).
- Unit (ctest 36; 15 new rank-N cases: jz, sasb_zz, tv_bond, kondo_j, ising1 (type 1), sp_type4 (type 4), nd_zero (r < N_P),
  asym_rn, gt_bonds (g_t), diag_e0, bonds_skip, bonds_illcond, rep_mixed, padcap, bigop; N_delay {1,2,7,32,1000}):
  delayed max ratio / G error 2.9e-12 / 2.3e-13 well-conditioned; bonds_illcond (Ndim 32, cond(1+A) 2.7e8 up / 6.8e7 down)
  7.8e-12 / 9.8e-13 up, 1.0e-11 / 8.4e-12 down (TOL_ILL 1e-8); kernel equivalence 0 toggle mismatches; pending-rank
  assertion CAP (removing the pre-flush fails 8 rows); mutations: no kind-2 hook → 101 rows fail, R^† instead of R^T → 4 rows.
- Trajectory: standard vs delayed IDENTICAL for all 14 sets (start, explicit, doped, density, continuous, tv, kondo,
  pv_adiab, z2gt (N_delay=5), tv_flux, tv_proj, z2, spin_peierls, jz) and tv/kondo/jz at N_delay 2 and 1000.
- Debug build (shadow G with standard rank-r updates): tv, kondo, z2, spin_peierls, jz, tv_proj, tv_flux: 3.2e5 checks,
  max relative G error 2.7e-13 (tv_proj), Mat error 1.6e-13; no termination.
- Statistical (delayed_s3 vs standard_s1, pre-registered rules): tv, kondo, jz, tv_proj, z2, doped: all PASS.
- Stability (paired, same seed => same trajectory; delayed_s3/standard_s3): Green mean 0.96–1.03, max 0.94–1.35, Phase max
  0.82–1.61 (kondo); all ≤ 2.
- Regress 13/13 IDENTICAL; ctest 41/41.
- Performance (Mac, 1 thread, Mz Hubbard, N_delay 32): the first rank-N version was 7–8 % slower for rank-1 (ZGEMM with
  n=1, Fields objects per hook); after fix round 1 (project-1 ZGEMV sequence for r = 1, no hot-path allocations) delayed is
  0.121 s/sweep (L=12) and 1.165 (L=20) vs 0.1227 / 1.178 in Phase 5 (project 1).
- Reviews: Task 1 (sonnet) approved; Tasks 2+3 (opus) Spec ✅ / approved with one Important (bonds_illcond too mild) and minors,
  all fixed in round 1 (re-review: all addressed).

## Phase 4 results: general submatrix, M-form (verified 2026-10-31, code c971a92a + fcc9f5d8 + 9b4eb424)
Command: `testsuite/lowrank/tiers.sh submatrix 4` (summary `runs/tiers_submatrix_summary2.txt`, logs `runs/tiers/submatrix/`).
- Unit (ctest 36, all rank-N cases, SUBMATRIX_RANKN = .true.): all rows pass; submatrix mz_illcond 3.6e-10 / bonds_illcond
  3.6e-11 (TOL_ILL 1e-8); mutation "no M transform in the hook" fails 142 rows (all submatrix).
- Trajectory: standard vs submatrix and delayed vs submatrix IDENTICAL for all 14 sets, plus tv/kondo/jz at N_delay 2 and 1000
  (34 comparisons).
- Debug shadow: 7 sets, ≤ 1.7e-14; padded sites are counted in the summary.
- **Bug found by the tiers (systematic debugging):** with K_ND = 8 the Kondo statistical run had Precision Green max 2.6e-9 vs
  4.0e-11 for standard and Phase max 7.7e-10 vs 2.8e-11 — on the identical trajectory (confout equal). Evidence: debug shadow
  over the full run (5.4e6 checks) max 1.2e-11 per block, no wrong update; N_delay = 2: 5.2e-11; K_ND = 4 / 8: 3.2e-9 / 2.6e-9;
  K_ND = 2 / 1: 4.9e-11. Root cause: the M-form's forward error with non-diagonal content in a block (cancellation in
  G0 − G0(:,S) M (…)), accumulating between stabilisations. Fix (9b4eb424, spec D15): K_ND = 1 — non-diagonal content is never
  kept in the M-form; unit tier passes.
- Statistical (submatrix_s4 vs standard_s1): tv, jz, tv_proj, z2, doped PASS; kondo FAIL on SpinZ_eqJK only (array rule),
  identical to standard_s4 vs standard_s1 (a property of the seed-4 chain; submatrix_s4's confout == standard_s4's); the
  pre-registered replicate (submatrix_s14 vs standard_s2) PASSES — kondo passes.
- Stability (paired): Green mean 0.99–1.11, max 0.53–1.32, Phase max 0.92–1.54 (all ≤ 2).
- Regress 13/13 IDENTICAL; ctest 41/41.
- Projector (Task 7): pv_adiab and tv_proj PASS for delayed (s3) and submatrix (s4); confout == standard with the same seed;
  pv_adiab Precision Green max 7.36e-8 (submatrix) vs 7.37e-8 (standard).
- Reviews: Task 5 (opus) Spec ✅ / approved (minors fixed).

## Phase 5: cross-slice delay prototype (Task 8; code c7071a4f + fde6d2cc, merged; decision 2026-10-31)
- Validation: unit variant (carry through the inter-slice similarity) all rows ≤ TOL; trajectories with the key on IDENTICAL for
  start (Symm), tv, kondo, tv_flux (complex hopping), pv_adiab (projector, reconstruction) × delayed/submatrix × default and
  one-slice measuring window (20/20); key on with standard: regress 13/13 IDENTICAL; refused with N_Global_tau > 0; debug runs
  without termination (shadow ≤ 1.6e-12); review found a K_ND bypass across the delayed→submatrix mode switch (fixed, with a
  debug assertion that fires on the unfixed code); kondo precision ratio with the key on 1.47.
- Measurement (`docs/benchmarks/cross_slice.md`, generated by `testsuite/lowrank/quiet_report.py xs` from
  `raw_cross_slice.txt`; Mz Hubbard Beta 2, checkerboard, U = 8, one measured slice, N_delay 32, 3 reps): key on/off = 1.074
  (delayed L=8), 0.967 (submatrix L=8), 1.016 / 1.017 (L=12) — no significant difference anywhere (rule: > 5 % and > 2 × spread).
- **Decision: keep as opt-in (default .F.), documented as "no measured gain, not recommended"; flagged for the user's review.**
  This deviates from my pre-registered rule (no gain → remove). Reason: removal would rework validated code paths (13 read
  points in main.F90, the automatic-N_delay skip, the timers) for no correctness or default-path benefit (priority order:
  correctness > default unchanged > speed > elegance); the default path is unaffected while the key is off. If you prefer
  removal, it is one reverted feature with the same validation tiers.

## D13: delayed vs submatrix for bond models (Mac, 1 thread, N_delay 32; `docs/benchmarks/d13_bond_models.md`)
Generated by `testsuite/lowrank/d13_report.py` from `docs/benchmarks/raw_d13_bond_models.txt` (bench.sh, 3 repetitions,
binary with K_ND = 1). Rule: submatrix slower if > 5 % and > 2 × the larger spread.
- tV: L=12 submatrix 1.64× delayed's time (slower); L=8 equal within noise. Kondo: 1.21× (L=8), 1.69× (L=12) (slower).
  Spin-Peierls: L=8 submatrix faster (0.87×), L=12 equal within noise.
- Speedups vs standard at L=12: delayed tV 1.59, Kondo 2.08, Spin-Peierls 1.30; submatrix 0.97 / 1.23 / 1.21.
- **Recommendation (D13): delayed for models with non-diagonal (bond) vertices** (tV, Kondo-J; Spin-Peierls either);
  submatrix for diagonal (Hubbard-type) models, where it is 2–5 % faster than delayed (project 1). Cause: K_ND = 1 (D15)
  flushes submatrix blocks at every non-diagonal update, so for bond vertices submatrix costs about as much as standard.
  Runs 30–54 (repetitions 2–3 of most points) may have overlapped intermittent 2-core builds of another task; the clear
  verdicts (tV L=12, Kondo: ≥ 21 %) are far above that noise; the Spin-Peierls verdicts are marginal and not relied on.

## Phase 5 (project 2): local parallel tier (Task 9; script `testsuite/lowrank/parallel_tier.sh`)
- First run on d374e204 (submatrix still rank-1 there): mpi 18/18 (1/2/4 ranks, half and tv, 3 schemes), threads 3/3,
  tempering 3/3, pp 3/3 PASS (some after the pre-registered replicate); scheme line once per Temp_<g>/info; max scalar |z| 2.44.
  Harness corrections after run 1 (independent seed for the second tempering reference; replicate judged per replica,
  the other replica still checked for hard failures) are recorded in the ledger; no final verdict depended on them.
- Review fixes (e6b5569d, b7ade70d): NaN = hard fail, tempering must show accepted exchanges, same-seed omp1/omp4 identity,
  single provenance (archived HEAD), concurrency budget.
- **Full tier on the merged code (HEAD bdf38fa6: rank-N delayed, M-form submatrix with K_ND = 1, cross-slice key, auto N_delay
  default 0)**, `testsuite/lowrank/parallel_tier.sh --rebuild` (summary `runs/parallel_full_summary.txt`): mpi 18/18 PASS
  (half and tv; standard/delayed/submatrix; 1, 2, 4 ranks), threads: same-seed omp1 vs omp4 bit-identical for all schemes and
  independent-seed comparison PASS (standard after its replicate), tempering 3/3 PASS (exchange acceptance 1.2–1.7 %),
  PARALLEL_PARAMS 3/3 PASS (standard after its replicate). Max scalar |z| 2.05. Wall: mpi 791 s, threads 580 s, tempering
  297 s, pp 288 s. Tempering with `N_delay=-1` (`runs/parallel_tune_summary.txt`): 3/3 PASS, statistics identical to the
  static run (N_delay does not change the chain), one `Auto N_delay` decision line per `Temp_<g>/info`, groups decide
  independently (submatrix: 16 and 8).

## Work packages (project 2)
| WP | Phase | Content | Done means |
|---|---|---|---|
| 1a | 1 | `Op_V_is_symmetric` fast path: all diagonal -> true; else block-local column updates of the forward/reverse products (O(n_op Ndim N_P^2) instead of O(n_op Ndim^3)) | unit test new == old on random sets; regress IDENTICAL; setup time at N = 3600 recorded |
| 1b | 1 | Validate `Sequential_scheme`, `N_delay` right after the parameter read | typo aborts before the Hamiltonian setup; regress IDENTICAL |
| 2 | 2 | Design + plan: rank-N_P diagonal/non-diagonal vertices (delayed, submatrix), multi-site submatrix blocks, projector, cross-slice delay | algorithm notes extended, independent re-derivation agrees, plan reviewed |
| 3 | 3 | Rank-N delayed | tiers 1-5 for tV, Kondo-J, Z2 or Spin-Peierls, Jz/SaSb-zz, type-1 with Delta != 0, sign-problem set |
| 4 | 4 | Rank-N submatrix | same tiers; decisions identical to delayed |
| 5 | 5 | Projector statistical validation; cross-slice prototype (keep/opt-in/remove decided by measurement); MPI rank-count independence; tempering and PP real runs | tier 3 for a projector set; tier 6 locally |
| 6 | 6 | Automatic N_delay (static rule from the papers, clamped; optional warm-up timing, rank 0 decides, broadcast; >5 % rule), fixed N_delay override, logged in `info` | auto value within noise of the best fixed value (Mac); chain identical except near-ties |
| 7 | after 6 | Cluster request (first stop) | message with evidence, open questions, benchmark matrix with core-hour estimates |
| 9 | any | Polish without cluster: `files.tex`, `default_parameters_generic.py`, §7 notices | docs build; generator knows the keys |

## Phase 6 results: automatic N_delay (functional part; verified 2026-10-31, code 3ed73f7e + 9ce96c27 + 6914a392)
- `N_delay = 0` (new default): static rule; `-1`: + warm-up timing (sweeps 2–3, group-rank-0 decision, > 5 % rule); `≥ 1` fixed;
  `< -1` rejected at the parameter read ("must be >= 1 (fixed), 0 (automatic) or -1 (automatic with warm-up timing)").
- Unit (ctest 43/43): static-rule table (56 cells + 10 boundary rows), set_N_delay down/up with pending rank reaching the set
  value, timing state machine and all log variants; two negative tests matched by their exact error messages.
- Trajectory: N_delay = 0 is byte-identical (traces included) to fixed n0 (start, tv × both schemes); standard vs scheme with
  N_delay 0 and −1 IDENTICAL for start, tv, kondo, jz, pv_adiab, tv_proj (both schemes); fixed 8 vs −1 IDENTICAL (kondo, jz).
- Edge cases (`runs/p63_edge`): NSweep 1–2 → "tuning incomplete", n0 kept; decision across a bin boundary; CPU_MAX truncation;
  Langevin (no low-rank init, no line); projector tV decides; N_delay = −2 rejected. (Nt_sequential_end < start is overridden by
  the Hamiltonians when Global_tau is off, so it cannot occur.) MPI 2 ranks and tempering 2 groups: see the parallel tier.
- Mac timing (`docs/benchmarks/auto_ndelay.md`, generated by `quiet_report.py nd` from `raw_auto_ndelay.txt`; binary 0f6a7651, 1 thread,
  3 reps; Mz Hubbard L = 12/16/20, tV L = 8/12, both schemes; fixed 8–128 vs 0 and −1): **all 20 verdicts ok** (within max(5 %, 2 × spread)
  of the best fixed). Static rule (0): worst +3.8 % (Hubbard L=20 delayed, n0 = 16 vs best 32), +3.7 % (tV L=8 delayed, 16 vs 8).
  Timing (−1): worst +2.4 %; it picked 32 where the static rule's 16 was slower (Hubbard L=16/20 delayed) and 8 for tV L=8 delayed.
  The −1 choice was the same in all 3 repetitions of every point. One outlier repetition (tV L=12 submatrix fixed 32: 4.18 s vs 0.37)
  is noted in the table; no verdict depends on it.

## Benchmark instrumentation (plan `docs/superpowers/plans/2026-10-31-bench-instrumentation.md`, all tasks complete 2026-11-01)
- T1 timers (`-DALF_TIMERS`, `testsuite/lowrank/timers_build.sh`), T3 campaign generator (`testsuite/lowrank/campaign/`), T2 BLAS ceiling
  microbenchmark (`testsuite/lowrank/blas_bench.sh`; MPI build with a barrier per point for srun/mpiexec, memory guard, `--probe` for the
  MKL provenance check; test `campaign/tests/test_blas_bench.sh`, 143 PASS; Homebrew mpiexec needs `FI_PROVIDER=tcp` on the Mac).
  Parked minors (cluster-only paths: barrier timeout/ABORT, Linux OMP_PLACES with inherited cpusets) are in the ledger; set
  `BLAS_BENCH_READY_TIMEOUT=3600` for N = 8000 jobs.
- Cluster request (first stop): `docs/cluster/2026-10-31-cluster-request.md`.

## Cluster-phase items (prepared locally, to be checked on the clusters)
- HDF5 provenance (P6-2): non-standard schemes write `Sequential_scheme`, `N_delay_requested` (at the parameter write) and
  `N_delay_used` (end of run) as attributes of `/` in `data.h5`; standard runs write nothing new (baseline data.h5 kept).
  Code only preprocessed and compile-checked against stub hdf5/h5lt modules locally (no HDF5 here); on the cluster: build with
  the HDF5 module, run delayed and standard, check `h5dump -A -g / data.h5`. Semantics: `N_delay_requested` is the
  parameter (≥ 1 fixed, 0 auto, −1 auto+tune); `N_delay_used` is 0 if no low-rank state was initialised (Langevin, HMC-only);
  attributes of an earlier non-standard run remain if a standard run later reuses the same data.h5 (review P6-2, approved).

- Large-N precision with the automatic default (final review, Important 2): at Ndim 1600–3600 the static rule gives n0 = 128
  (delayed) / 128–256 (submatrix), never run locally beyond Ndim 400 in the statistical tiers. Phase-7 step: Mz Hubbard L = 40 and
  60, U = 4 and 8, Beta = 4, same seed: standard vs auto vs fixed 32, compare Precision Green mean/max and Phase (paired, ≤ 2×).
  Mac pre-check done (`docs/benchmarks/largeN_precision.md`, `precision_report.py`): L = 32 (U = 4, 8) and L = 40 (U = 4), Beta 2, same seed,
  auto (n0 64; submatrix 128 at L = 40) and fixed 32 vs standard: same acceptance everywhere, Green max ratio 1.000, mean 0.99–1.02,
  Phase max ≤ 1.61. At U = 8 standard itself reaches Green max 8.8e-2 with Nwrap = 10 → C0 needs a shorter Nwrap at U = 8.
- Intel toolchain: ifx/ifort + MKL builds never made (B3/B4); first build on each machine compiles all variants (noMPI, MPI,
  TEMPERING, HDF5, debug, trace, timers) and runs ctest + regress + a trajectory subset.
- Untested combinations (low risk by construction — G flushed and asserted at every WRAPGRUP/WRAPGRDO exit and between sweeps;
  one short trajectory/stat run each in Phase 7): HMC + Sequential with a low-rank scheme; flavour reconstruction (N_FL_eff < N_FL)
  with non-diagonal rank-N vertices; N_delay = -1 with Global_tau; cross-slice under MPI/tempering. Global_moves cannot be tested
  (no shipped model implements Global_move). Debug/trace builds under MPI write one summary/trace file for all ranks (diagnostic
  builds only).

## Open questions for the user (collected for the cluster request)
00. Cross-slice prototype: no measured gain on the Mac; kept opt-in instead of removed (deviation from my own pre-registered
   rule, reason in "Phase 5: cross-slice"). Keep or remove?
0. Helma + fakeintel: ALF's `HELMA` preset in configure.sh links with `-static` (INTELLLVMOPTFLAGS); `LD_PRELOAD` cannot act on a
   static binary, so the fakeintel library cannot take effect. Proposal: build the Helma benchmark binary without `-static`
   (dynamic MKL) and keep fakeintel; the campaign jobs abort if fakeintel is requested for a static binary. Which do you prefer,
   and how is fakeintel meant to be loaded on Helma?
1. Submatrix accuracy guard (D15, final value `K_ND = 1`): non-diagonal content is never kept in the M-form (found by the Kondo
   precision tier: 65× worse Precision Green max with K_ND = 8; K_ND = 1 restores standard level). Consequence (D13): for bond
   models (tV, Kondo-J) submatrix costs about as much as standard for the bond vertices; delayed is recommended there. If you
   want a fast submatrix for bond models, the option is a hybrid that routes non-diagonal blocks to the delayed form (not built).
