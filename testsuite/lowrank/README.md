# testsuite/lowrank

Validation tools for the delayed/submatrix local updates (branch `feature/delayed-submatrix`).

## make_params.py

`make_params.py [--template DIR] OUT_DIR KEY=VALUE ...` copies `Scripts_and_Parameters_files/Start/{parameters,seeds}`
(or `DIR/{parameters,seeds}`) to `OUT_DIR` and sets each `KEY = VALUE` in every namelist that contains `KEY` (case-insensitive; the
line's comment is kept). A key present in no namelist (e.g. `Sequential_scheme`) is appended to
`&VAR_QMC`. `SEED=<int>` replaces the first seed in `seeds`. Run it with `python3 -I`.

## regress_vs_baseline.sh

Bit-identity check of the branch's default path against the unmodified baseline commit `f81b5cbc`.

```bash
testsuite/lowrank/regress_vs_baseline.sh [KEY=VALUE ...]
BRANCH_EXTRA="Sequential_scheme=0" testsuite/lowrank/regress_vs_baseline.sh   # branch-only keys
```

- Creates the worktree `../ALF_baseline` (detached at `f81b5cbc`) if missing and builds it; rebuilds the branch.
- Each side starts from its own tree's `Start` template (`--template`), since the branch template carries keys the baseline cannot parse.
- The nine sets are defined in `sets.sh` (sourced by this script and by `compare_trajectory.sh`).
- Runs nine short sets (`NSweep=5 NBin=2`, `OMP_NUM_THREADS=1`) with identical parameters and seeds:
  `start`, `explicit` (`start` plus branch-only `sequential_scheme='STANDARD' N_delay=7`: lower-case key name, single-quoted upper-case value, non-default `N_delay`, which must be irrelevant for the standard scheme), `doped` (`ham_chem=-0.5 L1=4 L2=4 Beta=4`), `density` (`Mz=.F.`), `continuous`
  (`Continuous=.T.`), `tv` (see below), `kondo`, `pv_adiab`, `z2gt` (see below). Outputs go to `runs/regress/<set>/{base,branch}`.
- `cmp`s every output file except `info*` (commit and time stamps) and the inputs `parameters`, `seeds`
  (they differ textually with branch-only keys); the `*_info` companions are compared. Exit 2 if a run
  produced no `confout_0` or `*_scal`, or if `../ALF_baseline` is not at `f81b5cbc`.
- Prints `<set> IDENTICAL|DIFFERENT`, exit 0 iff all identical. Differing files are listed on stderr.
- Sanity check of the harness: `BRANCH_EXTRA=NSweep=6 ...` must report `DIFFERENT` for every set.

Model coverage sets (all run on the baseline binary in seconds):
- `kondo`: Kondo bilayer 4x4, Beta=2; eligible and ineligible vertices (16 of 32 are not diagonal rank-1) mix in one time slice.
- `pv_adiab`: Hubbard_Plain_Vanilla with `Projector=.T. Adiabatic=.T. Theta=2 Beta=1`: time-dependent couplings `g_t`, and
  particle-hole symmetry gives `N_FL_eff = 1 < N_FL = 2` (flavour reconstruction).
- `z2gt`: Z2_Matter 4x4, `ham_U=1`, `Global_tau_moves=.T. N_Global_tau=4 Symm=.F. Checkerboard=.F.`: Global_tau moves with
  low-rank updates pending (run the trajectory check with `N_delay=5`).

The `tv` set does not use the Start default `ham_V=4`: it aborts at baseline (positive V gives an
imaginary coupling, the run is flagged unstable). It uses `ham_V=-1`, `N_SUN=1`, `L1=L2=4`, `Beta=4`
(as in `testsuite/test_branch_parameters.json`).

## sets.sh

Sourced file defining `SETS=("name|common keys|branch-only keys" ...)`, shared by `regress_vs_baseline.sh` and
`compare_trajectory.sh`.

## compare_trajectory.sh, trace_diff.py, compare_runs.py

```bash
testsuite/lowrank/compare_trajectory.sh [--rebuild] SET SCHEME_A SCHEME_B [KEY=VALUE ...]
python3 -I testsuite/lowrank/trace_diff.py --selftest
```

Two runs with identical seeds and parameters (`NSweep=3 NBin=2`) that differ only in `Sequential_scheme` must make the
same accept/reject decisions. `SET` is one of the `sets.sh` sets (its common keys only).

- `Prog/Trace_updates_mod.F90` is compiled in every build; `Upgrade2` calls `Trace_update_write` (FINAL mode, after the
  decision) only under `#ifdef ALF_TRACE_UPDATES`, and stores the `ranf_wrap()` draw in a local there (the `#else`
  branch is the original line). The trace line is `nt n_op Re(R) Im(R) Weight ranf toggle`, appended to
  `trace_updates.txt` (not MPI-safe: one rank only).
- The script builds the trace binary (`ALF_FLAGS_EXT=-DALF_TRACE_UPDATES`, `make -C Prog clean`, `make program`), copies it
  to `builds/ALF_trace.out` (git-excluded) and rebuilds the normal binary, so `Prog/ALF.out` stays the normal one.
  It is reused while it is newer than every `*.F90`, `*.h`, `Makefile` under `Prog/` and `Libraries/` (generated files
  `git.h`, `git_status.h`, `Hamiltonians_case.h`, `Hamiltonians_interface.h`, `*_read_write_parameters.F90` are ignored);
  `--rebuild` forces a rebuild. Both schemes use that one binary. `ALF_TRACE_BIN=path` substitutes a binary (testing hook).
- Runs live in `runs/traj/<SET>/<scheme>` (`<scheme>_2` for the second of two identical schemes).
- `trace_diff.py` walks both traces in step (identical RNG order: line i is the same proposal) and prints the max relative ratio
  difference over the common prefix up to the first differing toggle. At that toggle it classifies by the margin
  `|Weight-ranf|/Weight` of both runs: both <= 1e-10 is `ROUNDOFF`, otherwise `BUG`. Out-of-step `(nt, n_op, ranf)` and
  different trace lengths are `BUG`. Last line `RESULT: IDENTICAL|ROUNDOFF k margin|BUG ...`; exit 0 for IDENTICAL and
  ROUNDOFF, 1 for BUG, 2 for unusable input. `--selftest` checks the classification on synthetic traces.
- `compare_runs.py` requires `confout_0` identical field by field and all bin files (`*_scal`, `*_eq`, `*_tau`) within
  relative tolerance 1e-8 plus an absolute floor of 1e-12 times the file's largest value (the floor overrides rtol for small entries, e.g. imaginary parts of order 1e-16; the summary reports the worst ratio |difference| / allowed tolerance over all entries, <= 1 passes; `stdout.log` and `trace_updates.txt` are not compared here, the trace is compared by `trace_diff.py`). For a
  non-standard scheme the run's `info` must contain a line starting with `Sequential scheme: <NAME>` (upper case) followed by
  `,` or the end of the line (`Lowrank_init` writes `Sequential scheme: DELAYED, N_delay = 7`), so a silent fallback to
  standard fails; for `standard` there is nothing to check.
- Outcome: no differing toggle -> pass iff `compare_runs.py` passes (exit 0/1). ROUNDOFF -> prints
  `DIVERGED-ROUNDOFF at proposal k (margin ...)`, the trajectories legitimately diverge, so only the `info` line is
  checked, exit 0. BUG -> exit 1. Exit 2: build, run or usage failure.
- `standard delayed` is the trajectory check of the delayed scheme (wired into `WRAPGRUP`/`WRAPGRDO`, Task 8);
  `standard submatrix` is the trajectory check of the submatrix scheme (kernel in Task 11, same routing).

## Debug build and self-tests

`ALF_FLAGS_EXT="-DALF_LOWRANK_DEBUG"` (then `make -C Prog clean; make program`) keeps a shadow G with standard rank-1 updates and
checks `Lowrank_current_G` against it after every accept, hook and flush (tolerance 1e-10); the counts and largest errors go to
`lowrank_debug_summary.txt` in the run directory. Adding `-DLOWRANK_DEBUG_SELFTEST` corrupts the first wrap hook (delayed: skips one
X-row scaling; submatrix: skips one `Ginv` row scaling) and the run must abort with `Lowrank_debug_check: stale G ...`;
`-DLOWRANK_DEBUG_SELFTEST_PROPOSE` perturbs a proposal. In Hubbard models the submatrix hook touches pending sites only in the
down-sweep (the un-wrap acts on the site just added), so the self-test fires there. Rebuild the normal binary afterwards.

## stat_runs.sh, stat_compare.py (statistical comparison of independent chains)

```bash
testsuite/lowrank/stat_runs.sh SET SCHEME SEED [KEY=VALUE ...]     # SET = half | doped | density
python3 -I testsuite/lowrank/stat_compare.py DIR_A DIR_B           # A = reference (standard), B = candidate
python3 -I testsuite/lowrank/stat_compare.py --selftest
```

- `stat_runs.sh` builds `runs/stat/<SET>/<SCHEME>_s<SEED>[$STAT_TAG]` from `STAT_SETS` in `sets.sh` plus
  `NSweep=50 NBin=40 Ltau=0 n_skip=5 N_rebin=1 Sequential_scheme="<SCHEME>" SEED=<SEED>` (extra `KEY=VALUE` args are
  applied last and override), runs `Prog/ALF.out` (`OMP_NUM_THREADS=1`) and `Analysis/ana.out` on `Ener_scal Pot_scal
  Kin_scal Part_scal SpinZ_eq Green_eq`. One job takes 25-55 s on the M1 Pro (6x6, 2000 sweeps), so no reduction of
  NSweep was needed. Start it with `run_in_background` / `&`; keep <= 8 jobs concurrent.
- Sets: `half` 6x6 U=4 beta=5 mu=0 (Mz); `doped` 6x6 U=4 beta=5 mu=-1.0 (Mz, average sign ~0.59, NBin=80 n_skip=10, ~105 s); `density` `Mz=.F.` 6x6 U=4 beta=4 mu=0.
- ana.out output layout: `*_scalJ` rows `idx value error` (the last row is the average sign/phase; the observable
  values are jackknife sign-reweighted); `*_eqJK`/`*_eqJR`: a coordinate line, then `orb orb Re errRe Im errIm`.
  `Pot_scal` (= U sum <(n_up-1/2)(n_dn-1/2)>-like) is used as the double-occupancy proxy; `Part_scal` is
  the particle number (constant at half filling and in the SU(N) density set; handled as constant, z = 0).
- `stat_compare.py [--no-ph] DIR_A DIR_B`: z = (a-b)/sqrt(sa^2+sb^2) for every scalar (Ener, Kin, Pot, Part, avg
  sign/phase) and every Re and Im entry of `SpinZ_eqJK` and `Green_eqJK`.
  **Pass rule**: all scalar |z| < 3; for each k-point array (`*_eqJK`), on the symmetry-UNIQUE entries only (n of them):
  fraction(|z|>2) <= 0.0455 + 3*sqrt(0.0455*0.9545/n) and max|z| < 4.5. Rationale: for independent Gaussian entries
  P(|z|>2) = 4.55 %, so the number of exceedances is binomial(n, 0.0455) and the limit is its mean + 3 sigma (n = 36
  gives 15 %, n = 20 gives 18.5 %). Lattice-symmetry copies of one value are not independent samples (one
  fluctuation would count up to 4 times), so they are merged. The rms z over the unique entries is printed.
  `*_eqJR` (Fourier transform of the same data) is printed for information only.
  **Uniqueness**: two entries are merged iff all four numbers (value and error in run A, value and error in run B) agree
  to 1e-12 relative (absolute floor 1e-14), or agree after flipping the sign of both values together (Im(k) = -Im(-k)).
  Both runs enter, so a genuine difference between A and B is never merged with an entry that differs. Entries that are
  numerically zero in both runs (|value|, error < 1e-10, e.g. Im of SpinZ) are skipped (column `skip0`).
  **Non-finite input** (NaN/inf in a value or error of either run) is a FAIL: scalar and sign rows get z = nan; an
  array entry with any non-finite number is counted in column `nonfin` (checked before the zero skip) and fails
  its array (also `*_eqJR`); non-finite precision values are flagged.
  **Particle-hole merge (default, `--no-ph` switches it off)**: entries q and r are also merged if q = c - r in BOTH runs
  with the same constant c and equal errors (rtol 1e-9). For Green_eq at half filling G(k)+G(k+Q) = 2 holds to
  ~1e-11 in the output, so such partners are exact duplicates; a genuine difference in either run breaks the
  match, so nothing real can be hidden. n for Green_eqJK drops from 36 to 26 (half) and from 54 to 36 (density); no
  effect for the doped set.
  It also reports `Precision Green Mean/Max` and `Precision Phase` of both `info` files and flags B/A > 10 for the maxima.
  Exit 0 PASS, 1 FAIL, 2 unusable input.
- Calibration (4 independent standard chains per set, seeds 1-4, 6 pairs per set, 18 comparisons): with `--no-ph`
  15/18 pass (3 false fails: half 1-3, density 2-3 and 3-4, all `Green_eqJK`); with the default merge 16/18 pass (2 false
  fails, density 2-3 and 3-4). No scalar failure and no max|z| >= 4.5 in any comparison. Pooled over unique entries the
  errors are calibrated: fraction(|z|>2) 1.7-7.4 %, rms z 0.91-1.09. The density `Green_eqJK` is mildly above
  calibration (pooled fraction 7.4 %, rms z 1.09); a rebinning check (N_rebin 1-8) showed no systematic error growth.

### Validation protocol for candidate schemes (pre-registered)

For `delayed`/`submatrix` (each set: candidate seed S, e.g. 1, against the standard reference chain):
1. Scalars decide: any scalar |z| >= 3 (energy, kinetic, Pot_scal double-occupancy proxy, Part_scal, average
   sign/phase) is a bug and triggers debugging; no replicate.
2. If only an array rule (`Green_eqJK`/`SpinZ_eqJK`) fails: run ONE replicate with fresh seeds, candidate seed+10
   against standard seed 2 (`stat_runs.sh SET SCHEME <seed+10>` vs the existing `standard_s2`). The scheme passes iff the
   replicate passes everything (scalars and arrays). An array failure that replicates is a bug and triggers debugging.
3. The expected false-fail rate of the array rule for a correct scheme is about 2/18 (11 %) per comparison (measured
   above), so a single array failure alone is not evidence of a bug; the replicate is.
