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
- The six sets are defined in `sets.sh` (sourced by this script and by `compare_trajectory.sh`).
- Runs six short sets (`NSweep=5 NBin=2`, `OMP_NUM_THREADS=1`) with identical parameters and seeds:
  `start`, `explicit` (`start` plus branch-only `sequential_scheme='STANDARD' N_delay=7`: lower-case key name, single-quoted upper-case value, non-default `N_delay`, which must be irrelevant for the standard scheme), `doped` (`ham_chem=-0.5 L1=4 L2=4 Beta=4`), `density` (`Mz=.F.`), `continuous`
  (`Continuous=.T.`), `tv` (see below). Outputs go to `runs/regress/<set>/{base,branch}`.
- `cmp`s every output file except `info*` (commit and time stamps) and the inputs `parameters`, `seeds`
  (they differ textually with branch-only keys); the `*_info` companions are compared. Exit 2 if a run
  produced no `confout_0` or `*_scal`, or if `../ALF_baseline` is not at `f81b5cbc`.
- Prints `<set> IDENTICAL|DIFFERENT`, exit 0 iff all identical. Differing files are listed on stderr.
- Sanity check of the harness: `BRANCH_EXTRA=NSweep=6 ...` must report `DIFFERENT` for every set.

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
  relative tolerance 1e-8 (plus an absolute floor 1e-12 of the file's largest value for noise-level entries). For a
  non-standard scheme the run's `info` must contain `Sequential scheme: <NAME>` (upper case), so a silent fallback to
  standard fails; for `standard` there is nothing to check.
- Outcome: no differing toggle -> pass iff `compare_runs.py` passes (exit 0/1). ROUNDOFF -> prints
  `DIVERGED-ROUNDOFF at proposal k (margin ...)`, the trajectories legitimately diverge, so only the `info` line is
  checked, exit 0. BUG -> exit 1. Exit 2: build, run or usage failure.
- Until the low-rank schemes exist only `standard standard` is meaningful; `standard delayed` currently fails the
  `info` check (non-standard values are silently ignored).
