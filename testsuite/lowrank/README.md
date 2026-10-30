# testsuite/lowrank

Validation tools for the delayed/submatrix local updates (branch `feature/delayed-submatrix`).

## make_params.py

`make_params.py OUT_DIR KEY=VALUE ...` copies `Scripts_and_Parameters_files/Start/{parameters,seeds}`
to `OUT_DIR` and sets each `KEY = VALUE` in every namelist that contains `KEY` (case-insensitive; the
line's comment is kept). A key present in no namelist (e.g. `Sequential_scheme`) is appended to
`&VAR_QMC`. `SEED=<int>` replaces the first seed in `seeds`. Run it with `python3 -I`.

## regress_vs_baseline.sh

Bit-identity check of the branch's default path against the unmodified baseline commit `f81b5cbc`.

```bash
testsuite/lowrank/regress_vs_baseline.sh [KEY=VALUE ...]
BRANCH_EXTRA="Sequential_scheme=0" testsuite/lowrank/regress_vs_baseline.sh   # branch-only keys
```

- Creates the worktree `../ALF_baseline` (detached at `f81b5cbc`) if missing and builds it; rebuilds the branch.
- Runs five short sets (`NSweep=5 NBin=2`, `OMP_NUM_THREADS=1`) with identical parameters and seeds:
  `start`, `doped` (`ham_chem=-0.5 L1=4 L2=4 Beta=4`), `density` (`Mz=.F.`), `continuous`
  (`Continuous=.T.`), `tv` (see below). Outputs go to `runs/regress/<set>/{base,branch}`.
- `cmp`s every output file except `info*` (commit and time stamps) and the inputs `parameters`, `seeds`
  (they differ textually with branch-only keys); the `*_info` companions are compared. Exit 2 if a run
  produced no `confout_0` or `*_scal`, or if `../ALF_baseline` is not at `f81b5cbc`.
- Prints `<set> IDENTICAL|DIFFERENT`, exit 0 iff all identical. Differing files are listed on stderr.
- Sanity check of the harness: `BRANCH_EXTRA=NSweep=6 ...` must report `DIFFERENT` for every set.

The `tv` set does not use the Start default `ham_V=4`: it aborts at baseline (positive V gives an
imaginary coupling, the run is flagged unstable). It uses `ham_V=-1`, `N_SUN=1`, `L1=L2=4`, `Beta=4`
(as in `testsuite/test_branch_parameters.json`).
