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
- The thirteen sets are defined in `sets.sh` (sourced by this script and by `compare_trajectory.sh`).
- Runs thirteen short sets (`NSweep=5 NBin=2`, `OMP_NUM_THREADS=1`) with identical parameters and seeds:
  `start`, `explicit` (`start` plus branch-only `sequential_scheme='STANDARD' N_delay=7`: lower-case key name, single-quoted upper-case value, non-default `N_delay`, which must be irrelevant for the standard scheme), `doped` (`ham_chem=-0.5 L1=4 L2=4 Beta=4`), `density` (`Mz=.F.`), `continuous`
  (`Continuous=.T.`), `tv` (see below), `kondo`, `pv_adiab`, `z2gt` (see below). Outputs go to `runs/regress/<set>/{base,branch}`.
- `cmp`s every output file except `info*` (commit and time stamps) and the inputs `parameters`, `seeds`
  (they differ textually with branch-only keys); the `*_info` companions are compared. Exit 2 if a run
  produced no `confout_0` or `*_scal`, or if `../ALF_baseline` is not at `f81b5cbc`.
- Prints `<set> IDENTICAL|DIFFERENT`, exit 0 iff all identical. Differing files are listed on stderr.
- Sanity check of the harness: `BRANCH_EXTRA=NSweep=6 ...` must report `DIFFERENT` for every set.

Model coverage sets (all run on the baseline binary in seconds):
- `kondo`: Kondo bilayer 4x4, Beta=2; diagonal rank-1 (U_c, U_f) and non-diagonal bond vertices (V_SUN/VJ_SUN, 16 of 32) mix in
  one time slice. Since rank-N delayed (5f493935) all are delayed-eligible, and since the submatrix M-form (Task 5) all are submatrix-eligible.
- `pv_adiab`: Hubbard_Plain_Vanilla with `Projector=.T. Adiabatic=.T. Theta=2 Beta=1`: time-dependent couplings `g_t`, and
  particle-hole symmetry gives `N_FL_eff = 1 < N_FL = 2` (flavour reconstruction).
- `z2gt`: Z2_Matter 4x4, `ham_U=1`, `Global_tau_moves=.T. N_Global_tau=4 Symm=.F. Checkerboard=.F.`: Global_tau moves with
  low-rank updates pending (run the trajectory check with `N_delay=5`).

Rank-N sets (multi-site, non-diagonal vertices; the validation targets of the rank-N delayed/submatrix work; each runs in
well under a second with the trajectory harness settings; all four also run on the baseline binary). The multi-site
vertices of each were verified to be proposed and accepted with the trace build (standard scheme, `NSweep=3 NBin=2`):
- `tv_flux`: tV 4x4, `ham_V=-1 N_SUN=1 Beta=4 N_Phi=1` (complex hopping; Symm families: 56 bond vertices, all N_P=2;
  26880 proposals, 10463 accepted).
- `tv_proj`: tV 4x4, `ham_V=-1 N_SUN=1 Beta=1 Projector=.T. Theta=2` (33600 proposals, 13010 accepted).
- `z2`: Z2_Matter 4x4, `ham_U=1 Beta=1 Symm=.F. Checkerboard=.F.` without Global_tau moves, so `Nt_sequential_start/end`
  default to `1 .. n_op` (81 vertices: 16 Hubbard, 32 gauge bonds, 32 matter bonds, 1 dummy site operator with g=0).
  Gauge bonds: 3840 proposals, 195 accepted; matter bonds: 3840 / 3601.
- `spin_peierls`: Spin_Peierls 4x4, `Beta=2 N_SUN=2 ham_Jx=ham_Jy=0.4 ham_U=1.5 ham_Lambda=0.1 ham_omega0=0.5`
  (48 vertices: 16 site, 32 type-4 bonds with `flip_protocol=2`; bonds 7680 proposals, 1507 accepted).

## sets.sh: TRAJ_ONLY_SETS and the test model LR_Jz_Test

`TRAJ_ONLY_SETS=("name|common keys|")` holds sets that only `compare_trajectory.sh` runs (it reads `SETS` and
`TRAJ_ONLY_SETS`; `regress_vs_baseline.sh` reads `SETS` only), because the baseline binary does not know their model.
Currently `jz`: the test-only Hamiltonian `LR_Jz_Test` (`Prog/Hamiltonians/Hamiltonian_LR_Jz_Test_smod.F90`, registered in
`Prog/Hamiltonians.list.d/lowrank_test`; not a production model): square lattice, `N_FL = 2`, `N_SUN = 1`, Mz-Hubbard sites
(`Predefined_Int_U_MZ`, `ham_U`) plus one `Predefined_Int_Jz` vertex (N_P = 2, `ham_Jz`) per nearest-neighbour bond,
`Symm = .F.`, one dense hopping operator, no projector. Operator order: `Ndim` site vertices, then per site the x-bond and
the y-bond. Namelist `&VAR_LR_Jz_Test`: `ham_T, ham_U, ham_Jz, ham_chem, Dtau, Beta` (plus `&VAR_lattice`). The Start
template has no such namelist: `make_params.py` appends it (`EXTRA_NAMELISTS`) when `ham_name` is set to `LR_Jz_Test`.
The `jz` set uses `L1=L2=4 ham_U=4 ham_Jz=1 Beta=2 Nwrap=5` (Nwrap=10 gives Precision Green max 1e-5, Nwrap=5 1.5e-9);
48 vertices per slice, 16 site (1583 of 3840 proposals accepted) and 32 Jz bonds (3029 of 7680).
Observables: Kin, Pot (U sum n_up n_do plus the Jz bond energy), Part, Ener; Green, SpinZ, SpinXY, SpinT, Den.

The `tv` set does not use the Start default `ham_V=4`: it aborts at baseline (positive V gives an
imaginary coupling, the run is flagged unstable). It uses `ham_V=-1`, `N_SUN=1`, `L1=L2=4`, `Beta=4`
(as in `testsuite/test_branch_parameters.json`).

## sets.sh

Sourced file defining `SETS=("name|common keys|branch-only keys" ...)`, shared by `regress_vs_baseline.sh` and
`compare_trajectory.sh`.

## compare_trajectory.sh, trace_diff.py, compare_runs.py

```bash
testsuite/lowrank/compare_trajectory.sh [--rebuild] SET SCHEME_A SCHEME_B [KEY=VALUE ...] [A:KEY=VALUE ...] [B:KEY=VALUE ...]
testsuite/lowrank/compare_trajectory.sh start delayed delayed A:N_delay=16 B:N_delay=-1   # fixed vs automatic, same scheme
python3 -I testsuite/lowrank/trace_diff.py --selftest
```

Two runs with identical seeds and parameters (`NSweep=3 NBin=2`) that differ only in `Sequential_scheme` must make the
same accept/reject decisions. `SET` is one of the `sets.sh` sets (its common keys only). Plain `KEY=VALUE` arguments apply to
both runs; `A:KEY=VALUE` / `B:KEY=VALUE` only to run A / B (after the plain ones, so they override), e.g. a fixed against an
automatic `N_delay` of the same scheme. Without `N_delay` both runs use the Start template's default `N_delay = 0` (automatic,
static rule). `MAKE_JOBS=n` sets the parallel make jobs of the trace/normal builds (default 8).

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
X-row scaling; submatrix: skips the first row transform of `M`, a row scaling or the `L` product of a non-diagonal wrap) and the run must abort with `Lowrank_debug_check: stale G ...`;
`-DLOWRANK_DEBUG_SELFTEST_PROPOSE` perturbs a proposal. The SELFTEST macros only take effect together with `ALF_LOWRANK_DEBUG`.
In Hubbard models the submatrix hook touches pending sites only in the
down-sweep (the un-wrap acts on the site just added), so the self-test fires there. A run whose hook never touches a pending entry
passes silently, which for Hubbard means the self-test fires only in the down-sweep. For submatrix the summary also lists the
padded sites (`padded sites (accepts)`: active sites added to S by accepts; `padded sites (hooks)`: sites P \ S padded by a
non-diagonal wrap that partially overlaps S).
`LOWRANK_DEBUG_SELFTEST_PROPOSE` has only been exercised for delayed. Rebuild the normal binary afterwards.

## stat_runs.sh, stat_compare.py (statistical comparison of independent chains)

```bash
testsuite/lowrank/stat_runs.sh SET SCHEME SEED [KEY=VALUE ...]     # SET = half | doped | density | tv | kondo | jz | tv_proj | z2 | pv_adiab
python3 -I testsuite/lowrank/stat_compare.py DIR_A DIR_B           # A = reference (standard), B = candidate
python3 -I testsuite/lowrank/stat_compare.py --selftest
```

- `stat_runs.sh` builds `runs/stat/<SET>/<SCHEME>_s<SEED>[$STAT_TAG]` from `STAT_SETS` in `sets.sh` plus
  `NSweep=50 NBin=40 Ltau=0 n_skip=5 N_rebin=1 N_delay=32 Sequential_scheme="<SCHEME>" SEED=<SEED>` (extra `KEY=VALUE` args are
  applied last and override; `N_delay=32`, the pre-Phase-6 default, is pinned so that the recorded results and the existing
  chains in `runs/stat` stay comparable now that the template default is `0` = automatic), runs `Prog/ALF.out` (`OMP_NUM_THREADS=1`) and `Analysis/ana.out` on `Ener_scal Pot_scal
  Kin_scal Part_scal SpinZ_eq Green_eq`. One job takes 25-55 s on the M1 Pro (6x6, 2000 sweeps), so no reduction of
  NSweep was needed. Start it with `run_in_background` / `&`; keep <= 8 jobs concurrent.
- Rank-N sets (4x4, standard runs 6-67 s, standard seeds 1/2 give `VERDICT: PASS` for all six): `tv` (N_SUN=1, `ham_V=-1`,
  `ham_chem=-1`, `Beta=2`, `NBin=80`; average sign 0.59(1); at half filling the N_SUN=1 sign is ~0.03), `kondo` (Bilayer 4x4,
  Beta=2), `jz` (`LR_Jz_Test`, `Beta=4 Nwrap=5`), `tv_proj` (projective tV with `N_SUN=2 ham_V=0.5 N_Phi=1 Theta=2 Beta=1
  NBin=80`: the N_SUN=1 projective sets have average sign ~0.05).
  `z2` (Z2_Matter, `ham_U=1 Beta=2 Symm=.F. Checkerboard=.F.`, no Global_tau, `NSweep=200`; average sign 1) and `pv_adiab`
  (Hubbard_Plain_Vanilla, `Projector=.T. Adiabatic=.T. Theta=2 Beta=1`, as in SETS; 6 s). Z2_Matter has no Ener/Pot/Kin
  observables (Part, Flux, X, Q), so `stat_runs.sh` analyses the standard observables that exist and, only when `Ener_scal`
  does not exist, every other `*_scal` file; `stat_compare.py` compares every `*_scalJ` present. With `NSweep=50` the `z2`
  errors were underestimated by the gauge-field autocorrelation (standard s1 vs s2: Green_eqJK 4 of 16 entries |z|>2, rms z 1.5,
  FAIL); with `NSweep=200` (8000 sweeps, 66 s) 5 of the 6 pairs among 4 standard seeds pass (s1/s2 PASS; s3/s4 FAIL on
  Green_eqJK alone, i.e. within the documented array-rule false-fail rate; the pre-registered replicate protocol applies).
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

## parallel_tier.sh (tier 6: MPI rank count, threads, tempering, PARALLEL_PARAMS)

```bash
testsuite/lowrank/parallel_tier.sh [--rebuild] [--no-replicate] [mpi] [threads] [tempering] [pp] [KEY=VALUE ...]
PT_SCHEMES="standard delayed" PT_SETS=half PT_NP="1 4" PT_TAG=_smoke testsuite/lowrank/parallel_tier.sh mpi NSweep=4 NBin=12   # smoke test
```

- Builds each variant once per source state (tree hashes of `Prog`, `Libraries`, `Analysis`, `Makefile`, `configure.sh` at HEAD) in a copy of
  the committed tree, `builds/<variant>/src` (`git archive HEAD`; the main tree's `Prog/` is never rebuilt): `noMPI`, `MPI`,
  `TEMPERING`, `PARALLEL_PARAMS` (`source configure.sh GNU <mode> NO-INTERACTIVE && make all -j4`, log `runs/parallel/build_<variant>.log`).
  `builds/noMPI/src/Analysis/ana.out` analyses every run. Runs and logs: `runs/parallel/<phase>/...`, independent noMPI reference
  chains in `runs/parallel/ref/` (cached, keyed by the build stamp and parameters). Summary lines: `runs/parallel/summary_<timestamp>.txt`.
- One summary line per check, `PASS|FAIL|FAIL-ARRAY <check>: scalars max|z|=.. (obs) sign z=..; Green_eqJK n=.. frac>2=..(lim ..) max|z|=.. rms=..;
  SpinZ_eqJK ...; wall=..` (`parallel_compare.py` wraps `stat_compare.py`, whose rules are unchanged; full output in `info_stat_compare.log` of the run directory).
  Exit 0 iff all cases pass.
- `mpi`: sets `half` (6x6, U=4, Beta=5: the `start` set) and `tv` (STAT_SETS), each scheme, 1, 2 and 4 ranks (MPI build), against one
  independent standard noMPI chain (seed family 1000). ALF with MPI runs one chain per rank and averages the per-rank bins
  (every bin is the mean over the ranks), so every rank count has the same number of bins (`NSweep=50 NBin=40 n_skip=5`, `tv`: `NBin=80`) and the
  N-rank run has N times the statistics. The rank r of an N-rank run reads line N-r of `seeds`, so the script writes 16 distinct seeds
  (`base + 7919 i`) per run.
- `threads`: noMPI, `half`, each scheme, two checks. (a) statistical: `OMP_NUM_THREADS=4` against `1` with independent seeds (the usual rules).
  (b) same seed, `OMP_NUM_THREADS=1` against `4`: `compare_runs.py` must find `confout_0` identical and every bin within rtol 1e-8 (the thread count
  must not change the chain; measured bit-identical bins; this is the sharper test, since (a) is then the same test as the 1-rank MPI case). A failure of (b) is a hard failure.
- `tempering` / `pp`: 2 replicas (`mpirun -n 2`, `mpi_per_parameter_set=1`) of `half` with `ham_U = 3` (Temp_0) and `4` (Temp_1), each scheme;
  every replica is compared with an independent noMPI standard run of the same parameters (`ham_U=3`, `4`). TEMPERING reads `&VAR_QMC` and
  `&VAR_TEMP` from the top-level `parameters` and the model namelists from `Temp_<g>/parameters`; PARALLEL_PARAMS reads `&VAR_QMC`/`&VAR_HAM_NAME`
  from `Temp_<g>/parameters`; the script writes the same keys into both files (`Start_Tempering` template, `N_exchange_steps=6`,
  `N_Tempering_frequency=2`). TEMPERING: `ham_U` 3 against 4 on a 6x6 lattice at Beta=5 exchanges rarely (acceptance 1-2 % per attempt;
  with the template frequency 10 only ~16 exchanges per run, hence frequency 2, ~80). PARALLEL_PARAMS has no exchange steps (`Acceptance
  Tempering` 0): its replicas are independent chains in one MPI job, and with the same seed a PP replica is bit-identical to the corresponding noMPI/MPI run.
  The `ham_U=3` and `ham_U=4` references have different seeds (reference seed, reference seed + 500).
- For every non-standard run, `info` (`Temp_<g>/info`) must contain `Sequential scheme: <NAME>` exactly once. TEMPERING additionally requires
  `Acceptance Tempering` in `Temp_0/info` to be > 0 (else the run is vacuous: hard failure); PARALLEL_PARAMS reports it (expected 0).
- Non-finite array entries (NaN/inf, `nonfin` column of `stat_compare.py`) are a hard failure, never replicated.
- Provenance: binaries, parameter templates, `sets.sh`, `make_params.py` and `compare_runs.py` all come from the archived HEAD copy `builds/noMPI/src`
  (the build stamp includes their tree hashes); `stat_compare.py`/`parallel_compare.py` come from the working tree and the script warns if
  `testsuite/lowrank` has uncommitted changes. The build subshell unsets `ALF_FLAGS_EXT`; `PT_MAKE_J` (default 4) sets `make -j`.
- Seeds: reference chains 1000-family, candidates 2000 (+100 delayed, +200 submatrix, so each scheme is its own chain), replicate
  candidates 3000 (+100/+200), replicate references 4000. A reference and a candidate never share a seed.
- Pre-registered protocol as above (a hard failure is a scalar |z| >= 3, non-finite data, a missing scheme line, a vacuous TEMPERING run, a failed
  same-seed thread check or a failed run): a scalar failure is a bug; an array-only failure (`FAIL-ARRAY`) triggers ONE replicate of the case with the
  fresh seeds (`--no-replicate` switches that off), and the case passes iff the replicate passes everything. In `tempering`/`pp` the replicate concerns the comparison (replica) that failed:
  in round 2 only that replica is judged fully (both replicas are rerun, one MPI job); the other replica is still checked (scheme line, scalars,
  non-finite data: a hard failure fails the case) but an array-only result for it is ignored (its round-1 verdict stands).
- Budget (M1 Pro, at most 3 cores busy in the `mpi` phase, 4 in `threads` for the OMP=4 run): a 6x6 run takes ~50 s whatever the rank count, a 4x4 `tv` run 23-38 s.
  In `mpi` the np1 and np2 runs of a (set, scheme) run concurrently (3 cores), np4 alone (4 ranks): per (set, scheme) ~100 s for `half`, ~60 s for `tv`,
  estimated ~9.5 min for the phase without replicates (first full run, serial: 1028 s with replicates, ~13 min without; the concurrent timing has to be confirmed by
  the next end-to-end run). First full run: `threads` 401 s (now +3 same-seed runs, ~7.5 min), `tempering` 254 s, `pp` 249 s; builds ~1.5 min per variant.
- Statistics of the tier: the single reference chain of a set is shared by all rank counts and schemes of that set (np1/np2/np4 candidates are also nested:
  the np1 chain is one of the np2/np4 chains), so the comparisons are correlated, and a fluke of the reference fails several cases at once (the replicate uses a second
  reference chain). In the first full run the `half` reference of seed family 1000 had an unusually low SpinZ at (pi,pi) (3.58 +- 0.34 against 4.2-5.0 for the four
  `runs/stat/half/standard_s*` chains and the candidates): the three standard cases and `pp/half/standard` failed an array rule only and passed their replicates.
  The array rule has ~11 % false-fail rate per comparison; with one replicate ~1.2 % per comparison, i.e. about 1 in 3 for the ~33 comparisons of the whole tier
  (less when they are correlated). A single case that fails after its replicate through an array rule alone is therefore a prompt for one more round with fresh seeds before debugging; a scalar
  failure, non-finite data or a failed same-seed thread check is not.

### Benchmarks of other models (project 2, Task 6a)
`testsuite/lowrank/bench.sh SCHEME L N_delay [THREADS] [KEY=VALUE ...]`: extra keys are applied after the Hubbard defaults,
e.g. `NSWEEP=20 testsuite/lowrank/bench.sh delayed 8 32 1 'ham_name="tV"' 'Model="tV"' ham_V=-1.d0 N_SUN=1` (set `NSWEEP`;
the per-L table is calibrated for Hubbard). Such lines carry two more fields, `model` (ham_name) and `keys` (the extra
arguments joined with ','); the run directory gets `_<model>` after the scheme. `bench_report.py` accepts 10, 13 or 15 fields;
the existing reports regenerate byte-identically.

### Automatic N_delay in the harnesses (project 2, Phase 6)
Since Phase 6 the Start template (and the `&VAR_QMC` default) has `N_delay = 0`: automatic, static rule (deterministic;
`Prog/Lowrank_update_mod.F90`, header "Automatic N_delay"); `N_delay = -1` adds warm-up timing in sweeps 2-3; `N_delay >= 1`
is fixed. `info` says `Sequential scheme: <NAME>, N_delay = auto (static n0 = <n0>, Ndim = <Ndim>)` (or `auto+tune`, then
`Auto N_delay: candidates ... s/proposal, chosen <n>` or `Auto N_delay: tuning skipped (<reason>), N_delay = <n0>`).
- `bench.sh SCHEME L N_delay`: `N_delay` 0 / -1 are accepted (run directory `nd0` / `ndm1`); the output line then ends with
  one more field `nd_used=<n>` (the value the run used, from `info`; `na` for standard). `bench_report.py` accepts 10, 13, 14,
  15 or 16 fields (`nd_used` is required exactly for `N_delay <= 0`), labels those rows `auto` / `tune` (with the value used),
  keeps them out of the best-fixed and "Default N_delay" computations and, in `--ndelay` mode, lists them against the best
  fixed value of the same L and scheme. The committed reports regenerate byte-identically.
- Audit of the old default 32: `compare_trajectory.sh`, `tiers.sh` (traj), `regress_vs_baseline.sh` and `parallel_tier.sh`
  check pass/fail properties that hold for any N_delay (same decisions as standard, bit-identity of the standard path,
  statistical agreement, same-seed identity; the static value is the same on every rank and thread count), so they now run
  the automatic default and exercise it (`tiers.sh` keeps tv/kondo/jz at the fixed 2 and 1000 and z2gt at 5).
  `stat_runs.sh` and the `debug` phase of `tiers.sh` pin `N_delay=32` (recorded statistical/stability and shadow-error
  numbers, existing chains in `runs/stat`). `bench.sh` always passes N_delay explicitly; the regression set `explicit` keeps
  `N_delay=7` with `standard` (ignored there).
