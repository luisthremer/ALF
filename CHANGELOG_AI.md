# Log of the AI features on `master_ai`

*Written with AI assistance on `master_ai`; not part of the ALF baseline. `CHANGELOG.md` is a
baseline file and stays as on `master`; entries for the AI features go here, in the same format.*

### 2026-11-20 LOWRANK: strict same-chain check, all models and ALF features, CI on every build mode
Author: L. Thremer, with AI assistance

1. `compare_schemes.py` now requires the final configuration `confout*` (fields and random number generator
   state) of the delayed and submatrix runs to be byte-identical to the standard run with the same seeds, i.e. the
   same Markov chain, and keeps the 1e-8 relative bound for all other output. Exception: sequential MALA and
   Global_tau MALA, whose proposal contains the Green-function force (observables to roundoff, <= 4e-11).
   N_delay is run automatic and fixed to 3.
2. Scenarios (`scenarios.py`): every Hamiltonian with lattices, N_SUN / N_FL variants, field types 1-4, flavour
   reconstruction; Global_moves, Global_tau, Propose_S0, Nt_sequential ranges, Checkerboard / Symm, restart, MPI
   (2 ranks), TEMPERING, PARALLEL_PARAMS. Test-only models (`testsuite/Lowrank.tests/models`, added by
   `add_test_models.sh`): multi-site Jz vertices, 10-site vertices (standard fallback), continuous and discrete
   fields in one run. Result: 264 serial scenario/scheme runs, observables <= 3.5e-11, confout identical except MALA.
3. CI (`master-ai.yml`): MPI tests through the MPI launcher; Intel MPI compiler name for CMake; IntelLLVM unit test:
   the bound on the drift of the standard path now scales with cond(1+A) (1.7e-10 at cond 5e7 was above the fixed
   5e-11); GNU additionally with TEMPERING, PARALLEL_PARAMS and HDF5; full scenario set on GNU serial.
4. `baseline_guard.py` sets the model name in both namelists (LRC and Z2_Matter were skipped before).
5. Documentation: `Documentation/lowrank.tex` (four verification levels, compatibility table); module headers of
   `Prog/Lowrank` explain the low-rank forms of G and the call order per vertex.

### 2026-11-20 AI branch CI gates; LOWRANK checked with MPI and tempering
Author: L. Thremer, with AI assistance

1. New workflow `.github/workflows/master-ai.yml` (branches master_ai, master_ai_cpu, ai/**): hook lint,
   baseline guard against ALF master, LOWRANK build + tests with GNU, PGI, Intel and IntelLLVM (noMPI, MPI).
2. `baseline_guard.py` and `compare_schemes.py`: option `--launcher` (e.g. "mpirun -np 2").
3. Checked locally (GNU, MPICH): MPI with 2 ranks, 17 scenarios, delayed/submatrix vs standard max rel.
   difference 3.4e-11; parallel tempering with 2 replicas (U = 3, 5) max 2.8e-11.

### 2026-11-20 LOWRANK: timed automatic N_delay (`N_delay = -1`)
Author: L. Thremer, with AI assistance

1. New value `N_delay = -1` in `&VAR_LOWRANK`: starts from the static value n0 (`N_delay = 0`), times the
   candidates {n0/2, n0, 2 n0} round-robin over the time slices of sweeps 2 and 3, and keeps the one with the
   least wall time per proposal if it beats n0 by more than 5 %, else n0. With MPI the sums are reduced over
   `Group_Comm`, group rank 0 decides and broadcasts. The decision is written to `info`
   ("Auto N_delay: candidates ... chosen n"). The arrays are sized for the largest candidate.
2. No new hook: sweeps are counted inside the feature (the WRAPGRDO call with NTAU = 1 ends a sweep), so
   `Prog/main.F90` stays as on master and `HOOKS.md` is unchanged.
3. `N_delay` changes only the summation order: `compare_schemes.py --n-delay -1` shows the same chain as the
   standard update up to roundoff. `36-lowrank-unit` checks the timing logic and the error exits of
   `Lowrank_set_N_delay` (CTest entries `36-lowrank-unit-pending`, `36-lowrank-unit-capacity`).
4. Not in this version (it was in the research branch): the "tuning incomplete" line for runs shorter than
   three sweeps, and the resolved N_delay in `data.h5`.

### 2026-11-20 LOWRANK: sequential MALA, Global_tau moves, CPU branch master_ai_cpu
Author: L. Thremer, with AI assistance

1. The low-rank schemes now work together with sequential MALA updates and with Global_tau moves
   (`Global_tau_moves`, `Global_tau_MALA_moves`); the run no longer stops. MALA vertices and the
   Global_tau moves run exactly as in the baseline, after the pending updates are applied.
2. The hooks in `Prog/Wrapgr_mod.F90` pass all arguments of WRAPGRUP/WRAPGRDO and the three Global_tau
   routines of Wrapgr_mod to the feature; see `HOOKS.md`.
3. Tests: unit test `36-lowrank-unit` (kernels against an exact Green function) and
   `compare_schemes.py`, both in CTest for builds with `LOWRANK` (hook in `testsuite/CMakeLists.txt`).
4. Branch `master_ai_cpu` = ALF master + this feature, built from the cleaned CPU research branch
   `feature/delayed-submatrix` (timing, trace and debug code removed).

### 2026-11-19 Delayed and submatrix sequential updates (feature LOWRANK)
Author: L. Thremer, with AI assistance

New build switch `LOWRANK` and optional namelist `&VAR_LOWRANK` (`Sequential_scheme`, `N_delay`);
see the documentation, section "Delayed and submatrix sequential updates". Nothing changes for builds
without the switch, nor for runs with `Sequential_scheme = "standard"` (the default).

1. Code in `Prog/Lowrank/`; hooks in `Prog/Wrapgr_mod.F90` listed in `HOOKS.md`.
2. Not supported together with a low-rank scheme (the run stops): sequential MALA updates and
   `Global_tau_moves`.
3. Compared with the research version on the fork branch `gpu-port`: the cross-slice delay prototype,
   the warm-up timing of `N_delay` (`N_delay = -1`), the debug shadow checks and the provenance
   output at the end of the run are not included; the parameters moved from `&VAR_QMC` to
   `&VAR_LOWRANK`.

Also new on `master_ai`: `HOOKS.md`, `testsuite/hook_lint.py` (diff of baseline files = listed
hooks) and `testsuite/baseline_guard.py` (bit-identical output against `master`).
