# Log of the AI features on `master_ai`

*Written with AI assistance on `master_ai`; not part of the ALF baseline. `CHANGELOG.md` is a
baseline file and stays as on `master`; entries for the AI features go here, in the same format.*

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
