# AGENTS.md — rules for AI agents working on ALF (branch `Code-modularity` and its children)

Read this file first. It says what goes where and how to work. Details are linked, so you do not
need to read the whole code base.

| You need | Read |
|---|---|
| What the current phase is, what the owner wants | [`docs/ai/phase2/OWNER_BRIEF.md`](docs/ai/phase2/OWNER_BRIEF.md), [`docs/ai/phase2/ROADMAP.md`](docs/ai/phase2/ROADMAP.md) |
| Which module does what, call graph, where data lives | [`docs/ai/CODEMAP.md`](docs/ai/CODEMAP.md) |
| How to branch, commit and merge | [`docs/ai/BRANCHING.md`](docs/ai/BRANCHING.md) |
| Build commands, tests, CI jobs | [`docs/ai/codemap/build.md`](docs/ai/codemap/build.md) |
| Tools that help you navigate the code | [`docs/ai/TOOLS.md`](docs/ai/TOOLS.md) |
| Known bugs and their status | [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) — read before working on stabilisation, CGR or large N |
| Why the code looks like this | origin ALF-QMC/ALF PRs #656 (Phases 0–5) and #657 (Phase 6): `gh pr view 656 -R ALF-QMC/ALF` |

## Working mode

- **Planning:** ponytail off. **Coding:** use the `ponytail` skill; a subagent that writes code
  uses it at level `lite` (owner rule, 2026-11-25).
- Every change that goes into `Code-modularity` uses a dev branch and a clean canvas branch
  ([`docs/ai/BRANCHING.md`](docs/ai/BRANCHING.md)). Never commit work-in-progress to `Code-modularity`.
- Never push to origin (ALF-QMC/ALF). The working remote is the fork (`luisthremer/ALF`).
- Report what you ran and what you did not run. "Bitwise identical" is claimed only after the
  check in [`docs/ai/codemap/build.md`](docs/ai/codemap/build.md) §2 (or the script it names).
- **Bugs found along the way** (owner rule, 2026-11-26): classify every bug at once and record it in
  [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) with class and status (open / fixed in <merge>) before anything else.
  - **A** — silently wrong results in the code path being worked on or in a configuration used in production
    (the run ends normally but the numbers are wrong, or a guard that hides such failures): interrupt the current
    step, write the entry, warn the production users (the session or card that runs production) at once, fix it as
    its own small dev/canvas feature, merge, return to the step.
  - **B** — a real bug outside the current code path (wrong results or crashes elsewhere, hangs, wrong measurements
    of features nobody uses now): entry + one line on the whiteboard card; fixed in one batch feature after each
    phase merge.
  - **C** — cosmetic or theoretical (messages, limits far beyond any size we run, small inefficiencies): entry
    only; fixed when that code is touched anyway.
  - Limits: one interruption at a time (an A found during an A fix comes next, no stack); at most two unmerged
    feature branches (the current step plus one A fix); unsure between A and B -> ask the owner, treat it as B
    meanwhile.
    While CI runs on a canvas, the next canvas may already be prepared (owner, 2026-11-26); merges stay one
    after another.
- **Status and ideas** (owner rule, 2026-11-26) live in `manager_tools/` on `Code-modularity`, `Code-modularity-CPU`
  and `Code-modularity-GPU` of the fork, the same on all three (public, owner's OK 2026-11-26; until then the private
  repo `luisthremer/alf-ops`). Commit them directly (they are not code: no dev/canvas) with `[skip ci]` in the
  subject on any of the three, then run `./manager_tools/sync.sh` to copy the folder to the other two:
  - `manager_tools/STATUS.md` is the only status page, rewritten at every step in the owner's format (Now / Waiting for you /
    Done since the last owner message / Next 3 steps / Open branches and running jobs), one screen at most, details
    only as links. Always `git pull --rebase` before writing, never force-push. Decisions that block the work are
    also sent through the calendar-invite channel.
  - `manager_tools/IDEAS.md` is the owner's inbox: read it at the start of a session, after every merge or finished
    step, and at least every two hours of work (e.g. while waiting for CI) (owner, 2026-11-26: "make sure the inbox
    is checked from time to time");
    triage each new line (fits a phase or CG-ALF idea -> add it there with a link; new research question -> own entry
    with what would decide it; unclear -> question under "Waiting for you"), then append "-> taken: <where> (<date>)".
    Never edit or delete the owner's lines. Ideas sent in the chat go there too unless the owner writes "urgent".

## What goes where

| Kind of code | Place | Rule |
|---|---|---|
| Orchestration of a run | `Prog/main.F90` | Calls only. No logic, no new state. |
| One-time setup, info file | `Prog/qmc_setup_mod.F90` | `Setup_model`, `Setup_fields`, `Setup_update_schemes`, `Write_info_*`. |
| Per-walker state (G, UDV stacks, phase) | `Prog/qmc_walker_mod.F90` (`qmc_walker_t`) | Large arrays are allocated once per run in `Walker_alloc`, not per call. |
| Sweeps | `Prog/qmc_sweep_mod.F90` | `Sweep_up`, `Sweep_down`. |
| Langevin / HMC / global MALA steps | `Prog/qmc_driver_mod.F90` | |
| Propagation of G through a slice, local updates | `Prog/Wrapgr_mod.F90` (`WRAPGRUP`, `WRAPGRDO`) | Local update: `Upgrade2` (`upgrade_mod`), or delayed/submatrix (`lowrank_update_mod`, selected by `Sequential_scheme`). |
| Stabilisation point | `Recompute_GR_phase` in `qmc_walker_mod` | The one routine a GPU stabilisation replaces. |
| Run-time parameters (`&VAR_QMC`, `&VAR_TEMP`) | `Prog/QMC_runtime_var_mod.F90` | Declaration, namelist, default, broadcast, `get_`/`set_` pair, `public`; also `Prog/default_parameters_generic.py`. Settings are locked after `check_update_schemes_compatibility`. |
| Model state (`Ndim`, `N_FL`, `Op_V`, `Op_T`, `WF_*`, …) | `Prog/Hamiltonian_main_mod.F90` | `PUBLIC, PROTECTED`: only Hamiltonian submodules write it. Tests use `Set_model_for_testing` + `Setup_flavor_map`. |
| A model | `Prog/Hamiltonians/Hamiltonian_<Name>_smod.F90` + `Prog/Hamiltonians.list` | Copy `Hamiltonian_##NAME##_smod.F90`. See [`codemap/model.md`](docs/ai/codemap/model.md) §2. |
| MPI communication | `Libraries/Modules/comm_mod.F90` | Use `is_root`, `comm_reduce_*`, `comm_bcast`; add an overload there if a type is missing. No new `#ifdef MPI` in `Prog/`. |
| BLAS / LAPACK | `Libraries/Modules/linalg_mod.F90` (`la_z*`) | No direct `call zgemm(...)` in `Prog/` (CI fails). Add a wrapper there if one is missing. |
| Errors | `runtime_error_mod`: `Terminate_on_error(ERROR_*, __FILE__, __LINE__)` | No `STOP` / `ERROR STOP` (CI fails). |
| A new core source file | `Prog/` (flat, no sub-folders) + `SRCS` in `Prog/Makefile` | Dependencies are generated by `gen_deps.py`. |
| Unit tests | `testsuite/<Area>.tests/`, registered in `testsuite/CMakeLists.txt` | Every new routine with logic gets a test; tests must build in all STAB/LOG and MPI variants. |
| User-visible changes | `CHANGELOG.md` | Only backward-incompatible changes and critical bugs. |
| Algorithms | `Documentation/*.tex` | `sampling.tex` Alg. 1 names the file of each QMC step. |

## Code rules (all checked by CI or by review)

1. Every `USE` in `Prog/` (outside `Prog/Hamiltonians/`) has an `ONLY:` list. Exceptions: `mpi`,
   `hdf5`, `h5lt`, `alf_hdf5`.
2. No implicit `SAVE`: never initialise a local variable in its declaration. Use `parameter` or an
   assignment.
3. No `STOP`; use `Terminate_on_error`.
4. BLAS/LAPACK only through `linalg_mod`; MPI only through `comm_mod`.
5. Do not change the interfaces of existing routines, the Hamiltonians, the input files or the
   output unless that is the task. A refactor must leave results bitwise identical; a new feature
   must leave the default path bitwise identical.
6. Keep the style of the file you edit (Fortran 2008, free form, `Kind=Kind(0.d0)`, `Implicit none`).
7. Doxygen header on every new module and public routine (`!> @brief`, `!> @details`).

## Commits (on canvas branches)

- One task per commit. Subject: `<area>: <new state>`, no period (e.g.
  `lowrank: delayed scheme behind Sequential_scheme`). Inside a numbered plan: `<step>: <area>: …`.
- Body: what moved where, why, why results can or cannot change, which test covers it, and what
  was left alone. End with the verification line ("Bitwise identical to Code-modularity." or the
  measured agreement) and the `Co-Authored-By` trailer.

## Checks before you say "done"

```bash
source configure.sh GNU noMPI NO-INTERACTIVE && make -j8 all
mkdir -p testsuite/tests && cd testsuite/tests && cmake .. && make -j8 && ctest -j6
```
Repeat with `MPI` (and `Tempering`/`PP` when MPI code changed). Run the three convention checks of
`.github/workflows/tests.yml` (job `report_use_only`, `report_save_state`) locally or wait for CI.
