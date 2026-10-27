# ALF GPU port – guide for agents

Branch `gpu-port` of a fork of ALF (Algorithms for Lattice Fermions, modern Fortran, GPL).
Goal: port the **finite-temperature DQMC** to GPUs with one source tree, two equal targets:

- **NVIDIA**: dev/test machine `ssh nividia` (hardware/toolchain: see `docs/gpu/MACHINES.md` once probed).
- **Intel**: SuperMUC-NG Phase 2 (2× Xeon 8480+, 4× Data Center GPU Max 1550 per node, oneAPI:
  ifx + oneMKL + Intel MPI). Not reachable by agents. We prepare build profiles, job scripts and
  a validation recipe; the owner runs them.

A build with `GPU=none` must behave **exactly** like upstream and must not get slower (GNU and Intel).
Later phases (only with the owner's go-ahead): delayed and submatrix updates (Sun & Xu), CPU first.

## Git

- `origin` = upstream `ALF-QMC/ALF` (fetch only; push URL is `no_push`). **Never push to origin,
  never open upstream PRs.** `fork` = `git@github-work:luisthremer/ALF.git` (public). Push freely.
- Work on `gpu-port`; one sub-branch per larger ticket if useful (`gpu-port/NNN-name`).
- Leave the fork's other branches (`master`, `luis_playground`) alone.
- Small, coherent, tested commits. Only the manager session commits and pushes; subagents never do.
- No secrets in code or commits.

## Where things live

| Path | Content |
|---|---|
| `Prog/` | DQMC program (main loop, wraps, `upgrade_mod.F90`, stabilization, Hamiltonians, observables) |
| `Libraries/` | `Modules/` (mat_mod, UDV, QDRP, …), `libqrref` |
| `Analysis/` | binning/error analysis |
| `configure.sh` | `source configure.sh MACHINE MODE STAB [HDF5] [Devel]`, then `make`. GPU switch `GPU=none\|nvidia\|intel` planned (Phase 1 ADR) |
| `testsuite/` | CMake unit tests (`matmod.tests`, `Prog.tests`, `Analysis.tests`), `test_vs_ed/` (pyALF), `test_branch_parameters.json` + `compare_dirs.py` |
| `docs/gpu/` | port documentation: `README.md` index, `HARD.md`, `PLAN.md`, `ADR-*.md`, `BENCHMARKS.md`, `FEATURES.md`, `DISCREPANCIES.md`, `tickets/` |
| `.claude/agents/` | `fortran-gpu-dev`, `test-runner`, `numerics-reviewer` |

Research leads (not ground truth; verify before relying on them, log mismatches in
`docs/gpu/DISCREPANCIES.md`): `/Users/luis/Documents/_WORK/Coding/CLAUDE/ALF_GPU_research/R{1,2,3}_report.md`.

## Build and test (CPU) – verified on nividia 2026-10-26

On `nividia`, environments live in `~/alf-gpu/env/{gnu,intel}.sh` (GNU 13.3 + OpenBLAS + Intel MPI's
gfortran wrapper; ifx 2025.3 + MKL + Intel MPI) plus the venv `~/alf-gpu/venv` (pyALF, PyYAML).
Helper scripts there: `build_ctest.sh`, `build_prod.sh`, `run_ed.sh` (see `docs/gpu/BASELINE.md`).

```sh
. ./configure.sh GNU devel noMPI NO-INTERACTIVE      # unit-test build, as in CI
make -j8 all                                          # 'all', not 'program': Analysis.tests need Analysis/
cd testsuite && mkdir -p tests && cd tests && cmake .. && cmake --build . -j8 && ctest -j8
. ./configure.sh GNU MPI [LOG] HDF5 NO-INTERACTIVE    # production; IntelLLVM for ifx
```
- Use **HDF5** for anything pyALF analyses: pyALF 2.6 cannot parse gfortran's plain-text complex
  output (`(a, b)` with a space). configure.sh self-builds HDF5 1.14.6 into `HDF5/<compiler>`.
- `test_vs_ed`: upstream `prep_test.py` + `analysis.py`, simulation via `run_ed.sh` (2 MPI ranks per Δτ).
- CPU references for the port: `testsuite/gpu_reference/` (`configs.py`, `prep.py`, `run.sh`,
  `summarize.py`, `compare.py` – provisional, see BASELINE.md); results summary in
  `docs/gpu/BASELINE.md`, full data on nividia `~/alf-gpu/runs/ref_*`. Copies of the nividia
  helper scripts: `docs/gpu/scripts/nividia/`.
- Branch comparison: run `test_branch_parameters.json` on reference and test branch, compare with
  `compare_dirs.py`.

## Correctness rules (non-negotiable)

- **Kernel level:** fixed HS fields (deterministic), GPU vs CPU for G after a wrap, after
  stabilization, after a sweep segment: relative agreement ≈ 1e-10 (documented per test).
- **Observable level:** energies, density/spin correlations, equal-time and τ-resolved, agree with
  the CPU reference within ≥3σ; report χ². `test_vs_ed` passes on GPU. `Control_PrecisionG` and
  phase precision stay at CPU levels and are reported.
- RNG streams differ between compilers: never compare Monte Carlo trajectories bin by bin across builds.
- **Never weaken a tolerance or skip a test to get green.** If a criterion can't be met: stop, explain.
- When research, docs or a subagent report conflicts with a measurement, trust the measurement and
  record the discrepancy.

## GPU design constraints (from the brief; ADR in Phase 1 decides the rest)

One source tree; GPU selected at configure time; no polymorphism or derived-type deep copies on the
device (plain contiguous arrays only); data stays device-resident across the sweep; every
host↔device transfer is documented. Target the default and `STABLOG` stabilization paths.

## Machine rules (`nividia`)

- **nividia is for testing only. All development happens on the owner's laptop** (this repo); code
  goes to nividia by rsync from the laptop or by cloning/pulling `fork`. Never edit source on nividia.
- **Confinement:** write only inside `~/alf-gpu/`. A ticket agent writes only inside its own
  `~/alf-gpu/work/tNNN/` (plus result binaries in `~/alf-gpu/bin/` where its ticket says so). Set
  `export TMPDIR=~/alf-gpu/work/tNNN/tmp` (configure.sh/HDF5 and CUDA tools write temp files) and
  never use `/tmp` directly.
- Never delete or modify anything outside your own directory, never touch other users' files or
  processes, never kill a process you didn't start. No `rm -rf` on paths built from variables or
  globs that could expand outside your directory.
- Each test/benchmark run ≲ 15 min; ask the owner before anything longer.
- Clean up processes and large files afterwards; don't hold the GPU between tasks.
- Install nothing system-wide without asking. User-space installs (e.g. oneAPI) only after the owner agrees.

## Tickets and delegation

- Tickets: `docs/gpu/tickets/NNN-short-name.md` (template `000-template.md`): goal, files allowed,
  acceptance (tests + tolerances), out of scope.
- One ticket per subagent, with all needed context. Independent tickets run in parallel.
- Every diff is reviewed (manager or `numerics-reviewer`) and acceptance tests are re-run by the
  manager before a ticket closes.
- Stuck after 2–3 different approaches: write down what was tried and ask the owner.
- Phases end at gates (owner review): summary, numbers, plan changes, open decisions, next tickets.
