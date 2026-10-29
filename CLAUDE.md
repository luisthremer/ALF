# ALF GPU port – guide for agents

Branch `gpu-port` of a fork of ALF (Algorithms for Lattice Fermions, modern Fortran, GPL).
Goal: port the **finite-temperature DQMC** to GPUs with one source tree, two equal targets:

- **NVIDIA**: dev/test machine `ssh nividia` (hardware/toolchain: see `docs/gpu/MACHINES.md` once probed).
- **Intel**: SuperMUC-NG Phase 2 (2× Xeon 8480+, 4× Data Center GPU Max 1550 per node, oneAPI:
  ifx + oneMKL + Intel MPI). Not reachable by agents. We prepare build profiles, job scripts and
  a validation recipe; the owner runs them. Intel test machine for agents: `ssh intel` (Arc B580
  workstation, not PVC; see `docs/gpu/MACHINES.md`).

A build with `GPU=none` must behave **exactly** like upstream and must not get slower (GNU and Intel).
Later phases (only with the owner's go-ahead): delayed and submatrix updates (Sun & Xu), CPU first.

## Git

- `origin` = upstream `ALF-QMC/ALF` (fetch only; push URL is `no_push`). **Never push to origin,
  never open upstream PRs.** `fork` = `git@github-work:luisthremer/ALF.git` (public). Push freely.
- Work on `gpu-port`; one sub-branch per larger ticket or lane (`gpu-port-NNN-name`; git cannot hold
  `gpu-port/NNN` next to the branch `gpu-port`). Rollback points: tags, see `docs/gpu/ROLLBACK.md`.
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
| `docs/gpu/VERIFICATION.md` | the verification standard (acceptance harness, tiers, evidence, review roles) |
| `docs/gpu/DECISIONS.md` | owner decisions, quoted and dated |
| `testsuite/gpu_reference/acceptance/` | acceptance harness: ship/build/run/driver/verdict scripts, `NNN.json` acceptance files |
| `.claude/agents/` | `fortran-gpu-dev`, `test-runner`, `numerics-reviewer`, `acceptance-verifier` |

Research leads (not ground truth; verify before relying on them, log mismatches in
`docs/gpu/DISCREPANCIES.md`): `R{1,2,3}_report.md` in the research folder next to this clone (outside the repo; the owner
has the path).

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
  When a criterion fails, never put relaxing it forward as the first or recommended option. Criteria live in the
  ticket's pre-registered acceptance file; changing one needs the owner's quoted, dated approval (in the file and in
  `docs/gpu/DECISIONS.md`) and fresh seeds afterwards. Never re-run a statistical test until it is green.
- When research, docs or a subagent report conflicts with a measurement, trust the measurement and
  record the discrepancy.

## GPU design constraints (from the brief; ADR in Phase 1 decides the rest)

One source tree; GPU selected at configure time; no polymorphism or derived-type deep copies on the
device (plain contiguous arrays only); data stays device-resident across the sweep; every
host↔device transfer is documented. Target the default and `STABLOG` stabilization paths.

## Machine rules (`nividia`, `intel`) – both are test sandboxes

- **Which machine:** Intel targets (`GPU=intel`, ifx offload, oneMKL/SYCL) are tested on `intel`, NVIDIA targets
  on `nividia`; CPU-only checks stay on `nividia`. On `intel` select the B580 with `ONEAPI_DEVICE_SELECTOR=level_zero:0`
  and don't run plain `sycl-ls` (driver abort, 110 MB crash report in `/var/crash`).
- **Both machines are for testing only. All development happens on the owner's laptop** (this repo);
  code goes there by rsync from the laptop or by cloning/pulling `fork`. Never edit source there.
- **Confinement:** write only inside `~/alf-gpu/`. A ticket agent writes only inside its own
  `~/alf-gpu/work/tNNN/` (plus result binaries in `~/alf-gpu/bin/` where its ticket says so). Set
  `export TMPDIR=~/alf-gpu/work/tNNN/tmp` (configure.sh/HDF5, CUDA and SYCL tools write temp files)
  and never use `/tmp` directly.
- Never delete or modify anything outside your own directory, never touch other users' files or
  processes, never kill a process you didn't start. No `rm -rf` on paths built from variables or
  globs that could expand outside your directory.
- **Never use `sudo`** or anything that needs root, even where the account could. No changes to system settings, groups, drivers or packages: ask the owner.
- Each test/benchmark run ≲ 15 min; ask the owner before anything longer. A T2 acceptance sequence of such runs
  (~2 h unattended) needs the owner's standing approval of the per-ticket GPU budget.
- **nividia GPU sharing:** one GPU process at a time (two time-sliced CUDA processes ran 13–15× slower each, 2026-10-28/29),
  except under the project's MPS daemon (approved 2026-10-29, DECISIONS.md): started and stopped only
  by the manager (the acceptance driver itself runs GPU jobs one at a time, without MPS), as the ordinary user, `CUDA_VISIBLE_DEVICES=0`, pipe/log directories under
  `~/alf-gpu/work/mps/`, at most 2 concurrent clients, never while another user's process is on the GPU; afterwards
  `echo quit | nvidia-cuda-mps-control` and check that no MPS process is left; after a killed client, quit the daemon.
  Measured 2026-10-29: two MPS clients each at solo speed, outputs bit-identical to a solo run.
  GPU runs of the acceptance driver take the lock `~/alf-gpu/gpu.lock`.
- **nividia cores and disk:** ≤ 4 cores per agent ; up to 3 single-threaded CPU runs may run alongside one GPU run: measured 2026-10-29, the GPU run took 66 s alone,
  66 s with 2 and 67 s with 3 CPU runs alongside, output bit-identical; keep ≥ 5 GB free on `/`
  (the acceptance driver refuses to start below); per ticket ≤ 1.5 GB live evidence.
- **Evidence retention:** never delete raw outputs (`info`, `data.h5`, `run.log`, `provenance.json`, verdicts) or
  binaries of an unclosed ticket; cleanup removes build trees only. At the owner's gate the manager archives them to
  `~/alf-gpu/archive/tNNN/` with a sha256 manifest.
- **Confinement exceptions (manager only):** `~/alf-gpu/backup/`, `~/alf-gpu/archive/`, `~/alf-gpu/refcache/` (CPU
  reference cache, written by the manager, read by everyone), `~/alf-gpu/work/mps_trial/`, `~/alf-gpu/work/mps/`,
  `~/alf-gpu/gpu.lock`, `~/alf-gpu/hdf5/` (prebuilt HDF5 for the harness), and the acceptance driver's attempts log
  `~/alf-gpu/archive/attempts/` (written by the driver for every party, including an implementer's T1).
- Clean up processes and large files afterwards; don't hold the GPU between tasks.
- Install nothing system-wide without asking. User-space installs (e.g. oneAPI) only after the owner agrees.
- `intel` is a shared workstation (other users, a desktop session): before a run check `uptime` and
  `ps -eo user,pcpu,comm --sort=-pcpu | head`; use ≤ 8 cores; if another user is busy, wait or ask.
  Its Arc B580 is a consumer Xe2 GPU: results show that ifx/oneMKL code runs on a real Intel GPU,
  but they are **not** PVC performance numbers (FP64 rate, memory, launch latency differ).

## Tickets and delegation

- Tickets: `docs/gpu/tickets/NNN-short-name.md` (template `000-template.md`): goal, files allowed,
  acceptance (tests + tolerances), out of scope.
- One ticket per subagent, with all needed context. Independent tickets run in parallel.
- **Acceptance follows `docs/gpu/VERIFICATION.md`:** a ticket closes only on a PASS `verdict.json` of the T2 harness
  run by the manager on a frozen candidate built from `git archive`, with all positive controls killed, an independent
  `acceptance-verifier` agreeing, and the `numerics-reviewer` approving code, acceptance file and scripts. The committed
  sha equals the accepted sha up to `docs/**` (`acc_check_commit.sh`). Logs without a matching provenance stamp are not
  evidence. Acceptance runs use `env -i` and refuse debug/fault variables outside declared control runs.
- **Roles:** implementers (`fortran-gpu-dev`, Sonnet by default; Opus per ticket for numerically delicate work on the
  manager's proposal) run tier T0 per iteration and T1 on their own frozen candidate; verification never runs in the
  implementer's context; corrections of a report go to a fresh verifier, not to the author. At ~250k tokens or a phase
  boundary an implementer writes a handoff note into the ticket and a fresh agent continues. Resume an agent only after
  an interruption; after a change to an agent definition, spawn new agents.
- **Lanes:** NVIDIA lane (2xx tickets, `nividia`, `~/alf-gpu/work/t2NN`) and Intel lane (3xx tickets, `intel`,
  `~/alf-gpu/work/t3NN`) under one manager, each on its own sub-branch with file ownership in `docs/gpu/PLAN.md`; the
  manager sequences edits of shared files; the `alf_la` interface is frozen while both lanes run.
- Stuck after 2–3 different approaches: write down what was tried and ask the owner.
- Phases end at gates (owner review): summary, numbers, plan changes, open decisions, next tickets.
