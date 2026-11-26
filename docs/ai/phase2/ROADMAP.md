# GPU Development Phase 2 — roadmap

Source: the owner's brief [`OWNER_BRIEF.md`](OWNER_BRIEF.md) (2026-11-25). Every phase ends with a
merge that follows [`../BRANCHING.md`](../BRANCHING.md) (dev branch → clean canvas branch →
merge, CI green, archive tags). Every line of code follows the Code-modularity rules in
[`../../../AGENTS.md`](../../../AGENTS.md).

**Working mode for every session of this phase:** plan with ponytail off; code with the
`ponytail` skill (subagents that write code use it at level `lite`).

Status values: `open`, `running`, `done (<date>, <commit>)`.

## Base line

| Branch (fork luisthremer/ALF) | Meaning |
|---|---|
| `master` | = origin ALF-QMC/ALF `master` (`4c690b04`) |
| `Code-modularity` | Fakher's Phases 0–6 + 6G (origin PRs #656, #657; tip `af22b976` on 2026-11-25) + our merges |
| `Code-modularity-CPU` | from step P2.3 on: CPU line |
| `Code-modularity-GPU` | from step P2.3 on: GPU line |
| `cm/dev/<x>`, `cm/feat/<x>` | dev and canvas branch of a feature `<x>` (see BRANCHING.md) |
| `archive/cm/<x>/{dev,feat}` | tags left behind after a merge |

## Rule: bugs found along the way (owner, 2026-11-26)

Every bug is classified at once and entered in `KNOWN_ISSUES.md` with class and status:
**A** (silently wrong results in the current code path or a production configuration, or a guard hiding them) ->
interrupt, entry first, warn production users, own small dev/canvas fix, merge, back to the step;
**B** (real bug outside the current path) -> entry + whiteboard line, fixed in one batch feature after each phase
merge; **C** (cosmetic/theoretical) -> entry only, fixed when the code is touched anyway.
At most one interruption at a time and at most two unmerged feature branches (current step + one A fix);
unsure between A and B -> ask the owner, treat as B meanwhile. Full text: `AGENTS.md`, Working mode.
Status page and idea inbox: `manager_tools/STATUS.md` and `manager_tools/IDEAS.md` on the three main lines (until
2026-11-26 the private repo `luisthremer/alf-ops`), see `AGENTS.md`.

## P2.0 — Sync (done 2026-11-25)
Fork `Code-modularity` = origin `Code-modularity_Phase6` (`af22b976`), a fast-forward of
`Code-modularity`. Fork `master` = origin `master`.

## P2.1 — Knowledge files for AI sessions (done 2026-11-25, merge `c771d98a`)
`AGENTS.md`, `CLAUDE.md`, `docs/ai/CODEMAP.md`, `docs/ai/BRANCHING.md`, this folder.
Feature name `ai-docs`. Done when merged into `Code-modularity` via `cm/feat/ai-docs`.

## P2.2 — Delayed and submatrix updates, native (done: merged 2026-11-25, `2bb0258a`; step 2 "fastest" merged 2026-11-26, `eaf94e27`)
Feature name `lowrank`. Source: fork `master_ai_cpu` (`8b18b95c`, LOWRANK in hook form:
`Prog/Lowrank/*`, `testsuite/Lowrank.tests/*`, `Documentation/lowrank*`).
1. `cm/dev/lowrank`: port without hooks. Scheme parameters (`Sequential_scheme`, `N_delay`) in
   `QMC_runtime_var`; scheme calls inside `WRAPGRUP`/`WRAPGRDO` (no duplicated copies of them);
   pending state in an owned type (no module-level allocatables); MPI only through `comm_mod`;
   BLAS/LAPACK only through `linalg_mod` (add wrappers there, e.g. `zgetrf/zgetrs`, with test);
   `USE, ONLY`; no implicit `SAVE`.
2. Then "absolute fastest" (owner decision): submatrix decides, the rank-k GEMM flush updates.
   Extend submatrix where it falls back today, one extension per task, each with test + timing:
   (a) non-diagonal vertices: replace the `K_ND = 1` guard by a representation without the
   cancellation loss, or switch that block to delayed; (b) rank-N_P > `N_P_LIM`; (c) block size
   rule per machine (static + timed). Delayed is used only where timing shows it is faster.
   Status 2026-11-25 (feature `lowrank-fastest`): measured on a Mac (GNU, 1 core): Hubbard 16x16
   b=6 standard 12.6 s, delayed 10.1 s, submatrix 9.8 s; tV 16x16 12.3 / 10.3 / 12.3 s; Kondo 12x12
   b=6 4.5 / 2.7 / 3.8 s. (a) is done as `Sequential_scheme = "auto"` (submatrix if every eligible
   vertex is diagonal, else delayed): always the faster of the two. A mixed submatrix/delayed block
   is postponed (low priority): it could gain at most the few percent by which submatrix beats
   delayed, and the factor order makes it numerically delicate. The larger time shares are
   measurements and stabilisation (P2.9).
3. Gates: all builds + ctest; standard scheme bitwise identical to `Code-modularity` on
   `alf_test_branch`; delayed/submatrix agree with standard to roundoff (`compare_schemes.py`)
   and with ED (`test_vs_ed`); fork CI green.
4. `cm/feat/lowrank`: one task per commit, tree equal to the dev tip, CI green, merge, archive.

## P2.3 — Split CPU / GPU (done 2026-11-26, from the merge of `cm/feat/p23-split`)
`Code-modularity-CPU` and `Code-modularity-GPU` both from the `Code-modularity` tip after P2.2
(that tip carries P2.2, the fixes det-phase and nan-guards, and the CG-ALF plan).
`Code-modularity` stays at the split point; new canvases target `-CPU` or `-GPU`.
The CPU line takes later CPU-only work; the GPU line merges the CPU line regularly
(`Code-modularity-CPU` → `-GPU`, never the other way round).

## P2.4 — GPU port onto Code-modularity-GPU (open; plan: [`P2_4_GPU_PLAN.md`](P2_4_GPU_PLAN.md))
Existing GPU work to bring over (fork, state 2026-11-25):
- `gpu-port` — general GPU code, backends `GPU=nvidia|intel|amd`, macros `ALF_LA_CUDA`,
  `ALF_LA_ONEMKL`, `ALF_LA_HIP`; docs in `docs/gpu/` (README, FEATURES, NOT_IMPLEMENTED,
  BUILD, MACHINES, BRANCHING).
- `gpu-port-unified` (+ `-tests`, `-maxent`) — point-defect model on top of it.
- `gpu-delayed-submatrix-defect-production*` (`-mpi`, `-mpi-amd`, `-mpi-pp`, `-mpi-pp-opt`,
  `-mpi-pp-helma`) — production variants; `fix/amd-checkerboard`; `gpu-port-apple-metal` (experimental).
Sub-phases, each a dev/canvas pair under `cm/dev/gpu-*`:
1. **Common scaffold** (`gpu-common`): device data layout of `qmc_walker_t`, backend dispatch
   inside Fakher's `linalg_mod` (the 6G wrappers are exactly the seam), persistent device
   workspaces, `GPU=none` bitwise identical. No vendor code yet.
2. **Per vendor** in parallel, each its own dev branch from the scaffold:
   `gpu-nvidia` (A30, H100/Helma), `gpu-intel` (Max 1550/SuperMUC-NG Phase 2, B580),
   `gpu-amd` (MI210). Sub-branches allowed per independent part.
3. **Integrate** the vendor branches on `cm/dev/gpu-integrate`; all three builds + CPU CI green;
   GPU vs `GPU=none` trajectory checks per machine (tolerances from `docs/gpu/`).
4. **Consolidate**: canvas `cm/feat/gpu` with one task per commit, tree equal to the integrate
   tip, tested on every GPU machine, merge into `Code-modularity-GPU`, archive all branches.

## CPU-line queue (D4, owner 2026-11-26)
gpu-port's unguarded CPU-path changes and other general features, each its own dev/canvas feature on
`Code-modularity-CPU` (then merged into the GPU line), in this order unless the owner says otherwise:
1. `ALF_UNSTABLE_THRESHOLD` (owner: first in the queue; small).
2. `prepare_sim.py`: general PARALLEL_PARAMS helper (Q4a), newest version from gpu-port
   `Scripts_and_Parameters_files/Start_PP/` (327 lines), with the Helma documentation.
3. `alf_exit` (clean exit path), `Trace_updates_mod`.
4. `Timers_mod` / `ALF_TEAM_TIMERS` as base of the compile-time timing flag (CG-ALF T8).
5. MaxEnt `N_boot` (with the MaxEnt branch).
6. `QR_PREPIVOT` variants 2/3 (variant 1 is a test-only flag in P2.4 1d/1e).
Elsewhere: `Op_V_is_symmetric` rewrite = CG-ALF T3; `Lowrank_cross_slice` = CG-ALF T6 (research branch now).

## P2.5 — Feature assessment (open, after P2.4)
Table in `docs/ai/phase2/GPU_FEATURES.md`: every ALF feature (models, observables equal-time and
time-displaced, global moves, Global_tau, Langevin/HMC/MALA, tempering, PARALLEL_PARAMS,
projector/trial wave function, STAB1/2/3/LOG, HDF5, symmetric Trotter, N_SUN, complex phases,
time-dependent hopping, …) × {works on GPU, does not work, untested} per vendor, each entry with
the test or run that shows it. Start from `docs/gpu/FEATURES.md` and `NOT_IMPLEMENTED.md`.

## P2.6 — Port the missing features (open)
One dev/canvas pair per independent feature (`cm/dev/gpu-feat-<name>`), sub-branches for
independent sub-features, consolidated and merged into `Code-modularity-GPU` one by one with CI
green and the GPU checks passed. The table of P2.5 is updated at every merge.
First (owner, Q4d): delayed and submatrix updates for every model (Kondo, tV, Spin_Peierls: non-diagonal and rank > 1
vertices, `Calc_Fl` on the device; ticket 407 / "401-A", tag `archive/tickets/gpu-port-delayed-submatrix-rank-nd`),
then a benchmark of which update scheme is fastest per model. Consistency across models has priority.

## P2.7 — Archive discipline (every phase)
After each merge: annotated tags `archive/cm/<x>/dev` and `archive/cm/<x>/feat` whose message says
what the branch did, which commit it was merged in, and where the tests/results are; a line in
`docs/ai/BRANCHING.md` (archive table). Merged `cm/dev`/`cm/feat` branches are deleted once their archive tags are on the fork (rule in BRANCHING.md); other branches only after the owner's yes.

## P2.8 — Break ALF: adversarial bug hunt (open, after P2.6, before CG-ALF)
Goal: find bugs before the speed work of P2.9 builds on them. Attack ALF (CPU and GPU lines) the
way a hostile user would, and record every finding in `KNOWN_ISSUES.md` (repo root). The model for
what we look for is the silent NaN phase of `fix/det-phase-underflow` (0b3ffba): return code 0,
precision checks fine, results wrong.
Rules: model-agnostic attacks only; every finding gets a minimal reproducer, an entry in
`KNOWN_ISSUES.md`, then its own dev/canvas fix branch with a regression test; fixes follow BRANCHING.md.
Attack list (each a task; split into parallel subagents):
1. **Size limits:** N from 1 to the memory limit (the underflow class: products of N numbers,
   exp/log ranges, integer overflow in N^2 index arithmetic, 32-bit sizes in BLAS/MPI calls).
2. **Parameter edges:** Beta -> 0 and large, Dtau large, Ltrot = 1, Nwrap = 1 and > Ltrot, U = 0
   and huge, N_SUN, N_FL, Nt_sequential ranges, N_delay extremes, Ltau with Projector, Theta = 0.
3. **Input fuzzing:** malformed or missing namelists, wrong types, unknown keys, absurd values,
   corrupted `confin`/`seeds`/HDF5 restart files: ALF must stop with `Terminate_on_error` and a
   clear message, never run on silently.
4. **Silent-failure guards:** NaN/Inf anywhere in G, phase, weights, observables; build with
   `-ffpe-trap=invalid,zero,overflow`, `-fcheck=all`, valgrind, AddressSanitizer.
5. **Differential testing:** standard vs delayed vs submatrix (byte-identical configurations),
   CPU vs GPU, MPI ranks 1/2/3/odd, Tempering/PP, STAB1/2/3/LOG, compilers, restart vs single run.
6. **Concurrency:** MPI rank/group layouts, CPU_MAX truncation at every point of a bin, OpenMP
   thread counts, GPU streams.
7. **Statistics:** long runs against ED (`test_vs_ed`) for every model; error bars must cover the
   exact value at the expected rate.
First switch the plugin back on (disabled 2026-11-26 until P2.8): `claude plugin enable claude-bughunter@elementalsouls`.
**Subagent with Claude-BugHunter** (owner instruction): one subagent loads the plugin
`claude-bughunter@elementalsouls` and applies its hunting method to ALF — `bb-methodology`
(5-phase workflow), `hunt-exceptional-conditions` (malformed input, fail-open), `hunt-race-condition`
(MPI/OpenMP), and `triage-validation` (its 7-question gate) to rate every finding before it enters
`KNOWN_ISSUES.md`. The plugin is written for web applications; its web-specific skills do not apply.
Done when: every attack task ran, every finding is in `KNOWN_ISSUES.md` with status, and every
fixed finding has a regression test in CTest.

## P2.9 — ENDPHASE: CG-ALF ("Crazy good ALF") (open, after P2.8)
R&D branches `rd/<topic>` (free to experiment; anything that survives goes through dev/canvas).
Loop: brainstorm → measure → implement → validate → iterate. Model-agnostic only; ALF stays modular.
1. **Theoretical limit.** Per GPU (A30, H100, MI210, Max 1550): peak FP64/FP32/tensor FLOPs and
   bandwidth; roofline of every hot kernel (flush GEMM, wraps, stabilisation QR, measurements).
   Target = distance to the roofline, not "GPU utilisation".
2. **Real utilisation.** Nsight Compute/Systems, rocprof, VTune/unitrace: SM occupancy, achieved
   FLOPs, stall reasons, launch gaps, host waits. Remove idle time (kernel fusion, graphs/streams,
   fewer host syncs, batched small ops).
3. **Precision.** Where FP32 (or FP32 split / tensor cores) is safe: candidate kernels, validation
   protocol against FP64 (observables and error bars, not only bitwise), opt-in switches.
4. **Numbers module** (start from Florian Goth's existing code): compile-time choice of real vs
   complex (real HS fields, real hopping → real arithmetic), integer hopping where it pays.
5. **Updates and stabilisation research:** newer stabilisation algorithms, larger submatrix
   blocks on GPU, fewer QR per sweep, mixed-precision stabilisation with checks.
6. **Throughput:** OpenMP inside a node, many independent Markov chains after warm-up (walker
   batches on one GPU), MPS/MIG sharing, multi-GPU.
7. **Largest system sizes:** memory footprint per walker, out-of-core/recompute trade-offs.
Each idea gets a one-line entry in `docs/ai/phase2/CG_ALF_IDEAS.md` with status and measured gain.
Plan, experiments E0-E4 and their results so far: [`CG_ALF_PLAN.md`](CG_ALF_PLAN.md); GPU hours used: [`CG_ALF_LEDGER.md`](CG_ALF_LEDGER.md).
