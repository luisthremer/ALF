# Branching rule: dev branch, canvas branch, merge

**Rule (owner, 2026-11-25).** Every change that is merged into `Code-modularity` (and later into
`Code-modularity-CPU` / `Code-modularity-GPU`) is written twice:

1. **Dev branch** `cm/dev/<feature>` — the real work. Commit as often as you like, iterate, fix,
   until the code is correct and CI is green.
2. **Canvas branch** `cm/feat/<feature>` — made fresh from the target branch after the dev branch
   works. Re-implement the change as a sequence of tasks, **one task per commit**, in the order
   a reviewer best understands it (a "red thread" from the start to the merge). Commits may be
   large, but each does exactly one thing. The canvas is written with full knowledge of the final
   code, so it can be cleaner than the history of the dev branch.

Only the canvas branch is merged.

## Before the merge

- [ ] `git diff cm/dev/<feature> cm/feat/<feature>` is empty (same files, same content). If the
      canvas is deliberately better (a comment, a name), port that back to the dev branch first.
- [ ] Every canvas commit builds (`git rebase -x 'make -j8 program' <target>` or equivalent).
- [ ] Fork CI is green on the canvas tip.
- [ ] One more test after CI (the feature's own check, e.g. a run or the scheme comparison).
- [ ] Merge with `git merge --no-ff cm/feat/<feature>` into the target; message
      `Merge cm/feat/<feature>: <one line>`.

## After the merge: archive

```bash
git tag -a archive/cm/<feature>/dev  cm/dev/<feature>  -m "<what was tried here; merged as <merge sha>; results in <path>>"
git tag -a archive/cm/<feature>/feat cm/feat/<feature> -m "Canvas of <feature>, merged into <target> as <merge sha>"
git push fork archive/cm/<feature>/dev archive/cm/<feature>/feat
```
Add a row to the table below. Then delete both branches, on the fork and locally
(standing yes of the owner, 2026-11-26), but only when all of these hold:
- the canvas is merged (its tip is an ancestor of the target branch);
- `archive/cm/<feature>/dev` and `archive/cm/<feature>/feat` are on the fork and point exactly at the
  branch tips (`git ls-remote fork 'refs/tags/archive/cm/<feature>/*'`);
- no worktree has the branch checked out and no open pull request uses it as head.

```bash
git push fork --delete cm/dev/<feature> cm/feat/<feature>
git branch -D cm/dev/<feature> cm/feat/<feature>
```
Every other branch (R&D, production, backups, anything without an archive tag) is deleted only after
the owner says yes for that branch.

Several dev branches per feature are fine (`cm/dev/<feature>-<part>`); larger efforts may use an
integration branch (`cm/dev/<feature>-integrate`) before the single canvas.

## Archive

| Feature | Target | Merge | Dev tag | Canvas tag | What it did |
|---|---|---|---|---|---|
| ai-docs | Code-modularity | `c771d98a` (2026-11-25) | `archive/cm/ai-docs/dev` | `archive/cm/ai-docs/feat` | AGENTS.md, CLAUDE.md, docs/ai (code map, branching rule, tools, phase-2 brief and roadmap), KNOWN_ISSUES.md |
| lowrank | Code-modularity | `2bb0258a` (2026-11-25) | `archive/cm/lowrank/dev` | `archive/cm/lowrank/feat` | native delayed/submatrix update (`lowrank_update_mod`, `Sequential_scheme`, `N_delay`), tests, documentation |
| det-phase | Code-modularity | `aef50743` (2026-11-26) | `archive/cm/det-phase/dev` | `archive/cm/det-phase/feat` | CPU determinant phase without underflow for large N (DET_C), phase guard; CI green on `eed184aa` |
| known-issues | Code-modularity | `ad861fc9` (2026-11-26) | `archive/cm/known-issues/dev` | `archive/cm/known-issues/feat` | bug-hunt findings in KNOWN_ISSUES.md with classes A/B/C; owner rules for bugs, status and ideas (docs only) |
| nan-guards | Code-modularity | `f8dbc251` (2026-11-26) | `archive/cm/nan-guards/dev` | `archive/cm/nan-guards/feat` | class-A fix H1: NaN/Inf checks of Green functions and UDV scales by bit pattern; CI green on `5defa12a` (run 45785448080) |
| lowrank-fastest | Code-modularity | `eaf94e27` (2026-11-26) | `archive/cm/lowrank-fastest/dev` | `archive/cm/lowrank-fastest/feat` | P2.2 step 2: `Sequential_scheme = "auto"` (submatrix if every eligible vertex is diagonal, else delayed), timings; CI green on `0865524b` (run 45791629146) |
| cg-alf-plan | Code-modularity | `c2efcce5` (2026-11-26) | `archive/cm/cg-alf-plan/dev` | `archive/cm/cg-alf-plan/feat` | CG-ALF plan (P2.9) with the results of E0-E4 on Helma H100, GPU-hour ledger (docs only; non-.md tree equals green `0865524b`) |
| p23-split | Code-modularity | `40b05e46` (2026-11-26) | `archive/cm/p23-split/dev` | `archive/cm/p23-split/feat` | ROADMAP: P2.2 and P2.3 done; split point of `Code-modularity-CPU` / `-GPU` (docs only) |
| branch-delete-rule | Code-modularity-CPU | `b64e88b2` (2026-11-26) | `archive/cm/branch-delete-rule/dev` | `archive/cm/branch-delete-rule/feat` | rule: merged and tagged `cm/*` branches are deleted right after the merge (docs only; `Code-modularity-GPU` fast-forwarded to it) |
| inbox-tasks | Code-modularity-CPU | `fc8561cf` (2026-11-26) | `archive/cm/inbox-tasks/dev` | `archive/cm/inbox-tasks/feat` | CG-ALF tasks T3 (fast `Op_V_is_symmetric`) and T4 (production guards) from the idea inbox; inbox read every two hours (docs only) |
| manager-tools | Code-modularity | `ef77fc51` (2026-11-26) | `archive/cm/manager-tools/dev` | `archive/cm/manager-tools/feat` | owner's folder `manager_tools/` (network graph tool, STATUS.md and IDEAS.md moved from alf-ops); docs point there |
| status-in-manager-tools | Code-modularity-CPU | `2cb655e2` (2026-11-26) | `archive/cm/status-in-manager-tools/dev` | `archive/cm/status-in-manager-tools/feat` | docs point to `manager_tools/` for status and ideas; plugin decisions (docs only) |
| gpu-plan | Code-modularity-CPU | `7782f8fe` (2026-11-26) | `archive/cm/gpu-plan/dev` | `archive/cm/gpu-plan/feat` | P2.4 plan `P2_4_GPU_PLAN.md`, probe of the GPU path on the host (docs only) |
| plan-decisions | Code-modularity-CPU | `038e442a` (2026-11-26) | `archive/cm/plan-decisions/dev` | `archive/cm/plan-decisions/feat` | P2.4 plan with the owner's decisions D1-D3 (docs only) |
| mt-all-lines | Code-modularity-CPU | `7868f455` (2026-11-26) | `archive/cm/mt-all-lines/dev` | `archive/cm/mt-all-lines/feat` | `manager_tools/` on all three lines, `sync.sh`; archive rows (docs only) |
| phase6-sync | Code-modularity-CPU | `095b21a6` (2026-11-26) | `archive/cm/phase6-sync/dev` | `archive/cm/phase6-sync/feat` | Fakher's `Code-modularity_Phase6` up to `fa21b3c0` (6G.E1-E7 linalg wrappers, STAB1 projector fix) |
| b-batch-1 | Code-modularity-CPU | `bfe5aeaf` (2026-11-26) | `archive/cm/b-batch-1/dev` | `archive/cm/b-batch-1/feat` | B bugs H3 H4 H6 H7 M1 M2 M3 M4 with regression tests |
| gpu-build | Code-modularity-GPU | `fe2aec86` (2026-11-26) | `archive/cm/gpu-build/dev` | `archive/cm/gpu-build/feat` | P2.4 step 1a: `GPU=none|host` build switch |
| cgalf-ideas | Code-modularity-CPU | `03457396` (2026-11-26) | `archive/cm/cgalf-ideas/dev` | `archive/cm/cgalf-ideas/feat` | CG-ALF ideas I11-I13 (docs only) |
| gpu-state | Code-modularity-GPU | `eb5b1022` (2026-11-26) | `archive/cm/gpu-state/dev` | `archive/cm/gpu-state/feat` | P2.4 step 1b: `gpu_state_mod`, device mirror of model and walker (Hubbard), CI job `gpu_host`; CI green on `984e3833` (run 45838961102) |
| gpu-wrap-update | Code-modularity-GPU | `0c222961` (2026-11-26) | `archive/cm/gpu-wrap-update/dev` | `archive/cm/gpu-wrap-update/feat` | P2.4 step 1c: `gpu_wrap_mod`, `gpu_update_mod`, `gpu_lowrank_mod` (device wraps, local updates, delayed/submatrix), Wrapgr hooks, `la_zgemm_cr`; CI green on `10d21073` (run 45854274445) |
| gpu-stab | Code-modularity-GPU | `36d17988` (2026-11-26) | `archive/cm/gpu-stab/dev` | `archive/cm/gpu-stab/feat` | P2.4 step 1d: `gpu_stab_mod` (device UDV stabilisation, LU phase as `LU_PHASE_SAFE`), `la_zqrpp`, `la_dznrm2`, `Control_PrecisionG_vals`, test-only `QR_PREPIVOT`; CI green on `fc0bcd0b` (run 45868504866) |
| owner-decisions | Code-modularity-CPU | `255251ff` (2026-11-26) | `archive/cm/owner-decisions/dev` | `archive/cm/owner-decisions/feat` | owner decisions of the QUESTIONS interview in the plans: D4, S1 (step 1f), CPU-line queue, CG-ALF T5-T11 (docs only) |
| p25-features | Code-modularity-GPU | `9ebe5897` (2026-11-26) | `archive/cm/p25-features/dev` | `archive/cm/p25-features/feat` | P2.5 GPU feature table, first draft (docs only) |
| cgalf-t6-t1 | Code-modularity-CPU | `e58a99bd` (2026-11-26) | `archive/cm/cgalf-t6-t1/dev` | `archive/cm/cgalf-t6-t1/feat` | CG-ALF T6 (cross-slice: sound, not faster on the CPU) and T1 (stabilisation survey) results (docs only) |
