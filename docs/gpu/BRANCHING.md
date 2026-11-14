# Branching scheme for the ALF GPU repository

**Status: approved by the owner and carried out 2026-11-14. 40 archived branches are deleted on `origin` (each kept as a tag), 16 branches remain. The mapping table below shows the state before; the section "Carried out" shows what was done and what is left.**

## Scheme

Trunk-based, with one long-lived variant.

```
 feat/*  fix/*  exp/*        short-lived, branch off the trunk, merge back, then archive tag + delete
        \    /     /
 gpu-port ---o--o--o--o--o--o--o----------> trunk: all general GPU work (NVIDIA, Intel, AMD, MPI, tests)
              \        \        \
 gpu-port-unified --M--------M--------M---> variant: trunk + point defect; only merges the trunk
                                 |
   tags:  prod/<date>-<run>  (state a production run was built from)
          archive/<name>     (last commit of a finished or abandoned branch)
```

| Branch | Role | Rules |
|---|---|---|
| `gpu-port` | trunk | Every general change lands here. CI must be green. No force-push. |
| `gpu-port-unified` | the only long-lived variant (trunk plus the point defect) | `git merge gpu-port` on `gpu-port-unified`, never the reverse. Defect-only work is committed here directly or on a `feat/defect-*` branch off it. |
| `feat/<name>` | new feature, off the trunk | short-lived, merged into the trunk, then `archive/feat/<name>` tag and delete |
| `fix/<name>` | bug fix, off the trunk | same |
| `exp/<name>` | experiment that may be abandoned (prototype, measurement) | same; if abandoned, tag and delete without merging |
| tag `prod/<date>-<run>` | production state, e.g. `prod/2026-11-14-smuc-L32` | Set when a production run is submitted from a commit. Production runs are built from the tag, never from a moving branch. |
| tag `archive/<name>` | keeps every commit of a finished or abandoned branch | Created before the branch is deleted. Restore with `git branch <name> archive/<name>`. |
| tag `archive/tickets/<name>` | the same, for branches of the original ticket-driven development (2026-10-26 to 2026-11-07: `gpu-port-207-*`, `208-*`, `212-*`, `302*`, `harness-*`, `process-lanes`, `gpu-port-delayed-submatrix*`) | Owner's addition, 2026-11-14. Ticket status: `TICKETS_STATUS.md`. |

Day to day:

1. `git switch -c feat/x gpu-port`, work, push, CI green, merge into `gpu-port` (fast-forward or `--no-ff`, owner's choice).
2. `git switch gpu-port-unified && git merge gpu-port`, push.
3. Before a production submission: `git tag prod/<date>-<run> <commit>` and push the tag.
4. After the merge: `git tag archive/feat/x origin/feat/x`, push the tag, then delete the branch.

Notes:
- GitHub shows tags without CI. A tag does not move, so a production state cannot be changed by accident (this is what made `smuc2-prod-2026-11b` impossible to fix; it should have been a tag from the start, which is now the rule).
- The one-time `-s ours` step in `docs/gpu/BRANCHES.md` has been done; `gpu-port-unified` contains `gpu-port`.
- Data and plans live outside the repo; the whiteboard (`~/Documents/_WORK/WHITEBOARD`) shows branch state live.

## Mapping of every existing branch (origin, state 2026-11-14)

Counts are commits not in `gpu-port` / not in `gpu-port-unified`. "Tag as archive and delete" is safe even for unmerged commits because the tag keeps them; branches with unmerged commits are marked so the owner can look first. Branches with a worktree need the worktree removed before deleting. Names with `<date>` are for the owner to fill in (which commit a run was built from is not recorded in one place).

| Branch | Action | Target name / tag | vs trunk / vs unified | Last commit | What it is (from BRANCHES.md) |
|---|---|---|---|---|---|
| `gpu-port` | trunk | — | +0 / +0 | 2026-11-13 | general GPU port, no defect (main branch) |
| `gpu-port-unified` | variant | — | +12 / +0 | 2026-11-13 | `gpu-port` plus point-defect model, all GPU lines (main branch) |
| `JAKOB_BACKUP` | keep (Jakob's data line) | — | +0 / +0 | 2026-11-09 | Jakob's defect backup, pin on first time slice |
| `luis_playground` | keep (your MaxEnt playground) | — | +7 / +7 | 2026-11-13 | your MaxEnt playground (Gaussian bootstrap) |
| `master` | keep (upstream mirror) | — | +0 / +0 | 2026-10-30 | upstream ALF 521ddfca |
| `gpu-delayed-submatrix-defect-production-mpi-amd` | rename | feat/amd-hip (+1 not in unified yet) | +1 / +1 | 2026-11-13 | `GPU=amd` HIP line |
| `gpu-port-apple-metal` | rename | exp/apple-metal (+1 not in unified yet) | +1 / +1 | 2026-11-13 | Apple Metal float-float ZGEMM prototype |
| `intel-batched` | rename | feat/intel-batched (+8 not in unified yet) | +11 / +8 | 2026-11-13 | Intel bottleneck study, batched rewrite begun |
| `production-mpi-intel-large-prompt-team-combined` | rename | feat/intel-devqr-team (+6 not in unified yet) | +6 / +6 | 2026-11-13 | Intel device QR, async, explicit-shape kernels |
| `s-wave-Pairing_correlations` | rename | feat/s-wave-pairing (+1 not in unified yet) | +1 / +1 | 2026-09-16 | on-site pairing correlation (2026-09) |
| `speedup/nvidia` | rename | feat/speedup-cuda-lasync (+4 not in unified yet) | +4 / +4 | 2026-11-13 | same three changes, separate commits |
| `gpu-delayed-submatrix-defect-production-mpi-pp` | tag as prod | prod/2026-11-10-a30-L24 | +0 / +0 | 2026-11-10 | parallel tempering on the GPU path |
| `gpu-delayed-submatrix-defect-production-mpi-pp-helma` | tag as prod | prod/2026-11-11-helma-L32 | +0 / +0 | 2026-11-11 | Helma H100 PP production chains |
| `production-mpi-intel-large-prompt` | tag as prod | prod/2026-11-12-smuc-large-prompt (and tag smuc2-prod-2026-11) | +0 / +0 | 2026-11-12 | Intel SuperMUC large-prompt production |
| `production-mpi-intel-unstable-guard` | tag as prod | no tag yet, branch still moves (tag smuc2-prod-2026-11b marks the submitted state) (+2 not in unified yet) | +2 / +2 | 2026-11-14 | `ALF_UNSTABLE_THRESHOLD`, SuperMUC scripts |
| `feature/delayed-submatrix` | tag as archive and delete | archive/feature/delayed-submatrix (check the 265 unmerged commits first) | +265 / +265 | 2026-11-01 | early delayed-submatrix work (SuperMUC tests) |
| `feature/delayed-submatrix-clean` | tag as archive and delete | archive/feature/delayed-submatrix-clean (check the 132 unmerged commits first) | +132 / +132 | 2026-11-05 | same, large text files removed from history |
| `fix/amd-checkerboard` | tag as archive and delete | archive/fix/amd-checkerboard | +0 / +0 | 2026-11-13 | AMD checkerboard hang fix (`!GCC$ novector`) |
| `gpu-delayed-submatrix-defect` | tag as archive and delete | archive/gpu-delayed-submatrix-defect | +0 / +0 | 2026-11-08 | first defect GPU port, delayed/submatrix |
| `gpu-delayed-submatrix-defect-fullgpu` | tag as archive and delete | archive/gpu-delayed-submatrix-defect-fullgpu | +0 / +0 | 2026-11-07 | defect port with full GPU path |
| `gpu-delayed-submatrix-defect-production` | tag as archive and delete | archive/gpu-delayed-submatrix-defect-production | +0 / +0 | 2026-11-09 | production line, JAKOB_BACKUP merged |
| `gpu-delayed-submatrix-defect-production-mpi` | tag as archive and delete | archive/gpu-delayed-submatrix-defect-production-mpi | +0 / +0 | 2026-11-09 | MPI with GPU path; base of `gpu-port-unified` |
| `gpu-delayed-submatrix-defect-production-mpi-pp-opt` | tag as archive and delete | archive/gpu-delayed-submatrix-defect-production-mpi-pp-opt | +0 / +0 | 2026-11-10 | PP profiling, `acc0` reduction (worktree in `/private/var/...wt_pp`) |
| `gpu-delayed-submatrix-defect-production-old-Ntau` | tag as archive and delete | archive/gpu-delayed-submatrix-defect-production-old-Ntau (check the 4 unmerged commits first) | +4 / +4 | 2026-11-09 | defect variant with old tau labelling of `SpinZ_localtau` |
| `gpu-port-207-2` | tag as archive and delete | archive/tickets/gpu-port-207-2 | +0 / +0 | 2026-11-01 | ticket 207-2: gate fix |
| `gpu-port-207-3` | tag as archive and delete | archive/tickets/gpu-port-207-3 | +0 / +0 | 2026-11-02 | ticket 207-3: device correct-side test |
| `gpu-port-207-4` | tag as archive and delete | archive/tickets/gpu-port-207-4 (check the 1 unmerged commits first) | +1 / +1 | 2026-11-02 | 207-4: beta=10 NVIDIA calibration scripts |
| `gpu-port-207-final` | tag as archive and delete | archive/tickets/gpu-port-207-final | +0 / +0 | 2026-11-02 | ticket 207 closing sections |
| `gpu-port-207-harness` | tag as archive and delete | archive/tickets/gpu-port-207-harness (check the 4 unmerged commits first) | +4 / +4 | 2026-11-02 | 207-5: harness QREF gates |
| `gpu-port-207-r` | tag as archive and delete | archive/tickets/gpu-port-207-r | +0 / +0 | 2026-11-01 | merge of 207-2 into 207-r |
| `gpu-port-207-tau` | tag as archive and delete | archive/tickets/gpu-port-207-tau | +0 / +0 | 2026-10-31 | ticket 207: quad-precision tau reference |
| `gpu-port-207-wip` | tag as archive and delete | archive/tickets/gpu-port-207-wip (check the 5 unmerged commits first) | +5 / +5 | 2026-10-29 | 207 evidence record, scripts (on hold) |
| `gpu-port-208-1` | tag as archive and delete | archive/tickets/gpu-port-208-1 | +0 / +0 | 2026-11-02 | ticket 208-1 closed |
| `gpu-port-208-2` | tag as archive and delete | archive/tickets/gpu-port-208-2 | +0 / +0 | 2026-11-02 | ticket 208-2 closed |
| `gpu-port-208-3` | tag as archive and delete | archive/tickets/gpu-port-208-3 | +0 / +0 | 2026-11-03 | ticket 208-3 closed (old `gpu-port` tip 43045813) |
| `gpu-port-212-hardening` | tag as archive and delete | archive/tickets/gpu-port-212-hardening (check the 2 unmerged commits first) | +2 / +2 | 2026-10-29 | 212: pre-registration only |
| `gpu-port-212-hardening-b` | tag as archive and delete | archive/tickets/gpu-port-212-hardening-b | +0 / +0 | 2026-10-30 | 212: T2 pass, closed |
| `gpu-port-302a-onemkl` | tag as archive and delete | archive/tickets/gpu-port-302a-onemkl | +0 / +0 | 2026-10-30 | 302a oneMKL, T2 pass |
| `gpu-port-302b` | tag as archive and delete | archive/tickets/gpu-port-302b | +0 / +0 | 2026-10-31 | 302b ticket, pre-registration |
| `gpu-port-302b-docs` | tag as archive and delete | archive/tickets/gpu-port-302b-docs (check the 2 unmerged commits first) | +2 / +2 | 2026-10-31 | 302b review and T1/T2 record |
| `gpu-port-302b-impl` | tag as archive and delete | archive/tickets/gpu-port-302b-impl | +0 / +0 | 2026-10-31 | 302b candidate 1 |
| `gpu-port-302b-impl2` | tag as archive and delete | archive/tickets/gpu-port-302b-impl2 | +0 / +0 | 2026-11-01 | 207 split into 207-1 ... 207-7, PLAN |
| `gpu-port-302c` | tag as archive and delete | archive/tickets/gpu-port-302c | +0 / +0 | 2026-11-01 | 302c-1: device tau pass, `GPU=intel` |
| `gpu-port-302c-2` | tag as archive and delete | archive/tickets/gpu-port-302c-2 | +0 / +0 | 2026-11-02 | 302c-2: B580 control ladders |
| `gpu-port-delayed-submatrix` | tag as archive and delete | archive/tickets/gpu-port-delayed-submatrix | +0 / +0 | 2026-11-06 | NOT_IMPLEMENTED.md, PARKING.md |
| `gpu-port-harness-302b` | tag as archive and delete | archive/tickets/gpu-port-harness-302b | +0 / +0 | 2026-10-31 | acceptance harness after review (c) |
| `gpu-port-harness-gpunone` | tag as archive and delete | archive/tickets/gpu-port-harness-gpunone | +0 / +0 | 2026-10-30 | harness exclusion list, `GPU=none` |
| `gpu-port-process-lanes` | tag as archive and delete | archive/tickets/gpu-port-process-lanes | +0 / +0 | 2026-10-29 | harness fixes, fourth review |
| `gpu-port-unified-maxent` | tag as archive and delete | archive/gpu-port-unified-maxent (check the 1 unmerged commits first) | +3 / +1 | 2026-11-13 | work branch: MaxEnt bootstrap merge, PP tool |
| `gpu-port-unified-tests` | tag as archive and delete | archive/gpu-port-unified-tests | +0 / +0 | 2026-11-13 | work branch: Run.tests, `run-local` |
| `intel-supermuc-ng2` | tag as archive and delete | archive/intel-supermuc-ng2 | +0 / +0 | 2026-11-10 | SuperMUC-NG Phase 2 session, benchmarks |
| `production-mpi-intel-large-prompt-team` | tag as archive and delete | archive/production-mpi-intel-large-prompt-team (check the 1 unmerged commits first) | +1 / +1 | 2026-11-13 | Intel team issue, staging-buffer check |
| `production-mpi-intel-large-prompt-team-async` | tag as archive and delete | archive/production-mpi-intel-large-prompt-team-async | +0 / +0 | 2026-11-12 | Intel async block path (phase 3b) |
| `production-mpi-intel-large-prompt-team-sycqr` | tag as archive and delete | archive/production-mpi-intel-large-prompt-team-sycqr | +0 / +0 | 2026-11-12 | Intel SYCL device Householder QR |
| `production-mpi-intel-unstable-guard-asyncblk-test` | tag as archive and delete | archive/production-mpi-intel-unstable-guard-asyncblk-test (check the 2 unmerged commits first) | +2 / +2 | 2026-11-13 | merge of unstable-guard and team-combined at older tips |
| `speedup/unified` | tag as archive and delete | archive/speedup/unified (check the 4 unmerged commits first) | +4 / +4 | 2026-11-13 | LU-free det phase (S1), `la_sync` removal on CUDA (S2), nsys script (S4) |

Not in the table: the six `defect/*` refs are remote-tracking refs of the sibling repo `../ALF_Defect-2`, not branches of `origin`; three are merged (`master`, `Point_Pin`, `Point_Pin_updated`), three hold unmerged commits (`localtau`, `localtau_Obs`, `Predefined_Hoppings_set_OPT-pinning`) and would get `archive/` tags in that repo. Existing tags `archive/gpu-delayed-submatrix-defect-fullgpu-old-Ntau`, `archive/gpu-delayed-submatrix-defect-production-old-Ntau-mpi`, `pre-lanes-2026-10-29`, `smuc2-prod-2026-11`, `smuc2-prod-2026-11b` stay; the two `smuc2-prod-*` tags could be copied to `prod/` names.

Counts of the action column: 1 trunk, 1 variant, 3 keep, 6 rename, 4 tag as prod, 41 tag as archive and delete.

## Carried out

### 2026-11-14: tags

- 41 `archive/` tags (24 of them `archive/tickets/`) and 3 `prod/` tags created at the branch tips and pushed.
  Check: `git ls-remote origin`, tag hash equal to branch hash for all 44.
- The `prod/` tags mark the tip of the production branch on the date in the name. Which commit a given binary was
  built from was not checked; where it is recorded it is in the tag message (`smuc2-prod-2026-11`: binary from `5667acda`).
- `smuc2-prod-2026-11` existed only locally; pushed.
- Three branches existed only in the old checkout `CLAUDE/ALF` and on no remote. Tagged and pushed from there:
  `archive/tickets/gpu-port-delayed-submatrix-rank-nd` (24c85e9e, WIP 401-A, general vertices on the device, **not in the trunk**),
  `archive/tickets/harness-e2e-local` (f02ee2d5), `archive/tickets/harness-e2e2-local` (f5f7ffee).

### 2026-11-14: deleted and renamed (owner's yes on the shown lists)

- Before deleting, each tag was compared again with the branch tip on `origin` (40 of 40 equal).
- Deleted on `origin`: the 32 branches without a worktree, then the 8 branches whose worktrees were removed first
  (`fix/amd-checkerboard`, `gpu-port-unified-maxent`, `gpu-port-unified-tests`,
  `production-mpi-intel-large-prompt-team`, `-team-async`, `-team-sycqr`,
  `production-mpi-intel-unstable-guard-asyncblk-test`, `speedup/unified`). 14 local copies deleted as well.
- Removed worktrees: `ALF_wt_fix_amd_cb`, `ALF_wt_unified_maxent`, `ALF_wt_unified_tests`, `ALF_wt_large_prompt_team`,
  `ALF_wt_lp_team_async`, `ALF_wt_lp_team_sycqr`, `ALF_wt_asyncblk_test`, `ALF_wt_speedup_unified` (all clean and pushed).
- Renamed: `s-wave-Pairing_correlations` is now `feat/s-wave-pairing` (ccab1c6f).
- Restore any deleted branch with `git switch -c <name> archive/<name>` (or `archive/tickets/<name>`).

### Branches on `origin` now (16)

| Branch | Role | Still to do |
|---|---|---|
| `gpu-port` | trunk | - |
| `gpu-port-unified` | variant (trunk plus point defect) | - |
| `master`, `JAKOB_BACKUP`, `luis_playground` | kept | - |
| `feat/s-wave-pairing` | feature | - |
| `gpu-delayed-submatrix-defect-production-mpi-amd` | AMD HIP line (worktree `ALF_wt_amd`) | rename to `feat/amd-hip` |
| `gpu-port-apple-metal` | prototype (worktree `ALF_wt_apple_metal`) | rename to `exp/apple-metal` |
| `intel-batched` | Intel speed work (worktree, running session) | rename to `feat/intel-batched` |
| `production-mpi-intel-large-prompt-team-combined` | Intel team (worktree, running session) | rename to `feat/intel-devqr-team` |
| `speedup/nvidia` | CUDA speedups (worktree `ALF_wt_speedup_nv`) | rename to `feat/speedup-cuda-lasync` |
| `gpu-delayed-submatrix-defect-production-mpi-pp` | production base, tagged `prod/2026-11-10-a30-L24` | delete when the A30 run is finished |
| `gpu-delayed-submatrix-defect-production-mpi-pp-helma` | Helma production (worktree `ALF_GPU_SUBMATRIX`), tagged `prod/2026-11-11-helma-L32` | delete when the Helma runs are finished |
| `production-mpi-intel-large-prompt` | SuperMUC production, tagged `prod/2026-11-12-smuc-large-prompt` | delete when the runs are finished |
| `production-mpi-intel-unstable-guard` | SuperMUC production, still moving (worktree `ALF_wt_large_prompt`) | `prod/` tag when it stops moving, then delete |
| `gpu-delayed-submatrix-defect-production-mpi-pp-opt` | archived as a tag; its worktree is in a temporary folder | remove the worktree entry, then delete |

Also left: the local branches and 25 worktrees of the old checkout `CLAUDE/ALF`, and the `defect/*` refs of `../ALF_Defect-2`.
