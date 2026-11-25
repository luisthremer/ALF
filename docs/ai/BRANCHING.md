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
Add a row to the table below. Delete the branches only after the owner says yes.

Several dev branches per feature are fine (`cm/dev/<feature>-<part>`); larger efforts may use an
integration branch (`cm/dev/<feature>-integrate`) before the single canvas.

## Archive

| Feature | Target | Merge | Dev tag | Canvas tag | What it did |
|---|---|---|---|---|---|
| ai-docs | Code-modularity | `c771d98a` (2026-11-25) | `archive/cm/ai-docs/dev` | `archive/cm/ai-docs/feat` | AGENTS.md, CLAUDE.md, docs/ai (code map, branching rule, tools, phase-2 brief and roadmap), KNOWN_ISSUES.md |
| lowrank | Code-modularity | `2bb0258a` (2026-11-25) | `archive/cm/lowrank/dev` | `archive/cm/lowrank/feat` | native delayed/submatrix update (`lowrank_update_mod`, `Sequential_scheme`, `N_delay`), tests, documentation |
