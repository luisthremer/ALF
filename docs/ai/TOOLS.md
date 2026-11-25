# Tools for understanding ALF without reading all of it

Checked 2026-11-25. In order of use:

1. **This folder.** `AGENTS.md` (what goes where), `docs/ai/CODEMAP.md` (module map, call graph,
   data ownership), `docs/ai/codemap/build.md` (build, tests, CI). Start here; it costs a few
   thousand tokens instead of the ~20 000 lines of `Prog/`.
2. **`gh` CLI** for pull requests and CI: `gh pr view 657 -R ALF-QMC/ALF`,
   `gh run list -R luisthremer/ALF --branch <branch>`, `gh run view <id> --log-failed`.
3. **`grep` on routine names.** Fortran names are unique enough; `grep -n "subroutine Upgrade2" Prog/*.F90`.
4. **Fortran language server `fortls`** (pip package `fortls`): go-to-definition, references, hover
   on Fortran 2008 modules and submodules. Usable through an LSP-aware agent tool:
   - **Serena** (MCP server, in the official Claude plugin marketplace as `serena`): symbol-level
     search and editing on top of language servers; its language list includes Fortran (`fortls`).
     Not yet tried on ALF.
   - or a small Claude Code plugin with an `.lsp.json` that starts `fortls` (there is no official
     Fortran LSP plugin; the official ones cover C/C++, Go, Rust, Python, …).
5. **Understand-Anything** (Claude Code plugin): `/understand` builds an interactive knowledge
   graph and dashboard of the code. Good for a human overview; heavy to build, so rebuild only after
   a phase.

Not useful here: `ctags` on macOS (BSD ctags has no Fortran parser); GitHub code search (does not
index the fork's branches).
