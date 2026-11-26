# manager_tools — the owner's folder

| File | What |
|---|---|
| `STATUS.md` | The only status page of GPU Development Phase 2 (format and rules: `AGENTS.md`, "Status and ideas"). |
| `QUESTIONS.md` | All open questions for the owner with background, options and recommendation. |
| `IDEAS.md` | The owner's idea inbox. Claude never edits or deletes the owner's lines; it only appends "-> taken: <where> (<date>)". |
| `netgraph.sh`, `netgraph.py`, `page.html` | GitHub-style network graph of this repo: `./manager_tools/netgraph.sh` opens `out/index.html` and rebuilds it on every ref change. |
| `sync.sh` | Copies this folder from the current branch to `Code-modularity`, `Code-modularity-CPU` and `Code-modularity-GPU` (one `[skip ci]` commit each, only `manager_tools/` changes) and pushes them to the fork. Run it after every change here. |
| `netgraph.json` | Which branches the graph shows (`only`; edit it when the active line changes), remotes hidden, archive tags that name the rows of merged branches. |

Lives on `Code-modularity`, `Code-modularity-CPU` and `Code-modularity-GPU` of the fork (public), the same on all
three (owner, 2026-11-26). Status and idea updates are committed directly with `[skip ci]` (they are not code), on
any of the three, then `./manager_tools/sync.sh` copies them to the other two. Until 2026-11-26 STATUS.md and IDEAS.md were in the
private repo `luisthremer/alf-ops`.
