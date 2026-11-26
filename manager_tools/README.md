# manager_tools — the owner's folder

| File | What |
|---|---|
| `STATUS.md` | The only status page of GPU Development Phase 2 (format and rules: `AGENTS.md`, "Status and ideas"). |
| `IDEAS.md` | The owner's idea inbox. Claude never edits or deletes the owner's lines; it only appends "-> taken: <where> (<date>)". |
| `netgraph.sh`, `netgraph.py`, `page.html` | GitHub-style network graph of this repo: `./manager_tools/netgraph.sh` opens `out/index.html` and rebuilds it on every ref change. |
| `netgraph.json` | Which branches the graph shows (`only`; edit it when the active line changes), remotes hidden, archive tags that name the rows of merged branches. |

Lives on branch `Code-modularity` of the fork (public). Status and idea updates are committed directly on
`Code-modularity` with `[skip ci]`; they are not code. Until 2026-11-26 STATUS.md and IDEAS.md were in the
private repo `luisthremer/alf-ops`.
