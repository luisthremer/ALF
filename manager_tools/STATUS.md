## Now
P2.3 done: `Code-modularity-CPU` and `Code-modularity-GPU` created at [40b05e46](https://github.com/luisthremer/ALF/commit/40b05e46). Session stopped here (stop point of the handoff); next is P2.4.

## Waiting for you
1. Upstream issues for H7 (Bilayer_honeycomb energy) and H6 (Tempering hang, odd replicas): texts drafted in the manager chat, post them? yes/no
2. Take speedup S1 (determinant phase without LU on CUDA, -12.6 % per sweep at L = 32, configurations byte-identical) into the GPU line? yes/no — [CG_ALF_PLAN.md, E4](https://github.com/luisthremer/ALF/blob/Code-modularity/docs/ai/phase2/CG_ALF_PLAN.md)
3. Switch off research-superpowers and research-desk in claude.ai (account plugins; Claude cannot), and set the GitHub MCP token (line in the manager chat).

## Done since the last owner message
- Status and ideas moved from alf-ops into `manager_tools/` on `Code-modularity` (public, your OK); understand-anything removed (traces in the Trash); claude-bughunter off until P2.8.
- The CPU determinant-phase fix is already offered upstream: PR ALF-QMC/ALF#660 (open).
- Merged into `Code-modularity`: nan-guards [f8dbc251](https://github.com/luisthremer/ALF/commit/f8dbc251) (CI green), lowrank-fastest [eaf94e27](https://github.com/luisthremer/ALF/commit/eaf94e27) (CI green; local: 57/57 tests, 16/16 byte-identical), cg-alf-plan [c2efcce5](https://github.com/luisthremer/ALF/commit/c2efcce5) (docs), p23-split [40b05e46](https://github.com/luisthremer/ALF/commit/40b05e46) (docs). Archive tags pushed.
- P2.3: `Code-modularity-CPU` / `-GPU` created; `Code-modularity` stays at the split point.
- Branch clean-up: 20 merged/obsolete branches deleted, `luis_playground` renamed to `PP-Upgrades_and_Max-Sac-Bootstrapping` (tag `stp/max-sac-bootstrapping`). New rule merged ([b64e88b2](https://github.com/luisthremer/ALF/commit/b64e88b2)): merged and tagged `cm/*` branches are deleted right after the merge. `JAKOB_BACKUP` and the three `production-mpi-intel-*` branches tagged `stp/…` and deleted.
- Inbox: two ideas now CG-ALF tasks T3 (fast `Op_V_is_symmetric`) and T4 (production guards from the Intel code); inbox is checked at least every two hours of work ([fc8561cf](https://github.com/luisthremer/ALF/commit/fc8561cf)).

## Next 3 steps
1. P2.4 step 1: plan the common GPU scaffold (`cm/dev/gpu-common` from `Code-modularity-GPU`).
2. Batch fix of the B bugs in KNOWN_ISSUES.md (after the phase merge).
3. CG-ALF tasks T3 (fast `Op_V_is_symmetric`) and T4 (production guards) after the P2.4 scaffold.

## Open branches and running jobs
- No unmerged feature branches. Baseline CI on `Code-modularity-CPU` fc8561cf ([run 45797219302](https://github.com/luisthremer/ALF/actions/runs/45797219302)).
- Cluster jobs: none. Plan: [ROADMAP.md](https://github.com/luisthremer/ALF/blob/Code-modularity-CPU/docs/ai/phase2/ROADMAP.md)
