## Now
Fakher's newest Phase-6 commits (6G.E1-E7 linalg wrappers, STAB1 projector fix) merged on `cm/feat/phase6-sync` (57/57 tests, 16/16 byte-identical) and B batch 1: both in CI, merged into the CPU line when green. P2.4 plan updated with your D1-D3; step 1a is redone for `linalg_mod` after the sync.

## Waiting for you
All questions in detail: [QUESTIONS.md](QUESTIONS.md).
1. Q1/D4: gpu-port's unguarded CPU-path changes later as separate CPU-line features? yes/no
2. Q3: speedup S1 into the GPU line? (memo for your professor in QUESTIONS.md)
3. Q4a-d older GPU-port questions; Q5-Q10 things only you can do.

## Done since the last owner message
- Upstream issues posted: ALF-QMC/ALF#662 (H7), #663 (H6). P2.4 plan: D1 linalg_mod, D2 GPU=host, D3 Hubbard only.
- B-bug batch 1 built and checked (8 bugs, regression tests); step 1b (GPU state) designed.
- `manager_tools/` now on `Code-modularity`, `-CPU` and `-GPU`, kept equal by `manager_tools/sync.sh`.
- P2.4 plan: gpu-port mapped (12 436 lines of common OpenMP offload code, vendor shims, hooks) and re-homed onto the refactored code; steps 1a-1e before vendor code.
- Status and ideas moved from alf-ops into `manager_tools/` on `Code-modularity` (public, your OK); understand-anything removed (traces in the Trash); claude-bughunter off until P2.8.
- The CPU determinant-phase fix is already offered upstream: PR ALF-QMC/ALF#660 (open).
- Merged into `Code-modularity`: nan-guards [f8dbc251](https://github.com/luisthremer/ALF/commit/f8dbc251) (CI green), lowrank-fastest [eaf94e27](https://github.com/luisthremer/ALF/commit/eaf94e27) (CI green; local: 57/57 tests, 16/16 byte-identical), cg-alf-plan [c2efcce5](https://github.com/luisthremer/ALF/commit/c2efcce5) (docs), p23-split [40b05e46](https://github.com/luisthremer/ALF/commit/40b05e46) (docs). Archive tags pushed.
- P2.3: `Code-modularity-CPU` / `-GPU` created; `Code-modularity` stays at the split point.
- Branch clean-up: 20 merged/obsolete branches deleted, `luis_playground` renamed to `PP-Upgrades_and_Max-Sac-Bootstrapping` (tag `stp/max-sac-bootstrapping`). New rule merged ([b64e88b2](https://github.com/luisthremer/ALF/commit/b64e88b2)): merged and tagged `cm/*` branches are deleted right after the merge. `JAKOB_BACKUP` and the three `production-mpi-intel-*` branches tagged `stp/…` and deleted.
- Inbox: two ideas now CG-ALF tasks T3 (fast `Op_V_is_symmetric`) and T4 (production guards from the Intel code); inbox is checked at least every two hours of work ([fc8561cf](https://github.com/luisthremer/ALF/commit/fc8561cf)).

## Next 3 steps
1. P2.4 plan for the GPU scaffold (docs), then your look at it.
2. Scaffold dev branch `cm/dev/gpu-common` (GPU=none bitwise identical), then per-vendor branches.
3. CG-ALF tasks T3 (fast `Op_V_is_symmetric`) and T4 (production guards).

## Open branches and running jobs
- `cm/feat/b-batch-1` and `cm/feat/phase6-sync` (CI running, CPU line); `cm/feat/gpu-build` (to be redone for D1).
- Cluster jobs: none. Plan: [P2_4_GPU_PLAN.md](https://github.com/luisthremer/ALF/blob/Code-modularity-GPU/docs/ai/phase2/P2_4_GPU_PLAN.md)
