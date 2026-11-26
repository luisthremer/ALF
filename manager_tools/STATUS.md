## Now
P2.4 step 1a `gpu-build` ready (canvas, 3 commits: GPU=none|host switch, `alf_la_mod` host backend, its tests): every commit builds, 69/69 tests, 16/16 scenarios byte-identical, CI green ([run 45804627437](https://github.com/luisthremer/ALF/actions/runs/45804627437)). Merge waits for your yes on question 4 (D1-D4).

## Waiting for you
All questions with background and recommendation: [QUESTIONS.md](QUESTIONS.md).
1. Q1 P2.4 decisions D1-D4 (merge of step 1a waits for this). yes/no
2. Q2 post the upstream issues H7 and H6? post/corrections
3. Q3 speedup S1 into the GPU line? (a) after step 1d / (b) CG-ALF / (c) no
4. Q4a-d older GPU-port questions: accept the proposed places? yes/no each
5. Q5-Q7 only you: delete alf-ops online, GitHub token, research plugins off in claude.ai

## Done since the last owner message
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
- `cm/feat/gpu-build` + `cm/dev/gpu-build` (step 1a, CI green, not merged: waits for question 4).
- Baseline CI of the new lines: green on `Code-modularity-CPU` fc8561cf ([run 45797219302](https://github.com/luisthremer/ALF/actions/runs/45797219302)).
- Cluster jobs: none. Plan: [P2_4_GPU_PLAN.md](https://github.com/luisthremer/ALF/blob/Code-modularity-GPU/docs/ai/phase2/P2_4_GPU_PLAN.md)
