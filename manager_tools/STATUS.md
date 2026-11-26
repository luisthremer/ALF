## Now
P2.4 step 1c (`gpu-wrap-update`: wraps, local updates, delayed/submatrix on the device mirror; Hubbard only) on `cm/feat/gpu-wrap-update`, 9 commits; CI running. Local: every commit builds with GPU=none and GPU=host; tests 67/67 (none), 86/86 (host, default and STABLOG); GPU=none 16/16 byte-identical. Hubbard 4x4, device path vs CPU path in the same binary: configurations byte-identical; energies agree to 2e-11 absolute (1.5e-12 relative) with Mz, 1e-14 for SU(2).

## Waiting for you
All questions in detail: [QUESTIONS.md](QUESTIONS.md).
1. Q1/D4: gpu-port's unguarded CPU-path changes later as separate CPU-line features? yes/no
2. Q3: speedup S1 into the GPU line? (memo for your professor in QUESTIONS.md)
3. Q4a-d older GPU-port questions; Q5-Q10 things only you can do.

## Done since the last owner message
- Merged step 1b (GPU state, device mirror of model and walker) into `Code-modularity-GPU`: [eb5b1022](https://github.com/luisthremer/ALF/commit/eb5b1022), CI run 45838961102 green (30 jobs, incl. the new `gpu_host` jobs); 81/81 tests GPU=host, 67/67 GPU=none, 16/16 byte-identical. Branches archived and deleted.
- Merged: Phase-6 sync (CPU+GPU), B batch 1 (8 bugs, CPU, then into GPU), step 1a `GPU=none|host` (GPU); combination green in CI.
- Ideas I11-I13 (Cholesky again, time-first sweep, random space-time path) in the CG-ALF plan.
- Phase-6 sync merged (6G.E1-E7 linalg wrappers, STAB1 projector fix): 57/57 tests, 16/16 byte-identical, CI green.
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
1. 1c (wraps, updates, lowrank on the mirror), tested with GPU=host.
2. 1d stabilisation (with the GPU det-phase fix; S1 if you say yes), 1e observables/tau.
3. Step 2: vendor backends inside linalg_mod (needs the clusters).

## Open branches and running jobs
- `cm/feat/gpu-wrap-update` + `cm/dev/gpu-wrap-update` (step 1c, CI run 45851444872).
- Side session "Interview owner on QUESTIONS.md decisions" (waits for you).
- Cluster jobs: none. Plan: [P2_4_GPU_PLAN.md](https://github.com/luisthremer/ALF/blob/Code-modularity-GPU/docs/ai/phase2/P2_4_GPU_PLAN.md)
