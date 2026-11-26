## Now
- P2.4 step 1e (`gpu-obs-tau`: observables and tau pass on the device, local h0 observable on the device instead of the silent host fallback) on `cm/dev/gpu-obs-tau`: coding agent working (observable part done, tau pass in work; builds with `-j2`).
- P2.4 step 1f (speedup S1, no LU): design done; starts on top of the 1e canvas. Commit 4 of 4 is the LU removal alone, so one `git revert` restores the LU.

## Waiting for you
Nothing.

## Done since the last owner message
- Merged step 1d (device stabilisation, determinant-phase fix for large N) into `Code-modularity-GPU`: [36d17988](https://github.com/luisthremer/ALF/commit/36d17988), CI run 45868504866 green. Local: tests 70/70 (GPU=none), 95/95 (GPU=host default, STABLOG, bounds-checked); 16/16 byte-identical; phase test at N = 3200 passes; device vs CPU on Hubbard 4x4: configurations byte-identical, energies within 1.2e-10 (Mz), 5e-14 (SU(2)). Fixed on the way: the NaN self-test fed a NaN to the host reference (trap in the CI debug job).
- P2.5 feature table, first draft merged into the GPU line (docs): [GPU_FEATURES.md](https://github.com/luisthremer/ALF/blob/Code-modularity-GPU/docs/ai/phase2/GPU_FEATURES.md) (vendor columns copied from gpu-port docs, not re-checked; GPU=host column from the code).
- All questions answered; decisions in the plans ([255251ff](https://github.com/luisthremer/ALF/commit/255251ff)): D4 (b) with a CPU-line queue, S1 as step 1f (same method on every vendor, LU removal as one revertible commit), Intel switches default with info report, P2.6 starts with delayed/submatrix for every model, CG-ALF tasks T5-T11.
- Merged step 1c ([0c222961](https://github.com/luisthremer/ALF/commit/0c222961)) and 1b ([eb5b1022](https://github.com/luisthremer/ALF/commit/eb5b1022)), CI green; G1 and two more problems found and fixed on the way.
- Earlier today: Phase-6 sync, B batch 1, step 1a, CG-ALF ideas, upstream issues, branch clean-up (history: `git log manager_tools/STATUS.md`).

## Next 3 steps
1. Review, canvas, CI and merge 1e.
2. 1f: S1 in four commits (tracking, check mode, tests, LU removal).
3. Merge what the parallel sessions deliver (CPU-line queue, CG-ALF research results); then step 2 vendor backends (clusters).

## Open branches and running jobs
- `cm/dev/gpu-obs-tau` (step 1e, local, coding agent).
- Parallel sessions reporting to the manager: cpu-queue (`ALF_UNSTABLE_THRESHOLD`, `prepare_sim.py`, CPU line), cgalf-research (T6 cross-slice, T1 stabilisation survey, `rd/` branches), g2-precision (G2, `rd/g2-precision`). Builds with `make -j2`, `ctest -j2` while many sessions run.
- Side session "Interview owner on QUESTIONS.md decisions": finished; can be archived on your yes.
- Cluster jobs: none. Plan: [P2_4_GPU_PLAN.md](https://github.com/luisthremer/ALF/blob/Code-modularity-GPU/docs/ai/phase2/P2_4_GPU_PLAN.md)
