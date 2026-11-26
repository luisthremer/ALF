## Now
- G2/G3 fixes on `cm/feat/gpu-hopd-optin` (CI running): dense hopping product `Hd` opt-in (`ALF_GPU_HOPD=1`), so the device path has the CPU's precision; the device precision check now honours `ALF_UNSTABLE_THRESHOLD` (G3, found when 1e met the parallel CPU-line merge). melmac: tests 76/76 (GPU=none), 119/119 (GPU=host, also bounds-checked with FP traps); the threshold run aborts as it should.
- P2.4 step 1f (speedup S1, no LU) on `cm/dev/gpu-s1`: coding agent working. Four commits; the fourth (LU removal) alone, so one `git revert` restores the LU.

## Waiting for you
1. Q13: `Lowrank_cross_slice` is sound but not faster on the CPU: drop it (recommended) or keep as opt-in?
2. Q14: measure the stabilisation candidates (CholeskyQR2, adaptive Nwrap) on the H100 (a few GPU-hours)? yes/no

## Done since the last owner message
- Merged step 1e (device observables incl. local h0 on the device, device tau pass): [b75f9d91](https://github.com/luisthremer/ALF/commit/b75f9d91), CI run 45876811655 green; all 9 step commits built and tested on melmac (Linux); device observables need no G download (counter 0); device vs CPU observables within roundoff.
- Merged on the CPU line (cpu-queue session): `ALF_UNSTABLE_THRESHOLD` (configurable G precision abort threshold) and `Start_PP/prepare_sim.py` (general PARALLEL_PARAMS helper, [954ba129](https://github.com/luisthremer/ALF/commit/954ba129)); both also in the GPU line.
- G2 explained (with session g2-precision): only with Apple's Accelerate BLAS is the dense `Hd` product less precise (Linux reference BLAS: not worse); fix pending above.
- G2 cause found (session g2-precision): the device multiplies with one dense precomputed product of all checkerboard factors (`Hd` in `gpu_wrap_mod`, on for Ndim <= 800) instead of factor by factor; with it off, Precision Green equals the CPU's (ratio 1.00-1.01 over 12 seeds; was 1.4-2.5). Fix planned as a small step after 1e: `Hd` opt-in, default the per-factor path.
- CG-ALF T6 and T1 results merged into the plan (both lines, docs): cross-slice update sound but not faster on the CPU (1.015x at L = 16); stabilisation survey of 10 methods (CG_ALF_PLAN.md section 7).
- Merged step 1d (device stabilisation, determinant-phase fix for large N) into `Code-modularity-GPU`: [36d17988](https://github.com/luisthremer/ALF/commit/36d17988), CI run 45868504866 green. Local: tests 70/70 (GPU=none), 95/95 (GPU=host default, STABLOG, bounds-checked); 16/16 byte-identical; phase test at N = 3200 passes; device vs CPU on Hubbard 4x4: configurations byte-identical, energies within 1.2e-10 (Mz), 5e-14 (SU(2)). Fixed on the way: the NaN self-test fed a NaN to the host reference (trap in the CI debug job).
- P2.5 feature table, first draft merged into the GPU line (docs): [GPU_FEATURES.md](https://github.com/luisthremer/ALF/blob/Code-modularity-GPU/docs/ai/phase2/GPU_FEATURES.md) (vendor columns copied from gpu-port docs, not re-checked; GPU=host column from the code).
- All questions answered; decisions in the plans ([255251ff](https://github.com/luisthremer/ALF/commit/255251ff)): D4 (b) with a CPU-line queue, S1 as step 1f (same method on every vendor, LU removal as one revertible commit), Intel switches default with info report, P2.6 starts with delayed/submatrix for every model, CG-ALF tasks T5-T11.
- Merged step 1c ([0c222961](https://github.com/luisthremer/ALF/commit/0c222961)) and 1b ([eb5b1022](https://github.com/luisthremer/ALF/commit/eb5b1022)), CI green; G1 and two more problems found and fixed on the way.
- Earlier today: Phase-6 sync, B batch 1, step 1a, CG-ALF ideas, upstream issues, branch clean-up (history: `git log manager_tools/STATUS.md`).

## Next 3 steps
1. Merge the G2/G3 fixes; review and merge 1f (S1).
2. CPU-line queue: next item `alf_exit` / `Trace_updates_mod` (cpu-queue session).
3. Merge what the parallel sessions deliver (CPU-line queue, CG-ALF research results); then step 2 vendor backends (clusters).

## Open branches and running jobs
- `cm/feat/gpu-hopd-optin` + `cm/dev/gpu-hopd-optin` (G2/G3, CI); `cm/dev/gpu-s1` (step 1f, local, coding agent).
- Parallel sessions reporting to the manager: cpu-queue (`ALF_UNSTABLE_THRESHOLD`, `prepare_sim.py`, CPU line), cgalf-research (T6, T1 done; waits for Q13/Q14), g2-precision (done). Builds with `make -j2`, `ctest -j2` while many sessions run.
- Side session "Interview owner on QUESTIONS.md decisions": finished; can be archived on your yes.
- Cluster jobs: none. Plan: [P2_4_GPU_PLAN.md](https://github.com/luisthremer/ALF/blob/Code-modularity-GPU/docs/ai/phase2/P2_4_GPU_PLAN.md)
