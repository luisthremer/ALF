## Now
- P2.4 step 1f (S1) is merged into the GPU line. Step 2 in work: gpu-nvidia, gpu-amd, gpu-intel have accepted design notes and write code (common rule S2-1: presence dispatch, one vendor module each).
- Also in parallel: b-batch-2 (CPU-line B bugs H5, M5, M6, MK), cpu-s1 (CG-ALF T5: tracked phase check, LU-free phase on the CPU line), cpu-queue (alf_exit, then Timers_mod), cgalf-research (H100 stabilisation scan).
- GPU line CI run 45890936664 on 0bc71e66 (first run of 1f together with Trace_updates_mod).
- Live progress page for you: [ALF GPU Phase 2](https://claude.ai/artifact/PmMk2VMb8iUzXtnEYUcHc3) (private; updated after every merge or push).

## Waiting for you
Answer with the buttons under "Needs you" on the [progress page](https://claude.ai/artifact/PmMk2VMb8iUzXtnEYUcHc3).
1. SuperMUC login for gpu-intel's device sessions (it asks on the progress page when ready); first cap 40 GPU-hours of the estimated 76.
("Needs you" on the page is live: sessions add a row when they need you.)

## Done since the last owner message
- Q13 decided: `Lowrank_cross_slice` stays off the CM lines, kept on branch `rd/lowrank-cross-slice`.
- Q15 answered yes: claude-bughunter for the Intel `Mz = T` investigation. Merged `alf_exit` (CPU 7466aa71, GPU 04c26bce).
- Merged P2.4 step 1f (S1: CGR phase from the tracked det(U), no LU) into the GPU line: [a3e1220f](https://github.com/luisthremer/ALF/commit/a3e1220f), CI run 45886167489 green (30/30). The LU removal is the single commit 789fe086 (`git revert 789fe086` gives the LU version exactly).
- Merged `Trace_updates_mod` (debug trace of every Metropolis decision, behind `ALF_TRACE_UPDATES`; default build unchanged) into the CPU line ([0440f859](https://github.com/luisthremer/ALF/commit/0440f859), CI run 45888075085 green) and the GPU line (0bc71e66).
- Q14 answered "Yes" on the progress page: session cgalf-research measures the stabilisation candidates on the Helma H100.
- Merged the G2/G3 fixes into the GPU line: [a576c1b6](https://github.com/luisthremer/ALF/commit/a576c1b6), CI run 45881567484 green.
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
1. Start the three vendor sessions of step 2 (NVIDIA on Helma, AMD on amd, Intel on SuperMUC), one at a time.
2. CPU-line queue: merge `alf_exit` when its rebased canvas is green (cpu-queue session).
3. CG-ALF: H100 measurement of the stabilisation candidates (cgalf-research, Q14 yes).

## Open branches and running jobs
- `cm/feat/alf-exit` (CPU line, CI run 45890618713). Step-2 branches `cm/*/gpu-{nvidia,amd,intel}`, `cm/*/b-batch-2`, `rd/cgalf-t5` + `cm/*/cpu-s1` (being created).
- GPU line CI run 45890936664 on 0bc71e66 (1f + Trace_updates_mod).
- Cluster jobs: cgalf-research on Helma (cap 4 GPU-h), gpu-nvidia on Helma (cap 10 GPU-h), gpu-amd on amd, gpu-intel on SuperMUC (needs your login). Plan: [P2_4_GPU_PLAN.md](https://github.com/luisthremer/ALF/blob/Code-modularity-GPU/docs/ai/phase2/P2_4_GPU_PLAN.md)
