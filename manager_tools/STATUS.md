## Now
- P2.4 step 1d (`gpu-stab`: UDV stabilisation on the device, determinant-phase fix for large N; Hubbard only) on `cm/feat/gpu-stab`, 9 commits. The first CI run failed in one debug test (the NaN self-test fed a NaN to the host reference, which traps); fixed, CI rerun 45868504866 running (27/30 jobs done, none failed). Local: tests 70/70 (GPU=none), 95/95 (GPU=host default, STABLOG and bounds-checked); GPU=none 16/16 byte-identical; phase test at N = 3200 (old code: NaN) passes; device vs CPU path on Hubbard 4x4: configurations byte-identical in all 12 cases, energies within 1.2e-10 (Mz) and 5e-14 (SU(2)).
- P2.4 step 1e (`gpu-obs-tau`: observables and tau pass on the device, local h0 observable on the device instead of the silent host fallback) on `cm/dev/gpu-obs-tau`: one coding agent working.

## Waiting for you
1. G2 (class B, open): the device path's Precision Green is 1.4-6x the CPU's (already in 1c; observables unaffected). Likely causes: B and B^-1 not exact inverses on the device; complex hopping on the device where the CPU uses real exponentials (`RealExpOpT`, `ZDSLSYMM`). Investigate now (about 1 h, local) or before step 2?

## Done since the last owner message
- All questions answered ([QUESTIONS.md](QUESTIONS.md)); decisions in the plans on both lines ([255251ff](https://github.com/luisthremer/ALF/commit/255251ff)): D4 (b) with a CPU-line queue in ROADMAP (ALF_UNSTABLE_THRESHOLD first, prepare_sim.py as general PP helper); S1 as P2.4 step 1f (same method on every vendor, LU removal as ONE revertible commit); Intel switches default with info report; P2.6 starts with delayed/submatrix for every model; CG-ALF tasks T5-T11.
- Merged step 1c (device wraps, local updates, delayed/submatrix) into `Code-modularity-GPU`: [0c222961](https://github.com/luisthremer/ALF/commit/0c222961), CI run 45854274445 green; device vs CPU on Hubbard 4x4: configurations byte-identical, energies within 2e-11 (Mz) and 1e-14 (SU(2)). Fixed on the way: G1 (stale device G after the self-test), `Sequential_scheme = "auto"` would have stopped GPU runs, bounds errors in ported kernels (CI debug job).
- Merged step 1b (GPU state, device mirror of model and walker): [eb5b1022](https://github.com/luisthremer/ALF/commit/eb5b1022), CI run 45838961102 green.
- Earlier today: Phase-6 sync, B batch 1 (8 bugs), step 1a `GPU=none|host`, CG-ALF ideas I11-I13, upstream issues ALF-QMC/ALF#662 and #663, branch clean-up, `manager_tools/` on all three lines (history: `git log manager_tools/STATUS.md`).

## Next 3 steps
1. Merge 1d when CI is green; review and merge 1e.
2. 1f: speedup S1 (LU removal as one revertible commit).
3. Step 2: vendor backends inside linalg_mod (needs the clusters; Helma/SuperMUC logins then).

## Open branches and running jobs
- `cm/feat/gpu-stab` + `cm/dev/gpu-stab` (step 1d, CI run 45868504866); `cm/dev/gpu-obs-tau` (step 1e, local, not pushed yet).
- Side session "Interview owner on QUESTIONS.md decisions": finished; can be archived on your yes.
- Cluster jobs: none. Plan: [P2_4_GPU_PLAN.md](https://github.com/luisthremer/ALF/blob/Code-modularity-GPU/docs/ai/phase2/P2_4_GPU_PLAN.md)
