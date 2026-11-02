# 208-3 state (single source of truth; read only this file + CLAUDE.md)

Process: decision 42. Branch `gpu-port-208-3` (from gpu-port d502b762 = 208-2 closed). Push only to `fork`.

## Goal (PLAN.md 208-3)
Tau-resolved observables (ham%OBSERT) on the device; remove 207's per-slice D->H of GT0, G0T, GTT (3 * 16 N^2 B per flavour and slice) and the
D->H of G00 after each CGR2_2.

## Design / code (Prog/gpu_obs_mod.F90, gpu_tau_mod.F90)
gpu_obs_tau_measure(nt, GT0_d, G0T_d, G00_d, GTT_d, PHASE, w): per observable a gather per output element (imj, no, no2) at slice nt into the device
accumulator tacc(N, Nt, Norb, Norb, nobs) (+ tacc0 for Obs_Latt0 of Den and SpinMz; added at every nt as the host routines do), formulas of
Predefined_Obs_tau_{Green,SpinSUN,SpinMz,Den}_measure (Green: sum_nf GT0 / N_FL); counters Obs_tau%N / Ave_sign once per pass (nt = 0).
gpu_obs_flush (before Pr_obs) adds tacc/tacc0 into Obs_tau. gpu_tau_mod: tobs_dev (set in tau_begin: device obs active and no TAUDIAG/QREF)
makes tau_obs call the device routine and skips tau_download and the G00 D->H. Same support gate as 208-2 (no Symm, no Projector, N_FL 1/2).
Lesson 208-2 kept: complex weight as two real scalars wr, wi inside the loops.

## Current step
Device tau observables written and built; first chain (work/t208-3) showed same-seed/obs numbers identical to 208-2 (nividia wall def 625 s,
LOG 696 s; ctest 53/53) but the self-test failed because tau_begin switched to device observables inside the self-test; fixed (tobs_dev needs
prec_report). One launch for everything left: docs/gpu/scripts/t208_3/launch_final.sh -> work/t208-3f on both machines, logs/SUMMARY.txt
(ends FINAL_DONE) = builds, T0, observables, transfers, ctest.

## Next step (owner says "continue")
Read logs/SUMMARY.txt on nividia and intel; if all pass (self-tests PASSED, same-seed <= 1e-9, sum t^2 9.2 nividia / 141.7 intel, ctest 53/53):
write results here, mark 208-3 done in PLAN.md, merge gpu-port-208-3 into gpu-port (ff via determined-burnell worktree, then main checkout), push to fork.

## Running
final.sh on nividia and intel (started 2026-11-03).
