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
DONE (2026-11-03, decision 42): 208-3 closed; merged into gpu-port, pushed to fork. Evidence: work/t208-3f/logs/SUMMARY.txt on nividia and intel
(one launch, docs/gpu/scripts/t208_3/launch_final.sh; the first chain in work/t208-3 had the same numbers but a self-test failure, fixed by
requiring prec_report for the device tau observables).

## Results
1 Builds: nividia 5, intel 5: PASS.
2 Self-tests (=4): PASS on both, def and LOG.
3 Same-seed hub6x6_cb_fd max dev (incl. the 5 tau-resolved datasets): nividia 1.70e-10 / 2.54e-10, intel 2.32e-10 / 2.90e-10 (def / LOG): PASS,
  identical to 207, 208-1, 208-2.
4 Observables hub12x12_cb vs CPU ref: nividia sum t^2 = 9.2 (p 0.60) def and LOG; intel 141.7 (the accepted 207 Szz(pi,pi) group): same as 207.
  Wall: nividia def 625 s (208-2: 703), LOG 707 s (841); intel 412 / 402 s (419 / 409).
5 GPU=none ctest (GNU, nividia): candidate 53/53, upstream master 39/39: PASS.
Transfers (obs_def): INTERIM 207 (tau pass) H->D 0 and D->H 0 (before: GT0/G0T/GTT per slice and G00 per stabilization point);
INTERIM 205 G copies H->D 0, D->H 1200 (398 MB per run, the host copy for the tau pass; was 31.85 GB in 208-1).

## Next step
Phase 2 remainder (PLAN.md): 209 validation suite (209-1 test_vs_ed on the GPU, 209-2 same-compiler trajectory test, 209-3 one runnable suite),
210 BENCHMARKS.md, 211 CPU regression gate, 213 STABLOG same-seed excess. Owner decides the order and the process (decision 42 covers 208 only).

## Running
Nothing.
