# 208-1 state (single source of truth; read only this file + CLAUDE.md)

Process: decision 42 = decision 41 style (functional checks only; no agents, no skills, serial, one step per turn, long jobs detached,
machine rules of CLAUDE.md). Branch `gpu-port-208-1` (from gpu-port-207-final = gpu-port 27b215a7+). Push only to `fork`.

## Goal (PLAN.md 208-1)
G device-authoritative: remove the 206 interim per-CGR G round trip (D->H new G, H->D wrapped G for the precision reduction,
gpu_stab_mod.F90 ~l.590-630) and the 205 interim H->D/D->H of G around WRAPGRUP/DO (gpu_update_mod.F90 ~l.52, 227, 241);
a host copy of G only where a host reader needs it (host Obser at measurement slices; 207's GR H->D into G00_d, gpu_tau_mod.F90 ~l.418).

## Current step
Survey: find every host reader/writer of GR/G (main loop, Obser, wrap, tau pass) and design the minimal change. Nothing edited yet.

## Next step
Read gpu_stab_mod.F90 l.560-640, gpu_update_mod.F90 l.40-60/220-250, gpu_tau_mod.F90 l.410-540 and the host call sites of GR; write the
design (5-10 lines) here; implement; build both GPUs.

## Functional checks (as 207: per GPU, default + LOG), scripts docs/gpu/scripts/t207_final/{build,t0,obs,ctest}.sh
1 builds - 2 self-tests (ALF_GPU_SELFTEST=4, *_pp build for default) - 3 same-seed vs GPU=none <=1e-9 - 4 production observables vs CPU ref
(compare2 sum t^2; Intel Szz(pi,pi) group deviation of 207 is known/accepted) - 5 GPU=none ctest 53/53 (GNU, nividia).
Also report H<->D bytes per sweep before/after (the point of the ticket).

## Running
Nothing.

## Open problems
None.
