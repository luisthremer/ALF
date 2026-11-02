# 208-1 state (single source of truth; read only this file + CLAUDE.md)

Process: decision 42 = decision 41 style (functional checks only; no agents, no skills, serial, one step per turn, long jobs detached,
machine rules of CLAUDE.md). Branch `gpu-port-208-1` (from gpu-port-207-final = gpu-port 27b215a7+). Push only to `fork`.

## Goal (PLAN.md 208-1)
G device-authoritative: remove the 206 interim per-CGR G round trip (D->H new G, H->D wrapped G for the precision reduction,
gpu_stab_mod.F90 ~l.590-630) and the 205 interim H->D/D->H of G around WRAPGRUP/DO (gpu_update_mod.F90 ~l.52, 227, 241);
a host copy of G only where a host reader needs it (host Obser at measurement slices; 207's GR H->D into G00_d, gpu_tau_mod.F90 ~l.418).

## Design (survey done)
Facts: device G is `G_d(Ndim,Ndim,nf)` (gpu_update_mod; g_in/g_out copy H<->D around every WRAPGRUP/DO); gpu_stab_cgr works on st_Gn (new G),
st_Gt (wrapped G, H->D for the precision check) and copies st_Gn D->H into host GR every call; gpu_tau_pass copies GR H->D into G00_d and
G00 D->H after each CGR2_2; host readers of GR: main.F90 ham%Obser (l.~843-847, ~873-880, same in the DO branch), Hop_mod_Symm,
ham%GR_reconstruction (all per slice in the measurement phase), plus host-path code (Global/Langevin/tempering are refused by gpu_supported).
Protocol (lazy sync, two flags g_dev_ok / g_host_ok, module vars in gpu_update_mod, only when device updates are active; with
ALF_GPU_NOUPDATE=1 everything stays eager as today):
 - g_in: H->D only if not g_dev_ok. g_out: no D->H; g_host_ok = .false.
 - gpu_stab_cgr: precision reduction against G_d(:,:,nf) (the wrapped G, no H->D of st_Gt); then G_d(:,:,nf) = st_Gn device-to-device
   (omp target loop); no D->H of st_Gn; g_host_ok = .false.; the initial CGR (with_prec = .false.) the same. check_cgr (chk_level>=1)
   and the eps_prod fault path keep their eager copies (debug only).
 - new public gpu_g_host(GR): D->H of all flavours if .not. g_host_ok, then g_host_ok = .true.  Called in main.F90 right before the host
   readers (Obser / Symm / GR_reconstruction sites) when the GPU path is active.
 - gpu_tau_pass: GR -> G00_d from G_d (D2D) if g_dev_ok, else H->D; at its end G00 -> G_d (D2D) and host copy as today.
Expected: no G traffic in equilibration; in measurement 1 D->H per slice (needed until 208-2) instead of 2 H->D + 2 D->H.

## Current step
Design written. Implementing in Prog (gpu_update_mod, gpu_stab_mod, gpu_tau_mod, main.F90 call sites).

## Next step
Implement; build nv (GPU=nvidia nvhpc noMPI HDF5) and in (GPU=intel); then checks 1-5; report interim byte counters before/after.

## Running
Nothing.

## Functional checks (as 207: per GPU, default + LOG), scripts docs/gpu/scripts/t207_final/{build,t0,obs,ctest}.sh
1 builds - 2 self-tests (ALF_GPU_SELFTEST=4, *_pp build for default) - 3 same-seed vs GPU=none <=1e-9 - 4 production observables vs CPU ref
(compare2 sum t^2; Intel Szz(pi,pi) group deviation of 207 is known/accepted) - 5 GPU=none ctest 53/53 (GNU, nividia).
Also report H<->D bytes per sweep before/after (the point of the ticket).

## Running
Nothing.

## Open problems
None.
