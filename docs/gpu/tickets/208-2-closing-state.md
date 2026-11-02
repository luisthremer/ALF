# 208-2 state (single source of truth; read only this file + CLAUDE.md)

Process: decision 42 (functional checks as 207, no agents/skills, serial, one step per turn, detached long jobs, machine rules).
Branch `gpu-port-208-2` (from gpu-port 6f883f4f = 208-1 closed). Push only to `fork`.

## Goal (PLAN.md 208-2)
Equal-time Hubbard observables (`Obser`) on the device, so the per-slice D->H of G (the remaining 96000 copies / 31.85 GB per production
run of 208-1) disappears.

## Survey (done)
- `Obser` (Hamiltonians/Hamiltonian_Hubbard_smod.F90 l.695): scalars Obs_scal(1..4) = kin, pot, rho, E (kin via
  Predefined_Hoppings_Compute_Kin: Z_kin = sum coef * GRC(I1,J1,nf) over bond terms (Z, conj Z) and T_Loc diagonal terms; pot = sum_i
  GRC(i,i,1) GRC(i,i,dec) U (bilayer: U2 for orbital 2); rho = N_SUN sum GRC(i,i,nf)); local obs only if |ham_h0| >= Eps_small;
  correlations via Predefined_Obs_eq_{Green,SpinMz|SpinSUN,Den}_measure (Predefined_Obs_mod.F90): each is a scatter-add over (I1,J1) into
  Obs_Latt(imj(I,J),1,no_I,no_J) with a weight ZP*ZS from GRC/GR products. GRC(i,j,nf) = delta_ij - GR(j,i,nf).
- Obs_scal/Obs_eq/Obs_tau are `public` module variables of the Hubbard module; the gpu_* modules compile AFTER the Hamiltonians
  (Makefile l.14), so the Hubbard code cannot call them: the hook sits in main.F90 (as gpu_g_host) and a new gpu_obs_mod (after
  gpu_tau_mod) reads the public Obs_* of the Hubbard module.
## Design
1. gpu_obs_mod: setup (once, after state setup): flatten the hopping into arrays (I1,J1,nf,coef) replicating Compute_Kin (bonds: (I1,J1,Z),
   (J1,I1,conj Z); T_Loc: (I1,I1,T)); device tables List (I1 -> site, orbital), jcell(I, imj) (the unique J with imj(I,J) = imj), Invlist;
   device accumulators with the shape of Obs_eq(k)%Obs_Latt and of Obs_scal(k)%Obs_vec.
2. gpu_obs_measure(PHASE, Ntau, Mc_step_weight): reads G_d (device), computes scalars (reductions) and the correlations by a GATHER per output
   element (imj,no_I,no_J) (deterministic, no atomics), adds ZP*ZS-weighted results to the device accumulators; host increments the counters
   (Obs%N, Obs%Ave_sign) exactly as Obser does.
3. Flush: before the host normalises (Hubbard Pr_obs call site in main / Init_obs reset) D->H the accumulators and ADD them into the host
   Obs arrays, then zero them (once per bin, not per slice).
4. Support gate gpu_obs_supported: Hubbard, N_FL in {1,2}, |ham_h0| < Eps_small, no Symm, no reconstruction_needed, no Projector, device updates on;
   else the host Obser stays (gpu_g_host before it, as 208-1).
5. Self-test (ALF_GPU_DEBUG=1, ALF_GPU_SELFTEST=5): on the first measurement slice run host Obser into scratch copies and the device path;
   compare to <= 1e-10 relative; planted fault (ALF_GPU_FAULT_OBS) must be caught.
## Current step
DONE (2026-11-02, decision 42): 208-2 closed; merged into gpu-port, pushed to fork.
Implemented as designed in Prog/gpu_obs_mod.F90 (setup from Alloc_obs of the Hubbard module, gpu_obs_measure from Obser, gpu_obs_flush before
Pr_obs in main.F90; Obs_scal/Obs_eq made public for ALF_GPU builds in Hamiltonian_main_mod). Fallback to the host Obser when unsupported or
ALF_GPU_NOOBS=1. No separate self-test (decision 42: the same-seed check compares all 'obser' datasets).
Lesson (nvfortran): a complex scalar (the sign/phase weight) inside the device loops was read as garbage (real part NaN); passed as two real
scalars it works. Complex private temporaries and reductions were fine.

## Results (work/t208-2 on nividia and intel; scripts docs/gpu/scripts/t208_2/)
1 Builds: nividia 5, intel 5: PASS (after the name clash inv/jc/uo with MyMats and the PGI-only public Obs_* were fixed).
2 Self-tests (=4): PASS on both, def and LOG.
3 Same-seed hub6x6_cb_fd max dev: nividia 1.70e-10 / 2.54e-10, intel 2.32e-10 / 2.90e-10 (def / LOG): PASS, identical to 207 and 208-1.
4 Observables hub12x12_cb vs CPU ref: nividia sum t^2 = 9.2 (p 0.60) def and LOG; intel 141.7 (the accepted 207 Szz(pi,pi) group): same as 207.
  Wall: nividia def 703 s, LOG 841 s (208-1: 696 / 812), intel 419 / 409 s (423 / 413): no speed gain at this size.
5 GPU=none ctest (GNU): candidate 53/53, upstream master 39/39: PASS.
Traffic (INTERIM 205 counter, obs_def): G copies D->H 96000 -> 1200 (the tau pass, once per sweep), H->D 0; bytes 31.85 GB -> 0.40 GB.

## Next step
208-3 (tau-resolved observables OBSERT on the device; removes 207's per-slice GT0/G0T/GTT D->H) when the owner says continue.

## Running
Nothing.
