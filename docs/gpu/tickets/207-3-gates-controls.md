# 207-3 – pointwise accuracy gate and graded controls (β = 4)

- **Phase / lane:** 2 (NVIDIA) · part of 207 (decision 34)
- **Owner agent:** fortran-gpu-dev (Opus: numerically delicate)
- **Review:** numerics-reviewer (code + calibration)
- **Depends on:** 207-1, 207-2
- **Status:** open

## Goal
Calibrated, on dev seeds at β = 4, the two accuracy gates of decision 33 and controls that show each gate is sharp:
- S_ens(pTT) ≤ 0.15 (dev5 numbers: device +0.03…+0.08 default, −0.06…−0.01 LOG; correct variants ≤ +0.12; NS=0 /
  TAUSWAP=8 controls ≥ +0.19);
- a pointwise gate S_ens(emax) or S_ens(eF) (dev5: device −0.14…−0.10 default, −0.07…−0.05 LOG; highest correct
  variant +0.04…+0.07 / +0.06…+0.10). The existing controls do not move it, so it needs its own planted control.

## Tasks
0. All QREF runs use stride 3 (decision 35: both spin flavours, all blocks and branches), lengthened to ~210 points
   with B. The dev5 numbers above are flavour-1-only; recalibrate everything (S_ens(pTT) bound included) on stride-3 runs
   and report whether flavour 2 behaves differently.
1. A graded fault on the device CGR2_2 output (`ALF_GPU_FAULT_TAU=eps:<x>`, refused in production runs as in 212,
   declared control runs only), and its predicted response.
2. Calibrate on dev seeds: device, correct variants, NS=0 / TAUSWAP=8, graded fault at a far and a near-threshold x
   (as decisions 26/27), for both gates and both paths.
3. Attribute the default path's top-of-band position (QVAR host-twin masks), as the review recommends.
4. Make τ checks fail loudly when the device τ pass is off (207-1 review nit): on every early return of
   `gpu_tau_setup` (GPU=intel per decision 22, Ltau=0, Projector, reconstruction, no `gpu_stab_on`), stop with an error if
   `ALF_GPU_SELFTEST=4`, `ALF_GPU_TAUTIME`, `ALF_GPU_FAULT_TAU` or any `ALF_GPU_TAUDIAG*` is set; today the self-test
   exits 0 and the fault prints "armed" but perturbs nothing, so an rc=0 cannot tell right from wrong. Show both:
   a run with the τ pass off and SELFTEST=4 exits non-zero; a normal run is unchanged.
5. Propose the pointwise bound and the control x values with their reasons. The manager brings them to the owner; they
   are fixed only in the pre-registration (207-6).

## Files allowed
`Prog/gpu_tau_mod.F90` (debug-gated fault hook and the task-4 refusals only; `gpu_state_mod.F90` if the refusal must live there, say why), `docs/gpu/scripts/nividia/t207_*`. Work dir
`nividia:~/alf-gpu/work/t207-3/`.

## Done when
Calibration tables from stamped runs; numerics-reviewer approves code and calibration.
