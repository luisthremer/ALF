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

## Result, main round (2026-11-01; dev seeds, not reviewed)
Commit `gpu-port-207-3` `428576dd`. Evidence `nividia:~/alf-gpu/work/t207-3/` (`logs/cal_*`, `logs/report_cal_*`, `gate/`).
- Code: `ALF_GPU_FAULT_TAU=eps:<x>` now also perturbs the four device CGR2_2 outputs in the sweep (1 + x·r, fixed r per
  element and block; debug only); `ALF_GPU_FAULT_TAU_AT=<nst>:<nfe>` localizes it; every early return of `gpu_tau_setup`
  refuses SELFTEST=4, TAUTIME, FAULT_TAU(_AT) and the TAUDIAG variables. Normal runs bit-identical to base (data.h5).
- β = 4, stride 3, KEEP 39, cover 1: proposed gates by decision 37's rule (correct-variant max / control min):
  S_ens(pTT) pooled ≤ 0.15 (def +0.113 per run / 2e-15 +0.204; LOG +0.061 / TAUSWAP +0.211, var8 ≥ +0.218);
  emax pooled ≤ 0.15 (1e-12 +0.25…+0.275); emax by_stratum ≤ 0.30 (localized 2e-12 @ (2,2) +0.56…+0.59);
  emax q95 ≤ 0.60. Not separable, report-only: pTT by_stratum, pTT q95, S_ens_max.
- **FAILED item:** on LOG the graded fault 2e-15 is not killed by pTT (+0.032); closing round runs 1e-14 and 5e-15.
- Var8 on LOG is a real defect at full run length (+0.218…+0.252); dev5's +0.123 came from the biased kept subset.
- Attribution: the device pTT excess (+0.066) comes from device rounding when Q is formed and applied (var4/var256
  bring it to cpu level); the device QR makes the device more accurate pointwise (var16).
- First genuine gate runs: approved gate 4118a1de end to end with KEEP 39 → PASS, re-processing 1034/1034 equal;
  by_stratum with 10 and 8 strata shown.
- Manager decisions (2026-11-01): localized fault stays `FAULT_TAU_AT` (gate table adapted in 207-5); `tau_off_refuse`
  also refuses TAUSWAP, TAU_NS, TAUTEST_V (closing round); keep rule adopted: gcd(m, P) = 1 and
  m ≤ (ndumps − 1)/(P − 1), m = cover·keep (β = 4: KEEP 39; β = 10: KEEP 11).

## Result, closing round (2026-11-01; dev seeds, not reviewed)
Commit `9a35bdd6`. Evidence `nividia:~/alf-gpu/work/t207-3b/` (`logs/summary_gates_final.{txt,json}`, generated).
- `tau_off_refuse` also refuses TAUSWAP, TAU_NS, TAUTEST_V (rc = 1 on τ-off input, both paths); normal runs
  bit-identical to the previous binaries.
- LOG near pTT control: eps 5e-15 killed (+0.201, 95 % interval above 0.15, margin +0.051), 1e-14 killed (+0.425);
  2e-15 not killed (+0.032). Default near control 2e-15 +0.204 (margin +0.054).
- emax q95 ≤ 0.60 misses the localized fault by 0.014 / 0.030; only by_stratum catches it.
- Open for the review and 207-6: which near controls to pre-register (must be killed reliably once on fresh seeds).

## Review (2026-11-01)
numerics-reviewer: (a) code **approve** (nits: comment placement, "pseudo-random" wording — the pattern is a
deterministic sawtooth, header sentence on SELFTEST=4 + eps). (b') calibration **changes required**:
1. near controls too thin for one kill on fresh seeds (predictive sd ≈ 0.027): def 2e-15 margin +0.054 (miss ≈ 2.3 %),
   LOG 5e-15 +0.051 (≈ 3 %), TAUSWAP=8 +0.061 (≈ 1 %, fixed); together ≈ 6 % that some control is not killed. Require
   margin ≥ ~0.10: LOG 1e-14, def ~3e-15–4e-15 confirmed on two seed sets; state TAUSWAP's miss risk and a control-miss
   budget in 207-6;
2. by_stratum / q95 correct-variant margins must be shown at the T2 pooling (per run up to +0.159 / +0.444 def,
   +0.190 / +0.474 LOG; the bounds 0.30 / 0.60 still pass); one minimum-stratum-size rule;
3. flavour gap: a var8-size defect in flavour 2 only gives pooled LOG S_ens(pTT) ≈ +0.10 and passes; a per-flavour max
   gate looks separable — calibrate it with a flavour-only control or state the gap to the owner.
Device values credible (attribution checked); generated numbers reproduced by an independent parser.
Fix round 207-3c running (FAULT_TAU_AT=0:<nfe> flavour-only fault; ≤ 55 GPU-min, counted against 207-4's pot so
207 calibration stays within decision 38's 2 × 4 h; 207-3 alone used ≈ 4.05 h).
Merge note: `gpu-port-302c-2` (08559f0b) merged 428576dd; merge the final 207-3 before Intel evidence counts as
final-code evidence.
