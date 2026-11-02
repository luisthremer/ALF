# 302c-2 – accuracy calibration of the device τ pass on the B580

- **Phase / lane:** 3 (Intel) · part of 302c (decision 34); feeds 207-6 (decision 36: 207 is judged on both GPUs)
- **Owner agent:** fortran-gpu-dev (Opus: numerically delicate, as 207-3)
- **Review:** numerics-reviewer (calibration)
- **Depends on:** 302c-1 (done), 207-2 pipeline; method parallel to 207-3
- **Branch:** `gpu-port-302c-2` (from `gpu-port-302c`)
- **Status:** open

## Goal
The QREF accuracy numbers of decision 33 (S_ens(pTT); emax/eF statistics as 207-3 defines them) for the device τ pass on
the B580, β = 4 and β = 10, default and LOG, with a reference ensemble that has an ifx/MKL host member in place of `nvh`;
the controls (`ALF_GPU_TAU_NS=0` default path, `ALF_GPU_TAUSWAP=8` LOG path, `var8`) shown to be caught.

## Tasks
1. Build on `intel` the QREF-capable `GPU=intel` binaries (default with `QR_PREPIVOT`, LOG), and the offline binary128
   program plus GNU oracle objects (gfortran 15.2 is on `intel`), with sha256 recorded.
2. Dev-seed QREF runs, stride 3 (decision 35), ~210 points with B per run, `hub12x12_cb` β = 4 (4 seed sets per path)
   and β = 10 (2 seed sets per path), each processed with `t207_d5_proc.sh` (KEEP 39) before the next run.
3. Name and justify the Intel host member (the in-run host CGR2_2 under ifx/MKL, `ALF_GPU_TAUDIAG_QHOST=1`), and the
   frozen Intel reference lists (which n-variants exist under ifx).
4. Tables from `t207_d5_agg.py` (generated, not typed): device, each correct variant, var8, controls; per seed set and
   pooled, leave-one-out; flavour 1 vs 2. Compare with NVIDIA (207-3) and say plainly whether one bound can serve both.
5. If the scripts need an Intel-specific change (member names, env allow-list), describe it; the gate code is being
   redesigned on `gpu-port-207-2` — do not edit `t207_d5_*`.

## Files allowed
`docs/gpu/scripts/intel/t302c/*` (new scripts). No `Prog/` changes (stop and report if one is needed). Work dir
`intel:~/alf-gpu/work/t302c-2/`, TMPDIR inside; unattended runs allowed (decision 36).

## Done when
Calibration tables from stamped runs (provenance.json, stamp.txt); numerics-reviewer approves the calibration.

## Result, round 1 (2026-11-01; dev seeds, development evidence, not reviewed)
Scripts on `gpu-port-302c-2` `6a5c499b` (`docs/gpu/scripts/intel/t302c/t302c2_*`), source `0e656fc8`, builds ifx 2026.2 +
oneMKL 2026.2 and gfortran 15.2 + OpenBLAS; 16 runs (β = 4: 4 sets per path + controls; β = 10: 2 sets per path +
controls), all rc 0, 0 QREFERR, 86 B580 GPU-minutes. Evidence `intel:~/alf-gpu/work/t302c-2/` (821 MB).

S_ens(pTT), reference list `ens`, pooled:

| | β4 def | β4 LOG | β10 def | β10 LOG |
|---|---|---|---|---|
| device | −0.032 | −0.117 | −0.122 | −0.235 |
| var8 | +0.034 | −0.008 | −0.057 | −0.051 |
| control (NS=0 / TAUSWAP=8) | +0.030 | −0.010 | −0.050 | −0.053 |
| highest correct member per set | 0.068–0.100 | 0.059–0.108 | 0.085, 0.087 | **0.174 (n5), 0.138** |

- **The device τ pass on the B580 is accurate:** below every correct host variant at every β and path.
- **The NVIDIA controls are no defect on oneMKL.** Device ZUNMQR without the Newton–Schulz step lies inside the correct
  band (NVIDIA: +0.19…+0.37); NS=0 raises Precision-tau by only 1.7–7.4 % (NVIDIA 17 %). On Intel no available control
  fails, so the Intel gates need a backend-independent planted defect acting on the device result in the dumps: the
  output fault of 207-3. The existing `ALF_GPU_FAULT_TAU=eps:<x>` acts only in the SELFTEST=4 comparison and is no
  such control.
- **A correct variant exceeds 0.15 at β = 10 LOG** (n5 +0.174, n4 +0.173 in one set; pooled n5 +0.156). Hypothesis (not
  verified): the row-permuted offline variants (flag 4) on row-graded STABLOG matrices; they are CPU members, so
  NVIDIA 207-4 should show the same. Under decision 37 the β = 10 LOG bound is set by the rule from these data, not
  fixed at 0.15; the reference list stays frozen.
- **Keep coverage at β = 10:** KEEP 39 keeps 6 of 234 dumps, fewer than the P = 20 positions (passed to 207-3/207-6).
- Script follow-ups (207-2/207-5): neutral name for the host member (label `nvh` hard-coded; Intel member `ifh`),
  `proc.sh` must not require `~/alf-gpu/env/gnu.sh`.

Next (round 2, after 207-3's output fault lands): rebuild at the merged sha, run the graded output fault on the B580
and calibrate per β and path by decision 37's rule; then numerics review.

## Result, round 2 (2026-11-02; dev seeds dev5–dev9, development evidence, not reviewed)
Scripts `d0faa7ac` (`t302c2b_*`), binaries from 08559f0b; evidence `intel:~/alf-gpu/work/t302c-2b/`
(`tab/summary.md`, `tab/gates_b{4,10}_{def,log}.*`, generated). 28 QREF runs, 139.6 GPU-min, no lock waits.
- T0 at 08559f0b passes on the B580 (self-tests max 8.4e-12; same-seed 6×6 vs ifx GPU=none 2.31e-10 / 2.90e-10; fault
  classes killed; the output fault acts linearly, 9.755e-8 per 1e-10). TAUSWAP / TAU_NS / TAUTEST_V not refused
  with the τ pass off at this sha: expected, 207-3's closing round was not in it.
- β = 4: the proposed gates separate except the near control 2e-15 (def +0.119, LOG −0.067, not killed); 1e-14 is
  killed on both paths (+0.625 / +0.332). The B580 device sits lower than NVIDIA's (def pTT −0.032, LOG −0.125).
- β = 10: **emax pooled / by_stratum / q95 are blind to eps 1e-12 and the localized 2e-12** (inside the correct band;
  only 1e-10 is caught); the localized fault is caught by pTT by_stratum only (+2.371 vs correct max +0.369).
  **β = 10 LOG: a correct variant exceeds pTT 0.15** (n100004 pooled +0.172, single run +0.207); the four highest
  correct members are always the flag-4 row-permuted variants (confirmed on fresh seeds). emax by_stratum / q95 at
  β = 10 LOG pass the correct variants with margins of only 0.018 / 0.015 (2 main sets).
- Fault runs kept KEEP 117 (gcd 1 with P; evidence budget); statistics use the full qref.txt, only the re-processing
  sample is smaller. Accepted for calibration evidence; T2 runs follow the keep rule.

### Manager decisions (2026-11-02)
- 207-3 final (b2c71200) merged into `gpu-port-302c-2` as `d3e82a36` (code changes identical to 428576dd..b2c71200).
- Control sizes are set per β (decision 37 covers control fault sizes): round 3 measures emax/localized ladders and
  larger pTT controls at β = 10, flavour-only controls, and two more β = 10 main sets per path, at d3e82a36, after a
  bit-identity check that round 2's main runs carry over to the merged code.
- The β = 10 LOG pTT bound is set from data by decision 37 (not 0.15) with the same rule on both GPUs; if no bound
  gives correct margin ≥ 0.05 and control margin ≥ 0.10, work stops for the owner.
