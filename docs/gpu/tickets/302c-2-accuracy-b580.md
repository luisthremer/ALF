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
