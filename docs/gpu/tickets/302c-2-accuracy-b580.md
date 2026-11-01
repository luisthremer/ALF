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
