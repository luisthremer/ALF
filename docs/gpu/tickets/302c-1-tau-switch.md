# 302c-1 – device τ pass on the B580 (`GPU=intel`)

- **Phase / lane:** 3 (Intel) · part of 302c (decision 34); runs in parallel with 207-3…207-5 (decision 36)
- **Owner agent:** fortran-gpu-dev (Sonnet)
- **Review:** numerics-reviewer (code, part a) · no own harness acceptance: the Intel items are judged in 207's
  `207.json` (207-6/207-7, decision 36)
- **Depends on:** 207-1 (207 code on `gpu-port`), 302b
- **Branch:** `gpu-port-302c` (from `gpu-port-207-r`)
- **Status:** T0 done (2026-11-01, commit `0e656fc8` on `gpu-port-302c`); numerics-reviewer (a) pending

## Goal
The 207 device τ pass (`Prog/gpu_tau_mod.F90`: CGR2_2 on the device, udvst kept on the device) runs under `GPU=intel`
on the B580 and agrees with the ifx `GPU=none` oracle as on NVIDIA.

## Context
- 207-1 added `tau_dev_allowed = .false.` under `ALF_LA_ONEMKL` (decision 22 until 302c). Before that guard, the
  unguarded merged code already ran the device τ pass on the B580 (207-1 run `i_def_6_try1`, `intel:~/alf-gpu/work/t207-1/`:
  200 passes, 1600 CGR2_2 calls, Precision tau mean 5.45e-12 vs 6.27e-12 host; equal-time datasets bit-identical; not
  checked against an oracle). So the alf_la routines it calls exist under oneMKL; whether each is right on the B580 is
  this ticket's question.
- The `alf_la` interface is frozen while both lanes run (CLAUDE.md). A needed backend fix inside the oneMKL
  implementation: stop and report first.
- 207-3 edits `gpu_tau_mod.F90` at the same time (fault hook, refusals on early returns). Keep your change to the guard
  region and the oneMKL-specific code; the manager merges.

## Tasks
1. Open the switch: the τ pass runs under `GPU=intel` (remove or flip the guard, info line as on NVIDIA: `ACTIVE`).
2. Check every alf_la call of `gpu_tau_mod.F90` on the oneMKL backend (workspace queries, QR with pivoting, the
   CGR2_2 pieces, precision reduction); fix what is wrong inside the allowed files.
3. T0 on dev seeds on the B580 (`ONEAPI_DEVICE_SELECTOR=level_zero:0`, load check before every run, ≤ 8 cores):
   `ALF_GPU_SELFTEST=4` all configs default and LOG (tolerance 1e-10 as on NVIDIA, report the largest values),
   fault controls (`perm`, `prop`, `prec`, `blocks`, `eps:1e-9` killed; `eps:1e-11` survives; invalid values and debug
   variables refused), same-seed `hub6x6_cb` and `hub12x12_cb` default and LOG vs the ifx `GPU=none` oracle
   (≤ 1e-9; LOG 12×12 compare with the known ticket-213 excess on NVIDIA and say plainly if it differs),
   `Ltau = 0` runs bit-identical to `gpu-port`, `GPU=none` (ifx and GNU) identical to `gpu-port`.
4. Report the Precision-tau means device vs host on the B580.

## Files allowed
`Prog/gpu_tau_mod.F90` (guard and oneMKL-specific code only), `Prog/gpu_stab_mod.F90` only if the udvst slots need a
oneMKL-specific change (say why). Work dir `intel:~/alf-gpu/work/t302c/` (TMPDIR inside). Unattended runs allowed
(decision 36).

## Done when
Report with run dirs; numerics-reviewer approves (a).

## Out of scope
Accuracy gates and QREF on the B580 (302c-2); harness (207-5); PVC performance.

## Result (2026-11-01, T0 on dev seeds, development evidence)
- Diff: 207-1's interim guard removed (−14 lines in `Prog/gpu_tau_mod.F90`); the alf_la oneMKL calls needed no change.
- Self-tests `ALF_GPU_SELFTEST=4`: 12/12 config × path pass on the B580, max 8.4e-12 (8 % of 1e-10). Fault controls
  `perm`, `prop`, `prec`, `blocks`, `eps:1e-9` killed on both paths; `eps:1e-11` survives; invalid values and 10 debug
  variables refused.
- Same-seed vs the ifx `GPU=none` oracle: 6×6 default 2.31e-10, LOG 2.90e-10; 12×12 default 4.81e-10; **12×12 LOG
  1.128e-9 FAILS ≤ 1e-9**, in the equal-time datasets (first bin 27); the host-τ base gives exactly the same value and its
  equal-time datasets are bit-identical, τ datasets 2.29e-10 (base 3.02e-10). This is the ticket-213 LOG excess
  (NVIDIA 1.214e-9), present on Intel too, not caused by the τ pass; it stays with 213 and is `report_only` in 207.json
  as on NVIDIA.
- Identity: `Ltau=0` and the `Ltau=1` equal-time datasets bit-identical to the base; runs deterministic; ifx `GPU=none`
  identical (GNU `GPU=none` not run: the only changed file is compiled under `ALF_GPU` only).
- Precision-tau means device/CPU-oracle 0.79–1.01 (β = 10 short runs 1.01 / 0.98).
- Evidence `intel:~/alf-gpu/work/t302c/` (no `provenance.json` in T0, as in 207-1); scripts `docs/gpu/scripts/intel/t302c/`.
