# 207-1 – 207 code on the current `gpu-port`

- **Phase / lane:** 2 (NVIDIA) · part of 207 (split by decision 34; context and history in `207-device-tau-pass.md`)
- **Owner agent:** fortran-gpu-dev (Sonnet)
- **Review:** numerics-reviewer (merge diff) · no harness acceptance (preparatory; 207 is accepted as a whole in 207-7)
- **Depends on:** –
- **Branch:** `gpu-port-207-r` = `gpu-port` `8efc6c31` + merge of the 207 WIP `c8d1c6f1` (manager; only docs conflicted)
- **Status:** done, pending numerics-reviewer approval of the merge diff (fix `2579bc5d` on `gpu-port-207-r`)

## Goal
The 207 device τ pass, unchanged in behaviour, builds and runs on top of the current `gpu-port` (212, 302a, 302b and the
harness merged since the WIP's base `3f8ec8fa`).

## Context
- The merge is automatic for all code (`gpu_stab_mod`, `gpu_state_mod` auto-merged); 302b added the `GPU=intel`
  device paths and a τ switch: under `GPU=intel`, `Ltau = 1` runs do `TAU_M` on the host (decision 22) until 302c.
  This must stay so after the merge.
- WIP state: kernel self-tests and same-seed tests passed on the old base (ticket table "Evidence status per acceptance
  item"); dev5 QREF runs in `nividia:~/alf-gpu/work/t207/dev5/`.

## Tasks
1. Build `GPU=nvidia` (default with `QR_PREPIVOT` oracle build, LOG), `GPU=none` (GNU), on `intel` `GPU=intel`.
2. T0 on dev seeds: acceptance-1 self-tests and fault controls (all configs, default and LOG); same-seed `hub6x6_cb`
   and `hub12x12_cb` vs the CPU oracle (≤ 1e-9); `GPU=none` preprocessed-source identity and `data.h5` identity vs
   `gpu-port`; on `intel` a 6×6 `Ltau = 1` run with the τ pass on the host (info line says so).
3. One QREF run (dev1, both paths, as in dev5): the S_ens numbers equal dev5's (expected bit-identical; if not, report
   the difference and why).

## Files allowed
Only fixes needed to make the merged code build and behave as before, in `Prog/gpu_*.F90`, `Prog/tau_m_mod.F90`,
`Prog/main.F90`, `Prog/control_mod.F90`, `Prog/Makefile`. Anything else: stop and report.

## Done when
Report with stamped runs (run dirs) for every task; manager check; numerics-reviewer approves the merge diff.

## Out of scope
Pipeline fixes (207-2), new gates and controls (207-3), β = 10 (207-4), harness (207-5).

## Result (2026-11-01, T0 on dev seeds, development evidence)
- **One fix needed:** the merged code ran the device τ pass under `GPU=intel` (seen on the B580, run `i_def_6_try1`),
  against decision 22 and 302b's `forbidden_lines`. `gpu_tau_setup` now returns early under `ALF_LA_ONEMKL` with the
  line `GPU tau pass (207): host TAU_M (GPU=intel: …, decision 22); udvst exported D->H (INTERIM 206)`; CUDA path
  unchanged. Commit `2579bc5d` (+13 lines in `Prog/gpu_tau_mod.F90`).
- **Self-tests** (`ALF_GPU_SELFTEST=4`): 12/12 config × path pass; table identical to the WIP's on the 212 base; largest
  forced-branch value 9.25e-11 (92 % of 1e-10, default b10, D1 ×1e-3). Fault controls: 5 faults × 2 paths rc=1,
  `eps:1e-11` rc=0, 4 invalid values rc=1, 10 debug variables refused in production.
- **Same-seed vs CPU oracle:** 6×6 default 1.49e-10, LOG 2.54e-10; 12×12 default 3.10e-10; **12×12 LOG 1.214e-9 FAILS
  ≤ 1e-9** (equal-time datasets, from bin 184; τ datasets 1.27e-10). Identical to the WIP value before the merge: the
  known LOG excess of ticket 213 (`report_only`), not a merge effect.
- **Identity:** merged vs `0df7217c` bit-identical (nividia 6×6; intel default and LOG 6×6); `GPU=none` 27/27
  preprocessed comparisons identical, `data.h5` 3 configs × 38 datasets identical, ctest 53/53 (the ticket table's
  "52" is stale; 53 is the base's count).
- **intel 6×6 `Ltau=1`:** host τ pass, info line as above; vs ifx `GPU=none` oracle default 1.45e-10, LOG 3.36e-10.
  The device τ pass on the B580 stays unverified until 302c.
- **QREF dev1 (both paths):** `data.h5`, all 280 dump hashes, `qref.txt` and aggregated S_ens identical to dev5
  (default dev +0.084 [+0.051, +0.118], LOG −0.012 [−0.043, +0.019]). The current `t207_d5_qref.F90` gives 0.00e+00
  difference on dev5's kept dumps, so the lost v2 source `b161b985` is not needed to reproduce dev5.
- **Evidence:** `nividia:~/alf-gpu/work/t207-1/`, `intel:~/alf-gpu/work/t207-1/` (indexes `logs/t0_index_*.json`;
  T0 scripts wrote no `provenance.json`, the QREF runs did). Scripts: `docs/gpu/scripts/nividia/t207_1/`.
