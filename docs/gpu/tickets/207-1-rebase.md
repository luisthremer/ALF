# 207-1 – 207 code on the current `gpu-port`

- **Phase / lane:** 2 (NVIDIA) · part of 207 (split by decision 34; context and history in `207-device-tau-pass.md`)
- **Owner agent:** fortran-gpu-dev (Sonnet)
- **Review:** numerics-reviewer (merge diff) · no harness acceptance (preparatory; 207 is accepted as a whole in 207-7)
- **Depends on:** –
- **Branch:** `gpu-port-207-r` = `gpu-port` `8efc6c31` + merge of the 207 WIP `c8d1c6f1` (manager; only docs conflicted)
- **Status:** open

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
