# 213 – STABLOG same-seed excess above 1e-9 at 12×12 (206 code)

- **Phase / lane:** 2 (NVIDIA lane)
- **Owner agent:** fortran-gpu-dev (Sonnet; Opus on the manager's proposal: numerically delicate)
- **Review:** numerics-reviewer (code, acceptance file, scripts) · **Verification:** acceptance-verifier
- **Depends on:** 212 (deterministic reductions, debug-variable gating)
- **Status:** draft (2026-10-29), not started. Split off from 207 by the owner's decision of 2026-10-29.

## Goal
Explain, and fix if it is a defect, why the full-device STABLOG path deviates from the CPU LOG oracle by more than 1e-9
(relative, per bin) on `hub12x12_cb`, 300 sweeps, same seeds — or show with a pre-registered, calibrated criterion that
the deviation is rounding amplified by conditioning and not a defect.

## Context
- Measured (206 acceptance and again in 207, 2026-10-28/29): first bin above 1e-9 = 184, maximum 1.214e-9, in the
  equal-time and scalar datasets only; the τ-resolved datasets stay ≤ 1.3e-10; no bin above 1e-4 (no divergent
  decision). The default path stays ≤ 3.1e-10 on the same config. 6×6 LOG: ≤ 2.5e-10.
- The same-seed bound 1e-9 stays (owner decision 3 of 2026-10-29); it is not relaxed here.
- Leads: the STABLOG decompose on the device (pre-pivoted QR with `k_lognorm` and the log-scale bookkeeping,
  `gpu_stab_mod.F90`); host ZGEQP3 in CGR; growth of the deviation with the sweep count (bin 184 of 300).
- 212 makes the reductions deterministic; re-measure first, then localise (per-sweep deviation curve, swap device
  steps for host ones as in 207's diagnostic, 1-ulp sensitivity of the CPU LOG chain itself).

## Files allowed
- `Prog/gpu_stab_mod.F90` (STABLOG decompose and CGR paths), `docs/gpu/scripts/nividia/t213_*`,
  `testsuite/gpu_reference/acceptance/213.json`. Anything else: stop and report.

## Acceptance (pre-registered in `213.json` before implementation)
- Same-seed `hub12x12_cb` 300 bins, LOG: no bin above 1e-9 (unchanged bound), all datasets; default path unchanged
  (bit-identical where the 212 determinism holds).
- Or, if the cause is conditioning: a criterion pre-registered **before** new data, calibrated with the CPU LOG chain
  against itself under a 1-ulp perturbation (null) and a planted defect; it needs the owner's quoted approval.
- Kernel self-tests and fault controls of 206 unchanged; GPU=none identity.

## Out of scope
The default path; the τ pass (207); performance.
