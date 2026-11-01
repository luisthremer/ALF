# 207-6 – pre-registration of `207.json`

- **Phase / lane:** 2 (NVIDIA) · part of 207 (decision 34)
- **Owner agent:** manager
- **Review:** numerics-reviewer (acceptance file, part b)
- **Depends on:** 207-1 … 207-5, 302c-1, 302c-2 (decision 36)
- **Status:** open

## Goal
`testsuite/gpu_reference/acceptance/207.json` pre-registered on the 207 branch (merged with `gpu-port` including
207-5), approved by review (b).

## Content
- Acceptance 1 (kernel self-tests, faults), 2 (same seeds, ≤ 1e-9), 4 (compare2 with τ-resolved observables),
  5 (`GPU=none` identity, ctest count) as in the ticket, re-derived for the new base.
- Acceptance 3 per decision 33: S_ens(pTT) ≤ 0.15 and the pointwise gate (bounds and control x from 207-3/207-4),
  1.5× Precision-tau mean rule, paired mean and count report-only; β = 4 and β = 10; frozen ensemble lists; control
  runs on fresh seed sets distinct from the main sets; false-alarm budget written down.
- Statement: WRAPUR accuracy rests on acceptances 1 and 2 plus the mean rule.
- Acceptance 6 (Intel), decision 36: 207 is judged on the B580 too, with the device τ pass from 302c-1: self-tests and
  faults, same-seed vs the ifx `GPU=none` oracle, the accuracy gates with an MKL host member in the reference ensemble
  (302c-2's calibration), compare2 on the B580.
- Bounds fixed by the manager by decision 37's rule; GPU time per decision 38.
- `207.draft.json` retired (the shipping check allows only an added `207*.json` between base and pre-registration).
- Run time and GPU budget (decision 11) from 207-4's numbers.

## Done when
Review (b) approves; PREREG_SHA recorded in the ticket.
