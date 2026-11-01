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
- Intel coverage (302c-1 review): include same-seed runs with N_FL_eff = 1, Symm and flux under `GPU=intel`.
- Bounds fixed by the manager by decision 37's rule; GPU time per decision 38.
- `207.draft.json` retired (the shipping check allows only an added `207*.json` between base and pre-registration).
- Run time and GPU budget (decision 11) from 207-4's numbers.

- Gate spec (207-2, approved 2026-11-01): `keep` must satisfy the keep rule (40 is refused at stride 3) and
  `min_reprocessed` ≥ P (number of block/flavour positions) so that every stratum is re-computed at least once per run;
  the spec pins seed hashes, so it is written after the fresh seeds are drawn; `--spec-sha` goes into the acceptance
  file.

- Keep rule (207-3, adopted 2026-11-01): gcd(m, P) = 1 and m ≤ (ndumps − 1)/(P − 1), m = cover·keep, so the kept
  dumps cover all P positions (β = 4: KEEP 39; β = 10: KEEP 11).
- Intel controls are the sweep output fault (`FAULT_TAU=eps:<x>`, `FAULT_TAU_AT`); var8 / NS=0 / TAUSWAP=8 are no defect
  on oneMKL (302c-2). Intel specs use host_member `nvh` (label fixed by the offline program).

- Spec instance fields (207-5): `binary_sha256`, `prog_sha`, `seeds` = `"@harness"`; `parameters_sha256` concrete.

## Done when
Review (b) approves; PREREG_SHA recorded in the ticket.
