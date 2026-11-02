# 207-7 – candidate, T1, T2, verification, merge

- **Phase / lane:** 2 (NVIDIA) · part of 207 (decision 34)
- **Owner agent:** fortran-gpu-dev (T1), manager (T2), acceptance-verifier, numerics-reviewer (a)
- **Depends on:** 207-6
- **Status:** closed 2026-11-02 (decision 41)

## Goal
207 closes under VERIFICATION.md: frozen candidate, review (a), T1 on dev seeds, manager T2 on fresh seeds with all
controls killed, audit of the T1 report, independent verifier agrees, merged into `gpu-port`, owner's gate, evidence
archived.

## Done when
The owner approves the gate.

## Closing (decision 41)
Closed 2026-11-02 under decision 41 (owner: "Accept." on the Intel check-4 result): functional checks only, no pre-registration/T2/verifier/review. Candidate branch `gpu-port-207-final`, merged into `gpu-port`. Results: docs/gpu/tickets/207-closing-state.md.

Closing result (2026-11-02, decision 41, evidence in work/t207-final/ on nividia and intel):
- Builds GPU=nvidia (nvhpc) and GPU=intel (ifx + oneMKL), default and LOG: pass. Self-tests (ALF_GPU_SELFTEST=4): pass on both.
- Same-seed GPU vs GPU=none, 6x6, max dev: NVIDIA 1.70e-10 / 2.54e-10, B580 2.32e-10 / 2.90e-10 (default / LOG), limit 1e-9.
- Observables, 12x12 beta 4, vs CPU reference (compare2.py): NVIDIA sum t^2 = 9.2 (chi2_11 p = 0.60), both paths; B580 sum t^2 = 141.7,
  Szz(pi,pi) group 6.3 sigma low (3 correlated observables), other 8 ok; owner accepted (error bars of the AF mode unreliable at 250 bins;
  NVIDIA with the same seeds is consistent). A bug search is deferred (decision 41).
- GPU=none ctest (GNU): 53/53 candidate, 39/39 upstream master and merge-base.
- Accuracy (QREF, report only, from 207-2..207-5 / 302c-2): GPU at or below the correct CPU variants everywhere; NVIDIA beta = 10 pTT
  pooled -0.061 (default) / -0.121 (LOG); B580 beta = 4 -0.032 / -0.125; LOG beta = 10 band widens (row-permuted CPU variants lose
  accuracy); var8 / TAUSWAP=8 about 0.2 decades worse than host row order on LOG.
