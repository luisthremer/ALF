# 207-5 – acceptance harness support for the QREF accuracy gates

- **Phase / lane:** 2 (NVIDIA) · part of 207 (decisions 33, 34)
- **Owner agent:** fortran-gpu-dev (Sonnet); manager integrates
- **Review:** numerics-reviewer (scripts, part c) · merged into `gpu-port` before 207-6 (the judging tools must be in
  207's base)
- **Depends on:** 207-2 (output format and parser)
- **Branch:** `gpu-port-207-harness`
- **Status:** open

## Goal
`acc_driver.py` / `verdict.py` can run and judge the QREF gates from a pre-registered acceptance file (reviewer
change 6, retention of decision 33).

## Tasks
1. Run type "QREF diagnostic": allowed debug variables only `ALF_GPU_TAUDIAG_QREF`, `_QHOST`, `_QVAR=8`, `_FROM`
   (and the 207-3 fault in declared control runs).
2. Build the offline program and the GNU oracle objects from the candidate's `git archive`, with provenance.
3. Process each run's dumps before the next run starts (nividia has ~13 GB free; ~1.6 GB per run): verdict inputs
   computed first, then retention as in 207-2 (manifest, full `qref.txt`, every 40th dump).
4. Gate type computing S_ens(pTT) and the pointwise statistic from `qref.txt` with the frozen ensemble list from the
   acceptance file; input checks (qself, qbr ≤ 1e-25; n0 identical to cpu; nvh present; device branch = host branch
   at every point; point count) → INPUT_ERROR.
5. Identity check: `data.h5` with QREF on and off, same seeds, identical.
6. `verdict.py --selfcheck` mutations: NaN, missing member, missing point, shifted fields, device and var8 swapped.

## Files allowed
`testsuite/gpu_reference/acceptance/` (scripts, README), `docs/gpu/VERIFICATION.md` (gate description).

## Done when
Self-checks pass, an end-to-end dry run on a dev build passes, review (c) approves, merged into `gpu-port`.
