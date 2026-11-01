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
1. Run type "QREF diagnostic": the env allow-list of `t207_d5_qparse.py` (`env_role`): main runs exactly `DEBUG=1`,
   `TAUDIAG_QREF`, `_QHOST=1`, `_QVAR=8`, `_FROM`; control runs that set plus `TAU_NS=0` / `TAUSWAP=8` / the 207-3 graded
   fault, in declared control runs only. Stride 3 (decision 35).
2. Build the offline program and the GNU oracle objects from the candidate's `git archive`, with provenance.
3. Process each run's dumps before the next run starts (nividia has ~11 GB free; ~1.6 GB per run): verdict inputs
   computed first, then retention as in 207-2 (manifest, full `qref.txt`, kept dumps by the gcd keep rule, KEEP 39).
   Reuse the 207-2 scripts (`t207_d5_proc.sh`, `t207_d5_qparse.py`, `t207_d5_agg.py`, `t207_d5_gate.py`) rather than
   re-implementing them; the harness calls them and records their hashes.
3a. Both machines (decision 36): QREF runs on `intel` (B580, `GPU=intel`, after 302c-1) with the offline program built
   there (gfortran 15.2 present on `intel`, real(16)); the Intel reference ensemble has an MKL/ifx host member in place
   of `nvh` (name fixed with 302c-2); keep the member list a parameter of the acceptance file.
4. Gate type computing S_ens(pTT) and the pointwise statistic from `qref.txt` with the frozen ensemble list from the
   acceptance file; input checks (qself, qbr ≤ 1e-25; n0 identical to cpu; nvh present; device branch = host branch
   at every point; point count) → INPUT_ERROR.
5. Identity check: `data.h5` with QREF on and off, same seeds, identical.
6. `verdict.py --selfcheck` mutations: NaN, missing member, missing point, shifted fields, device and var8 swapped.

## Files allowed
`testsuite/gpu_reference/acceptance/` (scripts, README), `docs/gpu/VERIFICATION.md` (gate description).

## Done when
Self-checks pass, an end-to-end dry run on a dev build passes, review (c) approves, merged into `gpu-port`.
