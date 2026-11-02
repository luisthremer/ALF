# 207 closing state (single source of truth; read only this file + CLAUDE.md)

Decision 41 (DECISIONS.md): 207 closes on functional checks only. No agents, no skills, one step per turn,
long jobs detached on the remote machine, then end the turn.

## Current step
Checks 1, 2, 3 DONE on both GPUs. Check 5 (GNU ctest, candidate vs upstream master, nividia) running detached:
scripts/ctest.sh -> ~/alf-gpu/work/t207-final/logs/ctest_cand.txt, ctest_up.txt (ends CTEST_DONE).

## Next step
Check 4: production-size observables per GPU and path vs CPU reference (compare.py / compare2.py, report chi2). Need to find
the reference data (testsuite/gpu_reference, refcache on nividia) and a production-size config; read narrowly.

## Running
nividia: ctest.sh (4 cores). intel: nothing.

## Functional checks (per GPU, default and LOG)
1. Build:            nividia PASS (4+1 built), intel PASS (4+1 built)
2. Self-tests (=4):  nividia PASS def+LOG; intel PASS def+LOG (hub6x6_cb, work/t207-final/runs/t0_st_*)
3. Same-seed ≤1e-9 (hub6x6_cb_fd, 60 bins, max dev): nividia def 1.70e-10, LOG 2.54e-10; intel def 2.32e-10, LOG 2.90e-10 -> PASS
4. Observables χ²:   nividia - (def/LOG), intel - (def/LOG)
5. GPU=none ctest = upstream (GNU, nividia): -

## Accuracy report (report only, from tickets; fill at close)
GPU at/below correct CPU variants everywhere; NVIDIA β=10 pTT pooled −0.061 def / −0.121 LOG; B580 β=4 −0.032 / −0.125;
LOG β=10 band widens (row-permuted CPU variants lose accuracy); var8/TAUSWAP=8 ≈0.2 decades worse than host order on LOG.

## Open problems
None.
