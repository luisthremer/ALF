# 207 closing state (single source of truth; read only this file + CLAUDE.md)

Decision 41 (DECISIONS.md): 207 closes on functional checks only. No agents, no skills, one step per turn,
long jobs detached on the remote machine, then end the turn.

## Current step
Check 5: candidate GNU devel ctest 53/53 pass; upstream master 39/39 pass; MISMATCH IN COUNT is explained by the diff
(candidate adds 14 alf_la.tests LA-layer tests from 302a, lacks MALA.tests because its base f81b5cbc predates PR #592).
Baseline ctest of the port's merge-base f81b5cbc running (scripts/ctest_base.sh -> logs/ctest_base.txt) to show 39 - MALA + 14 = 53.
Check 4 running detached on both: scripts/obs.sh (hub12x12_cb beta 4, 250x70, Ltau 1, seeds rotated 37; GPU def then LOG; then
compare2.py vs CPU ref ~/alf-gpu/runs/ref_gnu-{DEFAULT,LOG}/hub12x12_cb, copied to intel work/t207-final/ref) -> logs/obs.txt (ends OBS_DONE).

## Next step
Read logs/obs.txt on both (compare2 verdict + chi2), logs/ctest_base.txt; fill checks 4/5; then merge to gpu-port, push, close tickets 207 and 207-7.

## Running
nividia: ctest_base.sh, obs.sh. intel: obs.sh.

## Functional checks (per GPU, default and LOG)
1. Build:            nividia PASS (4+1 built), intel PASS (4+1 built)
2. Self-tests (=4):  nividia PASS def+LOG; intel PASS def+LOG (hub6x6_cb, work/t207-final/runs/t0_st_*)
3. Same-seed ≤1e-9 (hub6x6_cb_fd, 60 bins, max dev): nividia def 1.70e-10, LOG 2.54e-10; intel def 2.32e-10, LOG 2.90e-10 -> PASS
4. Observables χ²:   nividia - (def/LOG), intel - (def/LOG)
5. GPU=none ctest (GNU, nividia): candidate 53/53 pass, upstream master 39/39 pass (count differs by design, see above)

## Accuracy report (report only, from tickets; fill at close)
GPU at/below correct CPU variants everywhere; NVIDIA β=10 pTT pooled −0.061 def / −0.121 LOG; B580 β=4 −0.032 / −0.125;
LOG β=10 band widens (row-permuted CPU variants lose accuracy); var8/TAUSWAP=8 ≈0.2 decades worse than host order on LOG.

## Open problems
None.
