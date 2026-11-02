# 207 closing state (single source of truth; read only this file + CLAUDE.md)

Decision 41 (DECISIONS.md): 207 closes on functional checks only. No agents, no skills, one step per turn,
long jobs detached on the remote machine, then end the turn.

## Current step
All five checks measured. Check 4 FAILS on Intel as measured (see below); NVIDIA is consistent. STOPPED for the owner (decision 41:
record and tell, no debugging campaign). Nothing merged into gpu-port, 207 NOT marked closed.

## Next step (owner's call)
(a) accept as is (compare2 error bars are under-resolved on every observable, the 3 failing ones are one strongly correlated AF-mode
group), or (b) one more Intel seed rotation (~7 min on B580, scripts/obs.sh with another offset) to see whether the 6 sigma is a fluke.
Then merge gpu-port-207-final -> gpu-port (ff main checkout as in the task), push fork, mark 207 and 207-7 closed.

## Running
Nothing.

## Functional checks (per GPU, default and LOG)
1. Build:            nividia PASS (4+1 built), intel PASS (4+1 built)
2. Self-tests (=4):  nividia PASS def+LOG; intel PASS def+LOG (hub6x6_cb, work/t207-final/runs/t0_st_*)
3. Same-seed ≤1e-9 (hub6x6_cb_fd, 60 bins, max dev): nividia def 1.70e-10, LOG 2.54e-10; intel def 2.32e-10, LOG 2.90e-10 -> PASS
4. Observables (hub12x12_cb beta 4, 250x70, rotated seeds, GPU vs CPU ref 16-rank GNU, compare2.py; runs/obs_{def,log}, logs/obs.txt):
   nividia def and LOG: sum t^2 = 9.2 vs chi2_11 (p = 0.60), no observable DIFF; compare2 verdict INCONCLUSIVE (all 11 flagged
   UNDER-RESOLVED, N=14 bins, the known plateau limitation); full-array chi2/n 0.3-1.4, max|z| <= 4.05. Wall 791 s / 920 s.
   intel def and LOG: sum t^2 = 141.7 (p = 0.000). DIFF (6.3 sigma): Szz(pi,pi) 2.33+-0.17 vs ref 3.86+-0.17, Szz(pi,pi,tau=b/2)
   1.29+-0.17 vs 2.86+-0.17, chi_zz(pi,pi) 5.98+-0.74 vs 12.33+-0.68 (one correlated AF-mode group); the other 8 observables ok
   (E, Kin, Pot, Den, G, Szz(pi/2,pi/2), tau G within 2.6 sigma). Wall 430 s / 418 s. NVIDIA with the same seeds gives
   Szz(pi,pi) 3.74+-0.70 (error bar 4x the Intel one): AF-mode error bars are not trustworthy at 250 bins, so fluke vs bug is open.
   def and LOG numbers are identical to all printed digits on both GPUs (same seeds, same-seed deviation ~1e-10).
5. GPU=none ctest (GNU, nividia): candidate 53/53 pass; upstream master 39/39; merge-base f81b5cbc 39/39. Candidate = base 39 + 14
   alf_la tests (302a); no MALA.tests (base predates PR 592). Pass rate 100% everywhere.

## Accuracy report (report only, from tickets; fill at close)
GPU at/below correct CPU variants everywhere; NVIDIA β=10 pTT pooled −0.061 def / −0.121 LOG; B580 β=4 −0.032 / −0.125;
LOG β=10 band widens (row-permuted CPU variants lose accuracy); var8/TAUSWAP=8 ≈0.2 decades worse than host order on LOG.

## Open problems
None.
