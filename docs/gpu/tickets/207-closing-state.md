# 207 closing state (single source of truth; read only this file + CLAUDE.md)

Decision 41 (DECISIONS.md): 207 closes on functional checks only. No agents, no skills, one step per turn,
long jobs detached on the remote machine, then end the turn.

## Current step
Check 1 (builds) PASSED on both machines (nv_def/nv_log/nvn_*; in_def/in_log/inn_*; logs/builds.txt). T0 (checks 2+3) started
detached on both: scripts/t0.sh (docs/gpu/scripts/t207_final/t0.sh) -> ~/alf-gpu/work/t207-final/logs/t0.txt (ends T0_DONE),
runs/t0_*/. Self-test uses the production builds (no -DQR_PREPIVOT): if it FAILS only because of that, rebuild *_pp.

## Next step
Read logs/t0.txt on both. Then check 4 (production-size observables vs CPU reference, compare.py/compare2.py) and check 5
(GNU ctest on nividia).

## Running
nividia and intel: t0.sh (self-test + same-seed, def and LOG).

## Functional checks (per GPU, default and LOG)
Scripts: docs/gpu/scripts/nividia/t207_3/t207_3_t0.sh; docs/gpu/scripts/intel/t302c/t302c2c_t0.sh + t302c2b_t0eval.py;
testsuite/gpu_reference/{prep,compare,compare2,summarize}.py.
1. Build:            nividia PASS (4 built), intel PASS (4 built)
2. Self-tests (=4):  nividia -, intel -
3. Same-seed ≤1e-9:  nividia - (def/LOG), intel - (def/LOG)
4. Observables χ²:   nividia - (def/LOG), intel - (def/LOG)
5. GPU=none ctest = upstream (GNU, nividia): -

## Accuracy report (report only, from tickets; fill at close)
GPU at/below correct CPU variants everywhere; NVIDIA β=10 pTT pooled −0.061 def / −0.121 LOG; B580 β=4 −0.032 / −0.125;
LOG β=10 band widens (row-permuted CPU variants lose accuracy); var8/TAUSWAP=8 ≈0.2 decades worse than host order on LOG.

## Open problems
None.
