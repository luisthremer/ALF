# 207 closing state (single source of truth; read only this file + CLAUDE.md)

Decision 41 (DECISIONS.md): 207 closes on functional checks only. No agents, no skills, one step per turn,
long jobs detached on the remote machine, then end the turn.

## Current step
Builds (check 1) running detached on both machines. Source = git archive of gpu-port-207-final (9c35b3b5) in
~/alf-gpu/work/t207-final/src; scripts docs/gpu/scripts/t207_final/{build,builds_nvidia,builds_intel}.sh.

## Next step
When logs/builds.txt ends with BUILDS_DONE and four `built` lines: T0 self-test (ALF_GPU_SELFTEST=4) and same-seed 6x6
vs GPU=none (nv*_n / inn_*), both paths. Then production-size observables runs; GPU=none ctest (GNU, nividia).

## Running
nividia: scripts/builds_nvidia.sh -> nv_def, nv_log, nvn_def, nvn_log (log ~/alf-gpu/work/t207-final/logs/builds.txt, make_*.log)
intel:   scripts/builds_intel.sh  -> in_def, in_log, inn_def, inn_log (same paths); 302c-2d queue also runs there (J=6)

## Functional checks (per GPU, default and LOG)
Scripts: docs/gpu/scripts/nividia/t207_3/t207_3_t0.sh; docs/gpu/scripts/intel/t302c/t302c2c_t0.sh + t302c2b_t0eval.py;
testsuite/gpu_reference/{prep,compare,compare2,summarize}.py.
1. Build:            nividia -, intel -
2. Self-tests (=4):  nividia -, intel -
3. Same-seed ≤1e-9:  nividia - (def/LOG), intel - (def/LOG)
4. Observables χ²:   nividia - (def/LOG), intel - (def/LOG)
5. GPU=none ctest = upstream (GNU, nividia): -

## Accuracy report (report only, from tickets; fill at close)
GPU at/below correct CPU variants everywhere; NVIDIA β=10 pTT pooled −0.061 def / −0.121 LOG; B580 β=4 −0.032 / −0.125;
LOG β=10 band widens (row-permuted CPU variants lose accuracy); var8/TAUSWAP=8 ≈0.2 decades worse than host order on LOG.

## Open problems
None.
