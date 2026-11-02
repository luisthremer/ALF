# 207 closing state (single source of truth; read only this file + CLAUDE.md)

Decision 41 (DECISIONS.md): 207 closes on functional checks only. No agents, no skills, one step per turn,
long jobs detached on the remote machine, then end the turn.

## Current step
DONE: candidate branch `gpu-port-207-final` built (no conflicts; merged 207-3, 207-r, gpu-port d6598bed onto 302c-2).
Was: build candidate branch `gpu-port-207-final` (from gpu-port-302c-2 4b625b52; merge gpu-port-207-3 c48e61e2,
gpu-port-207-r 6f452aa8, gpu-port-302b-impl2/gpu-port). Keep 302c-1's removal of the GPU=intel guard block.

## Next step
Start builds + T0 (self-test, same-seed) on both machines: on both machines (nividia: GPU=nvidia nvhpc; intel: GPU=intel ifx+oneMKL), work dir ~/alf-gpu/work/t207-final/.

## Running
Nothing.

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
