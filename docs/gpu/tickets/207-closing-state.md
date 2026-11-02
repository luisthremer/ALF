# 207 closing state (single source of truth; read only this file + CLAUDE.md)

Decision 41 (DECISIONS.md): 207 closes on functional checks only. No agents, no skills, one step per turn,
long jobs detached on the remote machine, then end the turn.

## Current step
Check 1 PASS both. T0: LOG self-test PASSED on both (B580 and NVIDIA). Default-path self-test needs a -DQR_PREPIVOT build
(refusal printed by the code), so nv_def_pp / in_def_pp are building (logs/build_pp.txt). Same-seed (t0.sh) still running on nividia
(bash -x chain, logs/t0.txt, dbg in logs/dbg.txt); intel chain results in logs/t0.txt were from a mixed-up attempt: redo intel.
Lesson: never overwrite a running script; never pkill -f without brackets; don't run scripts in the foreground.

## Next step
1. When *_pp built: run default self-test (ALF_GPU_DEBUG=1 ALF_GPU_EXIT=1 ALF_GPU_SELFTEST=4) with nv_def_pp / in_def_pp (t0.sh needs the
   st_def line to use ${P}_def_pp; edit only when no chain is running).
2. intel: wipe runs/t0_*, rerun same-seed + selftest cleanly. 3. check 4 and 5.

## Running
nividia: t0.sh chain (same-seed etc.); nv_def_pp build. intel: in_def_pp build.

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
