# 207 – how it would close under the verification standard (draft, 2026-10-29)

207 is on hold (ticket 207 on branch `gpu-port-207-wip`). The existing `mgr207*` and agent results are
**expectations only**, not acceptance. Order:

1. **Owner decisions** (gate of 2026-10-29), each recorded in DECISIONS.md and the acceptance file:
   - the acceptance-3 criterion: the original max rule, the rewritten rule (approved 2026-10-29), or a new
     pre-registered criterion calibrated with null and known-bias runs, applied to fresh seeds;
   - the default-path paired FAIL (sd ratio 1.707 > 1.25): (a) find the cause first (`k_colnorms` pivot order in
     `la_qrpp`; a fixed-order fix in 212 may make it pass without a criterion change), (b) a new pre-registered
     criterion, (c) keep 207 open;
   - Intel acceptance 6 recorded as BLOCKED (oneAPI version mix) or not;
   - the GPU-time budget for T2 (~2 h unattended).
2. **Freeze the candidate:** a commit on `gpu-port-207-wip` (or its successor) that contains the code, the ticket and
   `testsuite/gpu_reference/acceptance/207.json` (from `207.draft.json` and the decisions). Record the sha.
3. **212 and the harness have landed** (212 changes code that 207 shares, so the 207 candidate is rebased on it).
4. **Positive controls first:** `ALF_GPU_TAUSWAP=8` (LOG) and `ALF_GPU_TAU_NS=0` (default) through the paired
   statistic on dev seeds; if they are not killed, the paired test cannot see the defect it was built for → owner.
5. **Full T2** (`acc_driver.py --who mgr`): deterministic gates in full, statistical gates once on fresh seeds, all
   controls; `acc_gpunone.sh`.
6. **Second review** (numerics-reviewer on Opus): (a) the code `64ecde72..candidate` including the Opus diagnostic and
   the review fixes, (b) `207.json`, (c) any script change. **acceptance-verifier** recomputes every verdict.
7. **Commit** (`acc_check_commit.sh`: only `docs/**` may differ from the accepted sha), push to `fork`; archive the
   evidence (`~/alf-gpu/archive/t207/`).
