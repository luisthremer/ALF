# 207 – how it closes under the verification standard (updated with the owner's decisions of 2026-10-29)

207 is on hold (ticket 207 on branch `gpu-port-207-wip`). The existing `mgr207*` and agent results are
**expectations only**, not acceptance. Decided at the gate of 2026-10-29 (DECISIONS.md):
- acceptance 3: the rewrite approved on 2026-10-29 stays in force, **including the sd condition** (decision 1);
- the default paired FAIL (1.707 > 1.25): find the cause first in ticket 212 (decision 2a); until then 207 stays open;
- same-seed bound 1e-9; the LOG 12×12 excess (1.214e-9) is an open finding of the 206 code (decision 3);
- acceptance 6 (Intel) is BLOCKED, never PASS (decision 4); the verdict then reads BLOCKED, not PASS;
- budget: T1 ≈ 90, T2 ≈ 120 GPU-min (decision 11).

- the LOG 12×12 same-seed excess (206 code) has its own ticket 213; 207 is not held open by it (decision 3b).

- 207 does not close while acceptance 6 is BLOCKED without the owner's explicit OK at that time (decision 4b); a
  BLOCKED verdict is reported to the owner, never committed as accepted.

Order:
1. **212 lands** (hardening, cause of the default-path spread, deterministic reductions). 212 changes code 207 shares.
2. **Rebase** the 207 code on the accepted 212 commit. Re-derive in `207.draft.json`: `base_sha`, the base builds,
   the GPU=none counts, the debug-variable switch that 212 introduces (control runs need it), `parameters_sha256`.
3. **Pre-register** `testsuite/gpu_reference/acceptance/207.json` in its own commit (recorded in the ticket and
   DECISIONS.md); review part (b) of `numerics-reviewer` before any candidate result is seen.
4. **Freeze the candidate** (a later commit on the 207 branch); the implementer runs T1 (`--who agent`, dev seeds).
5. **Positive controls first** inside T2: `ALF_GPU_TAUSWAP=8` (LOG) and `ALF_GPU_TAU_NS=0` (default) through the
   paired statistic; if they are not killed, the paired test cannot see the defect it was built for → owner.
6. **Full T2** (`acc_ship.sh … mgr …`, `acc_driver.py --who mgr`): deterministic gates in full, statistical gates once
   on fresh seeds, all controls, GPU=none.
7. **Second review** (numerics-reviewer on Opus): (a) code `base..candidate` incl. the Opus diagnostic and the review
   fixes, (c) any script change; **acceptance-verifier** recomputes every verdict.
8. **Commit** (`acc_check_commit.sh`: only `docs/**` outside `docs/gpu/scripts/` may differ from the accepted sha),
   push to `fork`; archive the evidence (`~/alf-gpu/archive/t207/`).
