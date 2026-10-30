# Owner decisions

Every owner decision that changes scope, criteria, process or rules, quoted verbatim where it changes a criterion,
with date and channel. Tickets link here; a decision that exists only in a chat does not count.

| Date | Channel | Decision (verbatim where it changes a criterion) | Applies to |
|---|---|---|---|
| 2026-10-28 | manager session | Go-ahead for ticket 207, then stop before 208. | 207 |
| 2026-10-28 | manager session | The `main.F90` hooks of 207 outside the three `Tau_m` call sites (setup, self-test hook, finalize; pure additions) are accepted, pending review. | 207 |
| 2026-10-28 | manager session | The STABLOG precision-tail investigation goes to a fresh Opus implementer; the Sonnet implementer finishes its runs. | 207 |
| 2026-10-29 | manager session | "Yes: for Precision tau (and the 206 Green max), replace the 1.5x-max rule with the paired same-trajectory test as you proposed; report the max as information. Keep the diagnostic in the code. Then review, re-run and commit 207." | 207 acceptance 3, 206 Green max |
| 2026-10-29 | owner instruction "Two lanes, a stricter verification standard, faster validation" | "do not drop the sd condition, and do not commit 207." The default-path paired spread (sd ratio 1.707 > 1.25) stays an open FAIL until the gate. Options to be put neutrally at the gate: (a) find the cause first; (b) a new, pre-registered criterion, calibrated with null and known-bias runs, only after owner approval; (c) keep 207 open. "From now on, when a criterion fails, never put relaxing it as the first or recommended option." | 207; all tickets |
| 2026-10-29 | same instruction | 207 is on hold ("on hold: awaiting owner decision and verification standard"). No raw output, run directory or binary of an unclosed ticket is deleted. The process change goes to its own branch and is not merged before the owner approves. No lane starts work and no ticket agent launches until a new go-ahead. | 207; process |
| 2026-10-29 | same instruction | Approved within limits: an NVIDIA MPS trial as an ordinary user on `nividia` (own pipe/log directories, no root, no compute-mode change); confinement exceptions for `~/alf-gpu/backup/`, `~/alf-gpu/archive/`, `~/alf-gpu/refcache/`, `~/alf-gpu/work/mps_trial/`. | machines |

## Gate of 2026-10-29: decisions

The owner answered the manager's recommendations with "ok" (2026-10-29, manager session), adopting all of them:

| # | Topic | Decision |
|---|---|---|
| 1 | 207 acceptance-3 criterion | Decided after item 2(a). Until then the rewrite approved on 2026-10-29, **including the sd condition**, stays in force. A new criterion only if the cause is understood and the condition still fails: pre-registered, calibrated beforehand with a null run (CPU vs 1-ulp-perturbed CPU) and a known defect (`ALF_GPU_TAUSWAP=8`), then judged once on fresh seeds. |
| 2 | Default paired FAIL (1.707 > 1.25) | (a) find the cause first (`k_colnorms` reduction / `la_qrpp` pivot order, ticket 212). If a fixed-order fix brings the ratio under 1.25 on fresh seeds, 207 passes without a criterion change. Until then (c): 207 stays open. |
| 3 | Same-seed bound | 1e-9 stays. The LOG 12×12 excess (1.214e-9 from bin 184, equal-time only) is recorded as an open finding of the 206 code; any other bound only pre-registered and calibrated against the CPU-vs-CPU same-seed spread. |
| 4 | 207 acceptance 6 (Intel) | Recorded as **BLOCKED** (never PASS); the nividia ifx compile/link is a stand-in only; re-run after the oneAPI fix. |
| 5 | Process diff `gpu-port-process-lanes` | Approved **after a second independent review** of its final state. MPS rule (≤ 2 clients, started/stopped by the manager or the driver only) and CPU-alongside rule (≤ 3 single-threaded runs) approved. Max-error statistics: information only; acceptance on the mean. The quantile option only once calibrated. |
| 6 | oneAPI on `intel` | Option (a): the admin installs `intel-oneapi-mkl-2025.3` and holds the ifx and MKL packages (`apt-mark hold`). Owner action. |
| 7 | MIG | Not now. |
| 8 | t104–t108 on nividia | Archive logs and summaries to `~/alf-gpu/archive/`, then delete raw runs and binaries **after the owner's OK on the listed contents**. |
| 9 | `module avail` on SuperMUC | The owner runs it; it decides the toolchain pin for both Intel machines. |
| 10 | Agent setup | Opus `acceptance-verifier`, Opus `numerics-reviewer`, Sonnet implementers (Opus per numerically delicate ticket on the manager's proposal), respawn and handoff rules as in VERIFICATION.md §10b. |
| 3b | Same-seed LOG 12×12 excess (1.214e-9, 206 code) | (2026-10-29, after the second review) "Own ticket for 206": the excess gets its own ticket (213); 207 is not held open by it. The pair stays in 207 as the report-only gate `same_seed_log12_reported`; the 1e-9 bound is unchanged. |
| 4b | Closing 207 with acceptance 6 BLOCKED | (2026-10-29) "Not without my OK": 207 may not close while acceptance 6 is BLOCKED unless the owner explicitly approves it at that time. |
| 11 | GPU budget and lanes | Standing budget per ticket: T1 ≈ 90 GPU-min, T2 ≈ 120 GPU-min unattended, every run ≤ 15 min. Lanes go ahead once the process diff is approved and merged. Order: NVIDIA lane 212 → 207 (controls, T2 on fresh seeds) → 208; Intel lane 302a once oneAPI is fixed. |
| 12 | compare2 positive-control shift, ticket 212 | (2026-10-29, manager session) "Approve shift 8": the compare2 positive control of 212 plants a shift of 8 standard errors instead of the 5 of VERIFICATION.md §5 (at m = 11, α = 0.01 Šidák a shift of 5 is caught only ~60 % of the time, 8 ~99 %). The compare2 criterion itself is unchanged. |
