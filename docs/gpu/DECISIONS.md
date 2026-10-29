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

## Pending at the gate of 2026-10-29

See the gate report: the 207 acceptance-3 rule, the default paired FAIL, Intel acceptance 6 as BLOCKED, the process
diff (incl. max → mean/quantile), the oneAPI fix, MIG, the old t104–t108 directories, `module avail` on SuperMUC,
the agent setup, the GPU-time budget per ticket, the lanes' go-ahead.
