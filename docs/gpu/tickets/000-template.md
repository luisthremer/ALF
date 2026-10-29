# NNN – short title

- **Phase / lane:** 2 (NVIDIA) | 3 (Intel) | …
- **Owner agent:** fortran-gpu-dev (Sonnet; Opus only on the manager's proposal) | test-runner | manager
- **Review:** numerics-reviewer (code, acceptance file, scripts) · **Verification:** acceptance-verifier
- **Depends on:** –
- **Status:** open | in progress | candidate frozen | review | on hold | done

## Goal
What should be true afterwards, in one or two sentences.

## Context
What the agent needs to know (files with line refs, decisions, prior measurements). Self-contained.

## Files allowed
- `path/to/file.F90`
Anything else: stop and report. Shared files of both lanes (`gpu_stab_mod`, `gpu_state_mod`, `gpu_wrap_mod`,
`PLAN.md`) only where this ticket owns the edit (PLAN.md, file ownership).

## Acceptance
Pre-registered in `testsuite/gpu_reference/acceptance/NNN.json`, committed with this ticket **before** any
implementation result is seen, and frozen with the candidate (`docs/gpu/VERIFICATION.md`). For each gate:
- oracle; statistic; threshold and its reason; the gated row/dataset; expected counts (configs, datasets, events, tests);
- tier (T0 every iteration / T1 pre-candidate / T2 manager acceptance);
- at least one positive control with its expected outcome (FAIL);
- statistical gates: α, multiple-comparison correction, the detectable effect size, the seed sets (dev seeds committed;
  acceptance seeds drawn fresh by the driver at T2).
Also: build profiles and machines; the `GPU=none` guarantee (preprocessed-source identity of changed files, `data.h5`
identity, ctest count); the report (diff summary, transfers, timings).

## Approvals
Owner approvals that change this ticket's criteria, quoted and dated (copy in `docs/gpu/DECISIONS.md`). Empty at start.

## Evidence
- Work dir `~/alf-gpu/work/tNNN/`; candidates in `cand/<sha12>/` (from `acc_ship.sh`), runs in `cand/<sha12>/runs/<who>/`.
- Never delete raw outputs or binaries of this ticket before the owner's gate; cleanup removes build trees only.
- Every log carries the provenance stamp; results without a stamp matching the candidate are not evidence.

## Handoff notes
Written by the implementer at ~250k tokens or at a phase boundary: state, decisions, what was tried, open problems,
evidence paths. A fresh agent continues from here.

## Out of scope
- …
