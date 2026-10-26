# NNN – short title

- **Phase:** 0 | 1 | 2 | 3
- **Owner agent:** fortran-gpu-dev | test-runner | manager
- **Depends on:** –
- **Status:** open | in progress | review | done

## Goal
What should be true afterwards, in one or two sentences.

## Context
What the agent needs to know (files with line refs, decisions, prior measurements). Self-contained.

## Files allowed
- `path/to/file.F90`
Anything else: stop and report.

## Acceptance
- Build: which `configure.sh` profile(s), which machine.
- Tests: exact commands, expected result, tolerance (e.g. rel. diff of G ≤ 1e-10).
- CPU path: `GPU=none` unit tests still pass (state which).
- Report: diff summary + measured numbers.

## Out of scope
- …
