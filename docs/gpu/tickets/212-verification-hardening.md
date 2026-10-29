# 212 – Verification hardening (code)

- **Phase / lane:** 2 (NVIDIA lane; the env-variable gating also applies to `GPU=intel`)
- **Owner agent:** fortran-gpu-dev (Sonnet)
- **Review:** numerics-reviewer (code, acceptance file, scripts) · **Verification:** acceptance-verifier
- **Depends on:** the owner's approval of `gpu-port-process-lanes` (harness)
- **Status:** draft (2026-10-29), not started

Must land before any further ticket closes (207 included).

## Goal
Debug and fault switches can no longer reach a production run unnoticed; every run states what it is; the default GPU
path becomes run-to-run reproducible if the owner decides so; two cheap gates guard determinism and non-interference.

## Context
- Today these are read without a self-test gate and change production numerics (measured or read in the code,
  2026-10-29): `ALF_GPU_TAUSWAP` (`gpu_tau_mod.F90`, diag_setup), `ALF_GPU_TAU_NS` (`gpu_tau_mod.F90:178`),
  `ALF_GPU_NS` (`gpu_stab_mod.F90:156`), `ALF_GPU_FAULT_STAB` (`gpu_stab_mod.F90:1401`, used at `:456, :511, :798`),
  `ALF_GPU_STABMIX` (`:1388–1391`). `ALF_GPU_FAULT_TAU` is gated since 207's review fix.
- The default path is not bitwise reproducible run to run (same binary twice, 12×12 β = 10: ≤ 2.4e-11 relative in
  observables); LOG is. Suspects: OpenMP device reductions whose combination order NVHPC does not specify —
  `k_colnorms` (feeds the `la_qrpp` pivot sort; ties flip pivots), `k_lognorm`, `k_prec`
  (`alf_la_mod.F90` ≈ 1014, `gpu_stab_mod.F90:1281, 1314`). cuBLAS/cuSOLVER are deterministic as configured.
- Gate report of 2026-10-29 (determinism test 3a) for the measured state; DISCREPANCIES.md row of 2026-10-29.

## Files allowed
- `Prog/gpu_state_mod.F90`, `Prog/gpu_stab_mod.F90`, `Prog/gpu_tau_mod.F90`, `Prog/gpu_update_mod.F90`,
  `Prog/gpu_wrap_mod.F90` (env handling, info lines, reductions)
- `Libraries/Modules/alf_la_mod.F90` (fixed-order `k_colnorms`; `la_init` library modes) — HOST path unchanged
- `testsuite/alf_la.tests/` (reduction determinism micro-tests)
- `testsuite/gpu_reference/acceptance/212.json`, `docs/gpu/scripts/nividia/t212_*`

## Work
1. One module-level routine reads every `ALF_GPU_*` variable once at setup. A production run (no `ALF_GPU_SELFTEST`,
   `ALF_GPU_UPDTEST` or the new `ALF_GPU_DEBUG=1`) refuses to start if any debug/fault variable is set, with a clear
   message. Self-test runs accept them as today.
2. Info lines at setup: "GPU debug env: none" or the list of variables read with non-default values; "fault X armed"
   whenever a fault is active; a provenance header (git hash from `git.h`, compile flags, backend).
3. `ALF_GPU_FAULT_TAU=eps:<x>` (graded fault: one element of GT0 at one slice × (1 + x)); `ALF_GPU_TAUDIAG_SEED`.
4. **Only if the owner decides (gate of 2026-10-29):** fixed-order reductions in `k_colnorms`, `k_lognorm`, `k_prec`
   (one thread per column, serial sum; or a fixed-shape two-stage tree into scratch); `cusolverDnSetDeterministicMode`
   and `cublasSetAtomicsMode(NOT_ALLOWED)` in `la_init`; `CUBLAS_WORKSPACE_CONFIG=:4096:8` in the run scripts.
5. Gates: (a) run-twice-identical — the same binary twice, same seeds, all `data.h5` datasets and full-precision
   `info` Precision lines bit-identical (6×6 fd and 12×12 β = 10, default and LOG); (b) the τ pass leaves the sweep
   unchanged — `Ltau = 1` vs the previous accepted binary with identical equal-time/scalar datasets.

## Acceptance (to be pre-registered in `212.json` before implementation)
- Production refusal: each listed variable set alone in a production run → the run refuses (exit ≠ 0, message names
  it); positive control: the same variable with `ALF_GPU_SELFTEST` → accepted.
- Graded fault: `eps:1e-9` FAILS the 207 self-test, `eps:1e-13` PASSES.
- Determinism (if item 4 is decided): gate (a) bit-identical; micro-tests of the reductions (1000 calls, tie cases)
  bit-identical; all 204–207 kernel self-tests and fault controls unchanged; same-seed vs CPU within the 207 bounds.
- `GPU=none`: preprocessed sources of changed files identical (`acc_gpunone.sh`), ctest count unchanged.
- Cost of item 4 reported (s/sweep, N = 144 and 400).

## Out of scope
New numerics; observables; anything in 207's open decisions.
