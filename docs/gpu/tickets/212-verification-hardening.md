# 212 – Verification hardening (code)

- **Phase / lane:** 2 (NVIDIA lane; the env-variable gating also applies to `GPU=intel`)
- **Owner agent:** fortran-gpu-dev (Sonnet)
- **Review:** numerics-reviewer (code, acceptance file, scripts) · **Verification:** acceptance-verifier
- **Depends on:** the verification standard (merged 2026-10-29)
- **Base:** `gpu-port` d53f9779 (206 state + harness; the 207 τ code is **not** in this base)
- **Branch:** `gpu-port-212-hardening`
- **Acceptance file:** `testsuite/gpu_reference/acceptance/212.json`, pre-registered in its own commit before any
  implementation (the commit sha is recorded below under Approvals; DECISIONS.md is pinned at that commit and is not
  edited on this branch)
- **Status:** open (2026-10-29); scope confirmed by the owner's decision 2(a) of 2026-10-29

Must land before any further ticket closes (207 included).

## Goal
Debug and fault switches can no longer reach a production run unnoticed; every run states what it is; the default GPU
path becomes run-to-run reproducible (bitwise, same binary, same seeds) by fixing the reduction order; a graded fault
shows that the stabilization self-test resolves errors near its tolerance.

## Context
- Measured 2026-10-29 (gate report): the default path is **not** bit-reproducible run to run (same binary twice,
  12×12 β = 10: 20 of 38 datasets differ, ≤ 2.4e-11 relative; Precision Green mean 6.188 vs 6.180e-10); the LOG path is.
  cuBLAS/cuSOLVER are deterministic as configured (single stream, atomics off, default cuSOLVER mode). The remaining
  suspects are OpenMP device reductions whose combination order NVHPC does not specify: `k_colnorms` in
  `Libraries/Modules/alf_la_mod.F90` (feeds the `la_qrpp` pivot sort; near-ties flip pivots — default path only),
  `k_lognorm` (`Prog/gpu_stab_mod.F90`, STABLOG selection sort) and `k_prec` (precision sums).
- The 207 default-path paired spread (1.707 > 1.25) is suspected to come from the same pivot-order sensitivity. The
  paired test needs the 207 code; it is re-measured after the 207 rebase onto 212 (207 closing plan), not here.
- Variables the base code reads (grep of `gpu-port`, 2026-10-29): `ALF_GPU_DISABLE, _EXIT, _FAULT, _FAULT_STAB,
  _NOUPDATE, _NS, _SELFTEST, _STABCHECK, _STABMIX, _STABTIME, _TIMING, _TL, _TR, _TRACE, _UPDTEST, _VB`.

## Work (exact formats: the acceptance file depends on them)
1. **One reader** of every `ALF_GPU_*` variable at setup (a small module or routine in `gpu_state_mod.F90`).
   - A run is a **debug run** if `ALF_GPU_SELFTEST`, `ALF_GPU_UPDTEST` or `ALF_GPU_DEBUG=1` is set; otherwise it is a
     **production run**.
   - In a production run each of `ALF_GPU_FAULT, _FAULT_STAB, _NS, _STABMIX, _STABCHECK, _TL, _TR, _VB, _NOUPDATE,
     _TRACE, _TIMING, _STABTIME, _EXIT` is **refused**: the program prints exactly
     ` GPU debug env refused: <NAME> (production run; set ALF_GPU_DEBUG=1 for a debug run)` on stdout and in `info`,
     and stops with exit code 1 before the first sweep. `ALF_GPU_DISABLE` stays allowed (it switches the GPU path off,
     which the ACTIVE lines show).
   - Every GPU run prints at setup, on stdout and in `info`, exactly one line ` GPU debug env: none` (no `ALF_GPU_*`
     variable set) or ` GPU debug env: <NAME>=<value>, <NAME>=<value>` (sorted by name, every `ALF_GPU_*` variable that
     is set).
   - Every armed fault prints ` GPU fault armed: <NAME>=<value>` (stdout and `info`).
2. **Build stamp line** at setup: ` GPU build: git <hash from git.h> backend <CUDA|ONEMKL|HOST>`.
3. **Graded fault** for the stabilization self-test: `ALF_GPU_FAULT_STAB=eps:<x>` (debug runs only, as every fault)
   multiplies one element of the device G after the device CGR of the `ALF_GPU_SELFTEST=3` comparison by `(1 + x)`.
   It must make the self-test fail at `x = 1e-9` (line ` GPU self-test FAILED: stab: ...`, exit code 1) and pass at
   `x = 1e-13`.
4. **Deterministic reductions:**
   - `k_colnorms`, `k_lognorm`, `k_prec`: a fixed combination order (e.g. one thread per column with a serial sum, or a
     fixed-shape two-stage tree into scratch); no OpenMP `reduction` clause on floating-point sums or maxima that feed a
     decision or a reported number;
   - `la_init` (CUDA backend only): `cusolverDnSetDeterministicMode(…, CUSOLVER_DETERMINISTIC_RESULTS)` and
     `cublasSetAtomicsMode(…, CUBLAS_ATOMICS_NOT_ALLOWED)`; the HOST path of `alf_la_mod.F90` stays byte-identical in
     its preprocessed form.
5. **Micro-test** `t_determinism` in `testsuite/alf_la.tests` (added to its `CMakeLists.txt`): `la_qrpp` on fixed inputs
   including matrices with exactly equal column norms (duplicated and permuted columns), 1000 calls, pivots and factors
   compared bit for bit between calls; prints ` FAILED t_determinism` and stops with `error stop 1` on a mismatch
   (la_test_util convention). Test-only switch `ALF_LA_TEST_PLANT_ULP=1`: changes one input element by 1 ulp in one of
   the calls, so the test must fail (positive control). Under the HOST backend the test runs too (and passes).

## Files allowed
- `Prog/gpu_state_mod.F90`, `Prog/gpu_stab_mod.F90`, `Prog/gpu_wrap_mod.F90`, `Prog/gpu_update_mod.F90`,
  `Prog/Wrapgr_mod.F90` (only where it reads `ALF_GPU_*`)
- `Libraries/Modules/alf_la_mod.F90` (CUDA paths and `k_colnorms`; HOST paths byte-identical when preprocessed)
- `testsuite/alf_la.tests/t_determinism.F90` (new), `testsuite/alf_la.tests/CMakeLists.txt`
- `docs/gpu/scripts/nividia/t212_*` (development only)
- Not allowed: `testsuite/gpu_reference/acceptance/**` (frozen), any `GPU=none` behaviour, `configure.sh`.

## Acceptance
Pre-registered in `testsuite/gpu_reference/acceptance/212.json` (gates, thresholds, counts, controls; read it). In short:
- **Kernel self-tests** 204 (`SELFTEST=2`), 205 (`UPDTEST=1`), 206 (`SELFTEST=3`) PASS on the listed configs, default
  and LOG; the 206 fault controls and the graded fault behave as specified.
- **Refusal:** every listed variable alone in a production run is refused with the exact line and exit code 1; with
  `ALF_GPU_DEBUG=1` the same variable is accepted.
- **Production line:** every production GPU run prints ` GPU debug env: none` and the build stamp line.
- **Determinism:** the same binary twice gives bit-identical `data.h5` and full-precision Precision lines, default and
  LOG, 6×6 (200 one-sweep bins) and 12×12 β = 10 (8 bins × 20 sweeps).
- **Same seeds vs the CPU oracle:** 6×6 default and LOG, 12×12 default: no bin above 1e-4, overall and τ line ≤ 1e-9.
  (12×12 LOG is ticket 213.)
- **Precision mean** (Green, τ) ≤ 1.5× the larger CPU mean, 12×12 β = 4, default and LOG.
- **compare2** GPU vs GPU=none on fresh seeds, with a CPU-vs-CPU null control.
- **`t_determinism`** passes under CUDA; its planted-ulp control fails.
- **GPU=none:** preprocessed changed CPU-path Fortran identical (`alf_la_mod.F90`, 9 macro sets), ctest 53/53
  (52 + `t_determinism`).
- Report: cost of the deterministic reductions (s/sweep at N = 144 and 400, default and LOG), with stamped runs.

## Approvals
- Pre-registration commit: `6f106abf` (2026-10-29) on `gpu-port-212-hardening`; base `d53f9779`. Any change of
  `212.json` after this needs the owner's quoted approval and fresh seeds (VERIFICATION.md §3), unless it is made
  before any implementation result exists and is pre-registered again in a new commit recorded here.

## Evidence
`~/alf-gpu/work/t212/`: development runs (T0) anywhere below it; candidates in `cand/agent/<sha12>/` (T1) and
`cand/mgr/<sha12>/` (T2). Never delete raw outputs or binaries before the owner's gate; cleanup removes build trees only.

## Handoff notes
(implementer, at ~250k tokens or at a phase boundary)

## Out of scope
207 (τ pass) and its paired test; 213 (STABLOG same-seed excess); performance work beyond the deterministic reductions;
Intel/oneMKL.
