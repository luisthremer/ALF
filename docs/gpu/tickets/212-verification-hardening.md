# 212 – Verification hardening (code)

- **Phase / lane:** 2 (NVIDIA lane; the env-variable gating also applies to `GPU=intel`)
- **Owner agent:** fortran-gpu-dev (Sonnet)
- **Review:** numerics-reviewer (code, acceptance file, scripts) · **Verification:** acceptance-verifier
- **Depends on:** the verification standard (merged 2026-10-29)
- **Base:** `gpu-port` 1a363070 (206 state + harness; the 207 τ code is **not** in this base)
- **Branch:** `gpu-port-212-hardening-b` (the old `gpu-port-212-hardening` keeps the superseded first pre-registration)
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
1. **One reader** of every `ALF_GPU_*` variable at setup (a small module or routine in `gpu_state_mod.F90`). The other
   files keep only calls to it; `Prog/Wrapgr_mod.F90` changes only inside its existing `#ifdef ALF_GPU` blocks.
   - "Set" means: defined and not empty, whatever the value (`ALF_GPU_NS=0` is set: it switches Newton–Schulz off).
     A run is a **debug run** if `ALF_GPU_SELFTEST`, `ALF_GPU_UPDTEST` or `ALF_GPU_DEBUG` is set to a value other than
     `0`; otherwise it is a **production run**. (Only for this debug/production decision does the value `0` of these
     three switches count as off; they are still listed in the debug-env line when set.)
   - **Allow-list:** in a production run every set `ALF_GPU_*` variable other than `ALF_GPU_DISABLE` and `ALF_GPU_DEBUG`
     is **refused** — the known debug/fault variables (`ALF_GPU_FAULT, _FAULT_STAB, _NS, _STABMIX, _STABCHECK, _TL,
     _TR, _VB, _NOUPDATE, _TRACE, _TIMING, _STABTIME, _EXIT`) and any unknown name (a typo, or 207's `_TAUSWAP`). The
     program prints exactly ` GPU debug env refused: <NAME> (production run; set ALF_GPU_DEBUG=1 for a debug run)` on
     stdout and in `info` (first refused name in sorted order) and stops with exit code 1 before the first sweep.
   - An unknown or unparsable **value** of a fault variable (e.g. `ALF_GPU_FAULT_STAB=xyz`) is an error in every run:
     ` GPU fault value invalid: <NAME>=<value>`, exit code 1 (no longer a warning).
   - Every GPU run prints at setup, on stdout and in `info`, exactly one line ` GPU debug env: none` (no `ALF_GPU_*`
     variable set) or ` GPU debug env: <NAME>=<value>, <NAME>=<value>` (sorted by name, every `ALF_GPU_*` variable that
     is set).
   - Every armed fault prints ` GPU fault armed: <NAME>=<value>` with the value exactly as given (stdout and `info`),
     for `ALF_GPU_FAULT` (204/205) and `ALF_GPU_FAULT_STAB` (206) alike.
2. **Build stamp line** at setup: ` GPU build: git <GIT_COMMIT_HASH, or none if git.h defines none> backend
   <CUDA|ONEMKL|HOST>` (harness builds come from `git archive` without `.git`, so `none` is expected there).
3. **Graded fault** for the stabilization self-test: `ALF_GPU_FAULT_STAB=eps:<x>` (debug runs only, as every fault)
   multiplies, in every compared CGR of the `ALF_GPU_SELFTEST=3` comparison, the element of largest |G| of the device G
   by `(1 + x)`. It must make the self-test fail at `x = 1e-9` (line ` GPU self-test FAILED: stab: ...`, exit code 1)
   and pass at `x = 1e-11` (0.1× the 1e-10 tolerance).
4. **Deterministic reductions:**
   - `k_colnorms`, `k_lognorm`, `k_prec`: a fixed combination order for floating-point **sums** (e.g. one thread per
     column with a serial sum, or a fixed-shape two-stage tree into scratch); max reductions are order-independent and
     may stay;
   - `la_init` (CUDA backend only): `cusolverDnSetDeterministicMode(…, CUSOLVER_DETERMINISTIC_RESULTS)` and
     `cublasSetAtomicsMode(…, CUBLAS_ATOMICS_NOT_ALLOWED)`; the HOST path of `alf_la_mod.F90` stays byte-identical in
     its preprocessed form, **including the integers `__LINE__` expands to** (`LA_FAIL` passes `__LINE__`; the harness
     masks only `"X.F90",N` line markers). Lines added or removed inside `ALF_LA_CUDA`-only blocks above HOST code must
     be compensated there (same line count per block, or a `#line` directive inside the CUDA-only region). The same
     holds for `Prog/Wrapgr_mod.F90` outside its `#ifdef ALF_GPU` blocks.
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
- **Kernel self-tests** 204 (`SELFTEST=2`) and 205 (`UPDTEST=1`) on the default build (their code does not depend on the
  stabilization variant), 206 (`SELFTEST=3`) default and LOG, PASS with their numbers ≤ 1e-10; the 204/205/206 fault
  controls and the graded fault behave as specified, and every armed fault shows its armed line.
- **Refusal:** every listed variable alone in a production run is refused with the exact line and exit code 1; with
  `ALF_GPU_DEBUG=1` the same variable is accepted (sampled: `NS`, `TL`). An unparsable fault value (`ALF_GPU_FAULT=xyz`,
  `ALF_GPU_FAULT_STAB=xyz`) gives the exact invalid-value line and exit code 1; a valid value does not.
- The 204/205 fault controls are sampled with one fault class each (`hop`, `upd`). The report quotes the CGR numbers of
  the `eps:1e-11` runs (expected ≈ 1e-11, showing that the fault was applied).
- **Production line:** every production GPU run prints ` GPU debug env: none` and the build stamp line (harness builds:
  ` GPU build: git none backend CUDA`).
- **Determinism:** the same binary twice gives bit-identical `data.h5` and full-precision Precision lines: 6×6 (200
  one-sweep bins) default and LOG, 12×12 β = 10 (8 bins × 20 sweeps) default and LOG.
- **Same seeds vs the CPU oracle:** 6×6 default and LOG, 12×12 default: no bin above 1e-4, overall and τ line ≤ 1e-9.
  (12×12 LOG is ticket 213, decision 3b; 213 re-measures it on top of 212.)
- **Precision mean** (Green, τ) ≤ 1.5× the larger CPU mean, 12×12 β = 4, default and LOG.
- **compare2** GPU vs GPU=none on fresh seeds (`s8_U4_b2_cb` shortened to 120 bins × 15 sweeps: the base GPU build needs
  0.31 s/sweep there, so the full 250 × 40 would take ~51 min; 120 × 15 ≈ 560 s leaves headroom under the 900 s
  timeout for the cost of the deterministic sums), with a CPU-vs-CPU null control.
- **`t_determinism`** passes under CUDA; its planted-ulp control fails.
- **GPU=none:** preprocessed changed CPU-path Fortran identical (`alf_la_mod.F90` and `Wrapgr_mod.F90`, 9 macro sets
  each = 18 comparisons; one build file, the `alf_la.tests` CMakeLists), ctest 53/53 (52 + `t_determinism`).
- Report: cost of the deterministic reductions — s/sweep at N = 144 and 400, default and LOG, from the declared timing
  runs `tm_*` (candidate) and `tmb_*` (base build) of `212.json` (no gate reads them: each number is quoted with its run's
  provenance stamp; 2 bins of 5–10 sweeps include setup time).

## Approvals
- First pre-registration `21e23a6a` (base `d53f9779`) was superseded before any implementation existed, after review
  (b) of the acceptance file (findings M1–M6, S1–S11). Second pre-registration `71bbcc93` was superseded before
  any implementation result existed, after its re-review (wrong base sha, the "set" rule vs `ALF_GPU_NS=0`, `__LINE__`
  in the HOST preprocessed output, an untested invalid-value rule, 204/205 numbers, compare2 timing margin). Current
  pre-registration commit: `PREREG_SHA` (2026-10-29) on
  `gpu-port-212-hardening-b`; base `1a363070`. Deviation to put to the owner: the compare2 positive control uses a
  shift of 8 standard errors, not the 5 of VERIFICATION.md §5 (at m = 11, α = 0.01 Šidák a shift of 5 is caught only
  ~60 % of the time); it needs the owner's quoted approval before T2. Any change of
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
