# 204 – Device B_T (checkerboard, dense) and Op_V wraps

- **Phase:** 2
- **Owner agent:** fortran-gpu-dev (review: numerics-reviewer)
- **Depends on:** 202 (`alf_la`, done), 203 (`gpu_state` flat mirror)
- **Status:** done (2026-10-28; reviewed, re-run by the manager)

## Goal
Device routines that apply the hopping exponentials and the Op_V operators of one time slice to device-resident
matrices, exactly as the CPU routines do, built only on the 203 flat arrays and `alf_la`. A self-test with fixed HS
fields shows agreement with the CPU routines to ≤ 1e-10 relative (expected ≈ 1e-15) and gives per-slice timings.
No sweep uses them yet (205–207 wire them in); results stay identical.

## Context
- CPU routines to mirror (one device routine each, same semantics, operating on a device matrix of shape
  (Ndim, ncol)):
  - hopping, `Prog/Hop_mod.F90:143–250`: `Hop_mod_mmthr` (left, elements Ncheck..1, `lmult`), `mmthr_m1` (left
    inverse, 1..Ncheck), `mmthl` (right, 1..Ncheck), `mmthlc` (left with `lmult`, 1..Ncheck), `mmthl_m1`. Element
    application `RealExpOpT/CmplxExpOpT` `lmult/rmult/...` (`Prog/OpTTypes_mod.F90:150–300`, `ZDSLSYMM` on the
    element's `P`; skipped when `g*g <= Zero`). Order of elements must match exactly (HARD H5).
  - Op_V, `Prog/Operator_mod.F90`: `Op_Wrapup` / `Op_Wrapdo` (743–951; diag, types 1/2, `N_type` 1 and 2),
    `Op_mmultR` / `Op_mmultL` (555–717; `cop` 'n' and 'c', `sign`).
  - Call sites that 205–207 will replace: `Wrapgr_mod.F90:110–145, 195–240` (G wraps around the updates),
    `wrapur_mod.F90:80–117`, `wrapul_mod.F90:85–124`, `tau_m_mod.F90:229–259`.
- **Kernels:**
  - checkerboard: elements with disjoint `P` that are consecutive in the application order touch disjoint rows
    (left) / columns (right), so they can run in one kernel launch and give the same result as sequential
    application. Compute these groups once at setup (extend `gpu_state` with group offsets; respect the order).
    Report launches per B_T and the group count for 6×6, 12×12, 20×20.
  - dense (one element with full `P`): real exp(T) → right multiplication via `la_gemm_cr` (DGEMM on the
    interleaved view, 101); left multiplication: choose between a complex copy of exp(T) with `la_gemm` and a
    transpose-based variant; justify with a timing at N = 400 and 1024. Complex exp(T): `la_gemm`.
  - Op_V diag (types 1/2, no `g_t`): row/column scaling with the 203 tables, indexed by the field value on the
    device (`nsigma_d`). For a full slice (all `n` in order, as in wrapur/wrapul/tau_m) provide a fused kernel that
    applies the factors of all operators to each row/column **in ALF's operator order**, so that overlapping `P`
    sets give bit-identical products; per-operator routines remain for the per-site wraps of 205.
- Lessons from the 203 review (read the `gpu_state_mod.F90` header):
  - the flat hopping matrices encode exactly what ALF reads: lower triangle filled from the upper for complex
    elements and for real nd ≤ 8; real nd > 8 keeps the full stored matrix, because `ZDSLSYMM` then calls `ZLARCM`
    on the full matrix (DISCREPANCIES 2026-10-28). nd ≤ 8 (checkerboard bonds) uses hand-written loops: aim for
    bit-identity. nd > 8 goes through `ZLARCM` (real) / `ZHEMM` 'U' (complex): device results ≈ 1e-15, not
    bit-identical; the bound is 1e-10;
  - Op_V: use `nint(real(nsigma_d(n,t)))` for the field index; honour the per-operator skip flag (tiny |g|) in the
    `Op_mmultL/R` kernels but **not** in the `Op_Wrapup/Wrapdo` kernels, as ALF does;
  - complex nd > 8 goes through `ZHEMM`, which ignores the imaginary part of the diagonal (the nd ≤ 8 loops use
    it; ≈ 1e-17): zero it for the device path or note the difference;
  - `Op_Wrapup`/`Op_Wrapdo` with `N_type = 2` apply nothing for diagonal operators;
  - comparisons must be NaN-safe (explicit loops, `ieee_is_finite`; never `maxval`/`max` for error measures).
- Rules (ADR-002 §1 amended): synchronous target regions only (no `nowait`/`depend`), single host thread; `alf_la`
  calls are ordered with them. No device allocation in these routines (scratch from `gpu_state`).
- New code in a new module `Prog/gpu_wrap_mod.F90`, compiled only with `ALF_GPU` (same Makefile mechanism as 203).

## Files allowed
- `Prog/gpu_wrap_mod.F90` (new), `Prog/gpu_state_mod.F90` (group arrays, scratch sizes, self-test hook),
  `Prog/Makefile` (conditional source), `Prog/main.F90` (only if the self-test hook needs it; guarded)
- `docs/gpu/scripts/nividia/t204_*.sh`
Anything else: stop and report.

## Acceptance
1. **Self-test** (`ALF_GPU_SELFTEST=2` or similar, `GPU=nvidia`, A30), after setup, fields = the initial `nsigma`
   (fixed, deterministic), random complex test matrices (seeded):
   - each hopping routine (5) on an Ndim×Ndim and an Ndim×N_part matrix: device vs the CPU routine, relative
     Frobenius difference ≤ 1e-10; report the measured maximum (checkerboard expected 0 or ≈ 1e-16);
   - `Op_Wrapup`/`Op_Wrapdo` (N_type 1, 2) and `Op_mmultR/L` ('n','c', both signs) for every operator and both
     field values: same bound;
   - a slice propagation as in `wrapur` (`Hop_mod_mmthr` + all `Op_mmultR`) over `Nwrap` slices, and the
     `wrapul` and `WRAPGRUP`/`WRAPGRDO` wrap sequences without updates (B G B⁻¹ for consecutive slices): device
     vs CPU ≤ 1e-10 relative after the full `Nwrap` (report growth per slice);
   - failure → nonzero exit with a message.
   Configs: `hub6x6_cb`, the dense 6×6 from 203, `hub12x12_cb`, and a 20×20 cb (N = 400), each default and LOG;
   Mz true and false if 203 covered both.
2. **Timing** (A30 vs the same routine on one CPU core in the same binary): µs per slice for the hopping
   (left/right), the fused Op_V slice, and the wrap sequence, at N = 144, 400 and 1024 (cb and dense). Compare
   with 101's kernel numbers (BENCHMARKS.md) and say whether the launches per slice match the group count.
3. **Unchanged physics:** `GPU=nvidia` without self-test vs `GPU=none` PGI on `hub6x6_cb`, same seeds, 20 sweeps:
   all `data.h5` datasets bit-identical.
4. **`GPU=none`:** GNU MPI `ALF.out` without the new module (nm), ctest 52/52 (GNU devel).
5. **`GPU=intel` compiles** on `intel` (build only; run the 6×6 self-test on the B580 if it works – report, not a
   blocker).
- Report: diff summary, kernel design (groups, launches, dense choice with numbers), self-test table, timings.

## nividia / intel rules
As in 203: develop on the laptop, rsync; `~/alf-gpu/work/t204/` on each machine, `TMPDIR` inside; ≤ 4 cores on
nividia, ≤ 8 on intel after checking load; GPU only while running; runs ≤ 15 min; scripts with `bash`; on intel
`ONEAPI_DEVICE_SELECTOR=level_zero:0`, never plain `sycl-ls`; keep the nividia work dir ≤ 4 GB; clean up.

## Out of scope
Updates and the per-site decision (205), stabilization (206), τ pass (207), measurements (208); `Symm` hopping for
observables; time-dependent hopping or couplings; any change of the CPU path or of results.

## Resolution (manager, 2026-10-28)
- Delivered: `Prog/gpu_wrap_mod.F90` (under `#ifdef ALF_GPU`): five hopping routines, per-operator Op_V
  (`Op_mmultR/L`, `Op_Wrapup/do`), fused Op_V slice kernels, composites for `wrapur`/`wrapul`/`WRAPGRUP`/`WRAPGRDO`,
  self-test (`ALF_GPU_SELFTEST=2`), timing (`ALF_GPU_TIMING=1`), fault injection (`ALF_GPU_FAULT`). `gpu_state`:
  checkerboard groups (maximal runs of disjoint elements; 4 groups per flavour for square cb, 7 with Symm), dense
  auxiliaries, per-site operator lists. Dense left: complex copy + ZGEMM for N ≤ 400, transposes + DGEMM above.
- Agent results (A30): self-test PASSED for 8 configs (cb, dense, complex flux cb and dense, Mz off, Symm, 12×12,
  20×20), default and LOG; max rel. diff ≤ 1.7e-16 per routine, wrap sequences over Nwrap = 10 ≤ 5.4e-16; the left
  hopping routines are bit-identical, right ones differ by FMA (1.1e-16). Timing per slice at N = 1024 cb: hopping
  280/174 µs (CPU 32.4/2.0 ms), fused Op_V 52 µs (CPU 11.8 ms), WRAPGR wraps 538 µs (CPU 47 ms); cb kernels
  1.2–1.5× 101's numbers, dense better than 101. B580: 6×6 cb self-test PASSED (`alf_la` parts skipped by the ONEMKL
  guard).
- numerics-reviewer: approve, no must-fix.
- Manager re-run (nividia): self-test PASSED for `hub6x6_cb`, `hub6x6_dense`, `hub12x12_cb`, `hub6x6_cb_mzoff`,
  `hub6x6_cb_symm` (default) and `hub6x6_dense`, `hub12x12_cb` (LOG); faults `hop`, `vexp`, `nan`, `order` exit 1;
  `GPU=nvidia` vs `GPU=none`, same seeds, 20 sweeps: `data.h5` identical, precision lines equal (3 configs); ctest
  52/52 (GNU devel).
- **Follow-ups (for 205 or a cleanup ticket):**
  - remove the experimental single-launch checkerboard kernel (`gpu_cb_fused`): slower on the A30 and 0.2 wrong on
    the B580 (cause not established; reviewer: legal OpenMP, suspect ifx). It is unreachable in production (default
    off, ifx disabled);
  - `ensure_la()` initialises `alf_la` lazily inside the first dense call: move `la_init` to `gpu_state_setup`;
  - fault magnitudes: `hop` 1e-9 is only just above the tolerance at N = 36; use 1e-6 as for `vexp`. The `nan` fault
    tests the up-front finite check, not the comparison;
  - the up/down round trip from a random G differs by ≈ 1.5e-10 device vs CPU, the same as a 1e-16 perturbation on
    the CPU (conditioning of B⁻¹(B G B⁻¹)B): reported, not a criterion;
  - dense elements with non-identity `P` (gather/scatter path) are untested (no config has one).
