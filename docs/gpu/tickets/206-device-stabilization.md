# 206 – Stabilization on the device: WRAPUR/WRAPUL, decompose, CGR, precision checks

- **Phase:** 2
- **Owner agent:** fortran-gpu-dev (review: numerics-reviewer)
- **Depends on:** 202, 203, 204 (104 and 108 decide the QR variants)
- **Status:** done (2026-10-28)

## Goal
With `GPU=nvidia` and the GPU path active, the stabilization steps of the sweep run on the device:
- the B-product propagation of the UDV states (`WRAPUR`, `WRAPUL`);
- their decomposition (`udv_state%decompose`), default and STABLOG;
- the Green function from two UDV states (`CGR`, both NVAR);
- the precision check against the wrapped G (`Control_PrecisionG`), as a device reduction with one scalar D→H.

The default stabilization uses the pre-pivoted QR (104, approved). STABLOG uses its own log-norm pre-sort plus an
unpivoted QR in decompose, and host ZGEQP3 in CGR via `la_geqp3_host` (108). `Precision Green` stays at CPU level.

## Context
- CPU code:
  - `Prog/wrapur_mod.F90`, `Prog/wrapul_mod.F90`: per Nwrap slices, B products on `udvr%U` / `udvl%U` via
    `Hop_mod_mmthR`/`mmthlc` and `Op_mmultR`, then `%decompose`, `%testscale`.
  - `Prog/udv_state_mod.F90` decompose (≈ 450–565): default = `QDRP_decompose` (ZGEQP3 with free pivoting);
    STABLOG = log-norm sort (`log(DZNRM2)+L`, `ZSWAP`, phase sign) then `QDRP_decompose` with `IPVT = 1`, i.e. an
    unpivoted QR.
  - `Prog/QDRP_decompose_mod.F90` (ZGEQP3, |R_ii|, `Pivot_phase`).
  - `Prog/cgr1_mod.F90` `CGR` (STAB3/STABLOG branch ≈ 180–640: QR with pivoting ≈ 440–500, ZUNMQR, ZTRSM, ZLAPMR/T,
    determinant/phase). The test-only macro `QR_PREPIVOT` (104/108) holds the CPU reference of the pre-pivoted
    variant 1.
  - Call sites: `Prog/main.F90` ≈ 637 (initial), ≈ 750–870 (sweep: `WRAPUR`/`WRAPUL`, `CGR`,
    `Control_PrecisionG`).
- Device building blocks: `alf_la` (`la_qrpp`, `la_geqrf/ungqr/ormqr`, `la_trsm`, `la_getrf/getrs`, `la_geqp3_host`,
  `la_gemm`), 204 slice propagation, `gpu_state` `udv_*` buffers and scratch.
- **Oracle for the default path:** the GPU QR differs from upstream ZGEQP3 by design. Compare it against a **CPU
  build with `-DQR_PREPIVOT`** (variant 1, the same algorithm), which must agree to ≤ 1e-10. Compare against
  upstream only at the precision level: `Precision Green` of the same order. STABLOG CGR uses host ZGEQP3, so it
  compares directly with upstream STABLOG.
- **Interim residency** (as 205): the parts not yet on the device copy at their boundaries. Mark each transfer
  `! INTERIM 206: removed by 207/208` and add it to the PLAN.md transfer table as interim. Ticket 205 runs
  concurrently and wraps `WRAPGRUP`/`WRAPGRDO` with the same scheme: G is host-authoritative between modules until
  208.
- Rules: ADR-002 §1 (synchronous target regions, single host thread, `alf_la` ordered), NaN-safe comparisons, no
  device allocation in loops; every transfer documented in the code.

## Files allowed
- `Prog/gpu_stab_mod.F90` (new), `Prog/wrapur_mod.F90`, `Prog/wrapul_mod.F90`, `Prog/udv_state_mod.F90`,
  `Prog/cgr1_mod.F90` (guarded early branches only), `Prog/main.F90` (guarded branches at the stabilization call
  sites only), `Prog/control_mod.F90` (guarded hook for the precision value only), `Prog/gpu_state_mod.F90`,
  `Prog/gpu_wrap_mod.F90` (small additions), `Prog/Makefile` (conditional source), `docs/gpu/PLAN.md` (transfer table
  rows only), `docs/gpu/scripts/nividia/t206_*`
Anything else: stop and report.

## Acceptance
1. **Kernel tests** (self-test mode, fixed fields = the initial configuration):
   - `WRAPUR`+decompose and `WRAPUL`+decompose over one Nwrap block: U, D (or L), V vs the CPU oracle ≤ 1e-10
     relative (default: oracle = `QR_PREPIVOT` CPU build or an in-binary host call of the same algorithm; STABLOG:
     upstream);
   - `CGR` (NVAR 1 and 2) from the same UDV states: G and phase vs the oracle ≤ 1e-10 relative;
   - the device precision reduction equals `Control_PrecisionG` on the same matrices to ≤ 1e-12 relative;
   configs `hub6x6_cb`, `hub6x6_dense`, `hub12x12_cb`, and one larger β (e.g. 12×12, β = 10) for scales;
   default and LOG.
2. **Sweep with device stabilization** (GPU path, interim copies), `hub12x12_cb` and a β = 10 case, default and
   LOG: `Precision Green`/`tau`/phase mean and max reported next to `GPU=none` and the `QR_PREPIVOT` CPU build;
   the GPU max must not exceed the larger of the two CPU maxima by more than 1.5× (104 criterion). Same seeds:
   report the first divergent sweep against the `QR_PREPIVOT` CPU build.
3. **Statistics:** one `configs_stat.py` config, GPU vs `GPU=none` with different seeds, `compare2.py` PASS.
4. **`GPU=none`:** GNU MPI `ALF.out` code sections unchanged vs HEAD (method of 203), ctest 52/52.
5. **`GPU=intel`** compiles on `intel` (build only).
- Report: diff summary, transfers per stabilization (bytes), s per stabilization step GPU vs CPU (N = 144, 400,
  1024; default and LOG), all tables above.

## nividia / intel rules
As in 204 (work dir `~/alf-gpu/work/t206/`, ≤ 4 cores, GPU only while running, ≤ 15 min per run, ≤ 4 GB disk,
`bash` scripts, guarded `rm -rf`; intel: load check, `level_zero:0`, no plain `sycl-ls`, no sudo). Ticket 205 runs
concurrently and edits `Wrapgr_mod`, `upgrade_mod`; don't touch them. Coordinate only through the files allowed.

## Out of scope
WRAPGR/updates (205), τ pass and CGR2_2 (207), measurements (208), removing interim copies, STAB1/STAB2 paths (the
GPU path refuses them in `gpu_supported` – add that check if missing), MPI.

## Resolution (manager, 2026-10-28)
- Delivered: `Prog/gpu_stab_mod.F90` (device WRAPUR/WRAPUL blocks via the 204 composites, decompose default and
  STABLOG, `testscale`, CGR NVAR 1/2, precision reduction; udvl/udvr/udvst device slots; self-test
  `ALF_GPU_SELFTEST=3`, timing `ALF_GPU_STABTIME=1`, faults `ALF_GPU_FAULT_STAB=perm|phase|nan`, debug aids
  `ALF_GPU_STABCHECK`/`STABMIX`), guarded branches at the stabilization sites of `Prog/main.F90` (pure additions; the
  `#else` text is HEAD's), `Control_PrecisionG_vals` in `Prog/control_mod.F90`, a `QR_PREPIVOT` (variant 1, test-only)
  block in `Prog/udv_state_mod.F90`; `gpu_supported` refuses STAB1/2/3. `wrapur_mod`, `wrapul_mod`, `cgr1_mod` are
  untouched (TAU_M keeps the host routines, so the device branches live in main). Interim transfers in PLAN.md.
- Numerics finding: cuSOLVER's Q is ≈ 5× less orthogonal than LAPACK's; the device decompose forms Q with
  `la_ormqr` on I and applies one Newton–Schulz step (`ALF_GPU_NS=0` switches it off; DISCREPANCIES.md). Cost at
  N = 400: WRAPUR+decompose 13.3 ms with, 12.8 ms without.
- numerics-reviewer: 2 must-fixes (NVAR/TEST lines reordered outside the guards in main; device slots left dirty after
  the self-test/timing hook) and should-fixes (NS switch, STABLOG faults inert, `sync_fields` redundant with 205,
  debug scaffolding) – all applied by the agent.
- Merge with 205: one Makefile line; `sync_fields` now only on the host-update path.
- Manager re-run (nividia, merged 205+206 tree, fresh builds):
  - self-test PASSED: `hub6x6_cb`, `hub6x6_dense`, `hub12x12_cb`, `hub12x12_cb` β=10, `hub12x12_dense` β=10,
    `hub6x6_cb_flux`, default (`-DQR_PREPIVOT` oracle) and LOG; faults `perm`, `phase`, `nan` exit 1 on both paths;
    a run continuing after the self-test gives `data.h5` identical to a run without it;
  - same seeds, updates **and** stabilization on the device vs the CPU oracle (`QR_PREPIVOT` build, resp. CPU LOG):
    hub6x6 200 sweeps max deviation 1.5e-10 / 2.5e-10 (default / LOG), hub12x12 300 sweeps 5.0e-10 / 1.2e-9 (LOG
    above 1e-9 from sweep 184), no divergent decision; host-update path with device stabilization 1.5e-10;
  - precision (`info`, one seed set), GPU max vs the larger CPU max (criterion ≤ 1.5×):

    | case | GPU Green mean / max | GPU=none | QR_PREPIVOT CPU | ratio |
    |---|---|---|---|---|
    | 12×12 β=4 default | 5.5e-11 / 2.27e-6 | 3.4e-11 / 1.23e-6 | 4.1e-11 / 1.57e-6 | 1.45 |
    | same, `ALF_GPU_NS=0` | 1.0e-10 / 3.80e-6 | | | 2.42 (fails) |
    | 12×12 β=10 default | 5.4e-10 / 4.48e-5 | 1.3e-10 / 3.26e-6 | 9.7e-10 / 1.40e-4 | 0.32 |
    | 12×12 β=4 LOG | 6.0e-11 / 2.11e-6 | CPU LOG 4.9e-11 / 1.42e-6 | | 1.48 |
    | 12×12 β=10 LOG | 9.3e-10 / 1.35e-4 | CPU LOG 8.8e-10 / 1.18e-4 | | 1.15 |

    Phase max 0 everywhere; tau precision at CPU level.
  - 205 update self-test and 204 wrap self-test still PASS on the merged build; `GPU=none` (GNU noMPI HDF5) `data.h5`
    identical to the 205 commit on 3 configs; ctest 52/52 (GNU devel). GNU MPI `ALF.out`: section sizes equal, but
    `.rodata`/`.data` offsets and ~10k `lea` displacements differ because the added blank lines shift the line numbers
    gfortran embeds in runtime-error strings. The preprocessed `GPU=none` sources of `main.F90`, `control_mod.F90`,
    `udv_state_mod.F90` are token-identical to the 205 commit apart from blank lines (method for main-loop tickets:
    compare `gfortran -E -P` output with `diff -B`, not the disassembly);
  - `GPU=intel` builds on `intel` (25 s); a 6×6 run completes with the 205 update path ACTIVE on the B580 and the
    stabilization INACTIVE (ticket 302), Precision Green 2.7e-11 / 8.0e-8.
- Agent results not re-run: compare2 PASS (`s8_U4_b2_cb`, 250 × 40 sweeps, device stabilization with host updates);
  timings (N = 1024 default: WRAPUR+decompose 52 ms vs 1.67 s CPU, CGR 58 ms vs 2.30 s; STABLOG CGR 412 ms, host
  ZGEQP3).
- **Open:**
  - the β = 4 maxima pass 1.5× only narrowly (1.45, 1.48; one seed set; the max is an extreme-value statistic, over 8
    seeds at 6×6 the mean of the maxima agrees to 2 %). Watch in 208 / Phase 3 validation with several seeds;
  - self-test outside the required set: `hub20x20_cb` β=10 default G NVAR 1 differs by 3.5–5.3e-9 (> 1e-10). The same
    test measures G's sensitivity to a 1-ulp change of U_R at 6e-10 (host) / 4.3e-9 (device): conditioning, not a bug;
    LOG agrees to 5e-12. Tolerance not changed; open question for the owner whether large-β kernel tests should be judged against the measured 1-ulp sensitivity;
  - two GPU processes on one A30 run 13× slower (time-slicing): GPU runs one at a time; MPS for production;
  - no combined compare2 (updates + stabilization on the device): one GPU run would exceed 15 min at the current
    speed; due with 208.
