# 206 – Stabilization on the device: WRAPUR/WRAPUL, decompose, CGR, precision checks

- **Phase:** 2
- **Owner agent:** fortran-gpu-dev (review: numerics-reviewer)
- **Depends on:** 202, 203, 204 (104 and 108 decide the QR variants)
- **Status:** open (2026-10-28; runs in parallel with the other one, separate worktrees)

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
