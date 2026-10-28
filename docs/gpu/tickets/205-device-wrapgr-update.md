# 205 – WRAPGR up/down on the device: ratio, host decision, rank-k update

- **Phase:** 2
- **Owner agent:** fortran-gpu-dev (review: numerics-reviewer)
- **Depends on:** 202, 203, 204
- **Status:** open (2026-10-28; runs in parallel with the other one, separate worktrees)

## Goal
With `GPU=nvidia` and the GPU path active, the sequential local updates of a time slice (`WRAPGRUP`, `WRAPGRDO`) run
on a device-resident G: per-site wraps and the rank-`N_non_zero` update on the device, the ratio's small block and the
accept/reject decision on the host, with the RNG draws in exactly ALF's order (ADR-001 §5, MVP). Same compiler, same
seeds: the GPU run follows the CPU trajectory flip by flip until round-off first changes a decision (report when).

## Context
- CPU code: `Prog/Wrapgr_mod.F90` `WRAPGRUP` (≈ 95–150) and `WRAPGRDO` (≈ 190–245): per slice B_T wrap; per site
  `Op_Wrapup`/`Op_Wrapdo` (N_type 1), proposal (`nsigma%flip`, `ham%S0`, `T0_proposal`, `ranf_wrap`), `Upgrade2`,
  then N_type 2 (nothing for diagonal operators). `Prog/upgrade_mod.F90` `Upgrade2` (105–302): ratio from
  `GR(P(n),P(m),nf)` (op_dim × op_dim block), `alpha`, `Gama`, `N_SUN`, phase; on accept the update via `u, v, x_v,
  y_v` (row `GR(P(n),:)`, column `GR(:,P(1))`), ZGERU for op_dim 1, ZGEMM otherwise; `nsigma%f` updated.
- Design (keep the CPU arithmetic for everything decided on the host):
  - guarded early branches (`#ifdef ALF_GPU`, `if (gpu_active)`) in `WRAPGRUP`/`WRAPGRDO` and in `Upgrade2`, calling
    a new module `Prog/gpu_update_mod.F90`;
  - in `Upgrade2`'s GPU branch: copy the op_dim × op_dim block of G per flavour D→H (documented per-site transfer),
    run **the same ratio/weight code** on it (refactor the shared host part into an internal routine only if that
    leaves the CPU path's object code unchanged – otherwise duplicate it, with a comment pointing to the original),
    draw `ranf_wrap` at the same point, and on accept apply the update on the device (`la_geru` / `la_gemm` or own
    kernels for u, v, x_v, y_v); update `nsigma%f` on the host and the device copy (one scalar H→D);
  - per-site `Op_Wrapup`/`Op_Wrapdo` and the slice B_T wraps: the 204 device routines; the field index for the
    per-site wraps comes from the host (no stale `nsigma_d` reads);
  - 203 header notes for 205: `Op_V%alpha` in the ratio (host side, fine), `nsigma_d` staleness, per-flavour N/P.
- **Interim residency (until 206–208 land):** the rest of the sweep (stabilization, τ pass, measurements) still runs
  on the host. On entry to the GPU `WRAPGRUP`/`WRAPGRDO` copy G H→D, on exit D→H. These are **interim** transfers:
  mark them in the code (`! INTERIM 205: removed by 206/208`) and add them to the PLAN.md transfer table as interim.
  The host `GR` must not be read between entry and exit.
- Rules: ADR-002 §1 (synchronous target regions, single host thread), NaN-safe comparisons, no device allocation in
  the loop.

## Files allowed
- `Prog/gpu_update_mod.F90` (new), `Prog/Wrapgr_mod.F90`, `Prog/upgrade_mod.F90` (guarded branches only),
  `Prog/gpu_state_mod.F90`, `Prog/gpu_wrap_mod.F90` (small additions), `Prog/Makefile` (conditional source),
  `docs/gpu/PLAN.md` (transfer table rows only), `docs/gpu/scripts/nividia/t205_*`
Anything else: stop and report.

## Acceptance
1. **One-slice kernel test** (self-test mode): from the same G and the same saved RNG state (`random_seed(get/put)`),
   run the CPU `WRAPGRUP` then `WRAPGRDO` for one slice and the GPU versions: identical accept/reject pattern;
   G relative difference ≤ 1e-10 after the slice (report the value); `Phase` equal to ≤ 1e-12. Configs `hub6x6_cb`,
   `hub6x6_dense`, `hub12x12_cb`, `hub6x6_cb_mzoff`; default and LOG.
2. **Trajectory test** (validation level 2): `GPU=nvidia` (GPU path active) vs `GPU=none`, both nvfortran, same seeds,
   `hub6x6_cb` and `hub12x12_cb`, 20 sweeps: report the first sweep/slice/site where a decision differs (or "none"),
   and the per-bin differences of `Ener_scal` before that point (expected ≲ 1e-10 relative). `info` precision lines
   reported side by side.
3. **Statistics:** one short-autocorrelation config from `testsuite/gpu_reference/configs_stat.py` (8×8, β = 2; pick
   the size that runs ≤ 15 min single-rank on the GPU path), GPU vs a `GPU=none` nvfortran run with **different
   seeds**, `compare2.py`: must PASS (not INCONCLUSIVE). If the run is too long, reduce bins and say so; never
   loosen the test.
4. **`GPU=none`:** GNU MPI `ALF.out` code sections unchanged vs HEAD (same method as 203), ctest 52/52.
5. **`GPU=intel`** compiles on `intel` (build only).
- Report: diff summary, the per-site transfer/launch count, s/sweep GPU vs CPU (N = 144, 400 with interim copies –
  expected slow, report honestly), all tables above.

## nividia / intel rules
As in 204 (work dir `~/alf-gpu/work/t205/`, ≤ 4 cores, GPU only while running, ≤ 15 min per run, ≤ 4 GB disk,
`bash` scripts, guarded `rm -rf`, intel: load check, `level_zero:0`, no plain `sycl-ls`, no sudo). Ticket 206 runs
concurrently in its own work dir and edits other files (`wrapur/wrapul`, `cgr1`, `udv_state`); don't touch them.

## Out of scope
Stabilization (206), τ pass (207), measurements (208), removing the interim copies, device-side decision (401–404),
global moves, MPI.
