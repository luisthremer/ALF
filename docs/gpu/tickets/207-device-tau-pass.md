# 207 – τ pass on the device: PROPR/PROPRM1, WRAPUR, CGR2_2, precision check

- **Phase:** 2
- **Owner agent:** fortran-gpu-dev (review: numerics-reviewer)
- **Depends on:** 204, 205, 206
- **Status:** open (2026-10-28)

## Goal
With `GPU=nvidia`, `Ltau = 1` and the GPU path active, the time-displaced pass `TAU_M` runs on the device:
- the four Green functions `GT0`, `G0T`, `GTT`, `G00` are device-resident during the pass;
- the per-slice propagation (`PROPR`, `PROPRM1`) runs on the device;
- at each stabilization point, the right propagation (`WRAPUR` + decompose, from the 206 blocks) and `CGR2_2` run on the
  device, followed by `Control_Precision_tau` as a device reduction.

Observables (`ham%OBSERT`) stay on the host until 208: per time slice the four matrices are copied D→H. That copy is
**interim**.

Like the stabilization (206), the default path uses the pre-pivoted QR (ADR-001 §4, approved), and STABLOG uses host
ZGEQP3 (`la_geqp3_host`, 108).

## Context
- **CPU code:**
  - `Prog/tau_m_mod.F90`: `TAU_M` (56–211), `PROPR` (215–236), `PROPRM1` (240–263).
    - At NT = 0: G00 = GT0 = GTT = GR, and G0T = −(1 − GR).
    - Per slice: `PROPR(GT0)`, `PROPRM1(G0T)`, `PROPRM1`+`PROPR(GTT)`, `Hop_mod_Symm` copies (if `Symm`), then `OBSERT`.
    - At `Stab_nt(NST)`: `WRAPUR(NTST, NT1, udvr)` on a fresh right state, `CGR2_2(..., udvr, udvst(NST))`, and four
      `Control_Precision_tau` calls against the propagated matrices.
  - `Prog/cgr2_2_mod.F90` `CGR2_2` (196–806):
    - It handles a 2N × 2N system with two branches, `D1(1) > D2(1)` and the other one.
    - It builds `HLPB1`/`HLPB2`, then calls `QDRP_decompose` (ZGEQP3), `solve_extended_System` (ZUNMQR, ZTRSM, D
      scaling, permutation) and `get_blocks`.
    - STAB3/STABLOG scale handling: lines 349–393.
    - The test-only `-DQR_PREPIVOT` (variant 1, the final `#else` of the `QR_PREPIVOT` blocks) is the CPU reference
      of the pre-pivoted algorithm.
  - `Prog/control_mod.F90` `Control_Precision_tau` (372–383): `COMPARE` max/mean, plus `NCG_tau`, `XMAX_tau`,
    `XMEAN_tau`.
  - Call sites: `Prog/main.F90` ≈ 725, 752, 994 (all three in a guard for the backward sweep or the final step).
- **Device building blocks:**
  - 204: `gpu_mmthr`, `gpu_mmthl_m1`, `gpu_Vop_R/L`, the slice composites.
  - 206: `gpu_stab_mod` WRAPUR blocks, decompose, the device `udvst` slots, and the precision reduction pattern
    (`Control_PrecisionG_vals`).
  - 203: `GT0_d`/`G0T_d`/`GTT_d` are already allocated when `Ltau = 1`. Add a `G00_d` or use scratch.
  - `alf_la`: `la_qrpp`, `la_geqrf`, `la_ormqr`, `la_trsm`, `la_getrf`/`getrs`, `la_geqp3_host`.
  - 206 found that cuSOLVER's Q needs one Newton–Schulz step to keep precision at CPU level. Use the same scheme in
    CGR2_2 if Q is formed explicitly; if Q is only applied (`ormqr`), measure first.
- **What 207 removes:** the 206 interim D→H export of the stored left states for `TAU_M` (`gpu_stab_mod` ≈ 322–340).
  `TAU_M` must read the device `udvst` slots directly. Keep the export for the host τ path (τ path inactive).
- **Interim residency:**
  - Per slice, D→H of `GT0`, `G0T`, `GTT`, `G00` (4 · 16N² B per flavour) for the host `OBSERT`.
  - `Symm` and the observables run on those host copies.
  - At the start, G (GR) comes from wherever 205/206 left it; document which copy is authoritative.
  - Mark every interim transfer `! INTERIM 207: removed by 208` and add it to the PLAN.md transfer table.
- **gpu_supported:** `Projector` (`tau_p`), Langevin/HMC and `reconstruction_needed` are already refused; keep it that
  way.
- **Rules:**
  - ADR-002 §1: synchronous target regions, single host thread, `alf_la` ordered.
  - NaN-safe comparisons.
  - No device allocation in loops.
  - Every transfer documented in the code.
  - GPU runs one at a time on the A30: two processes are 13× slower each (DISCREPANCIES.md).
  - For `GPU=none` identity after `main.F90` edits, compare the preprocessed sources (`gfortran -E -P`, `diff -B`).
    The disassembly differs, because blank lines shift the line numbers embedded in runtime-error strings (206).

## Files allowed
- `Prog/gpu_tau_mod.F90` (new)
- `Prog/tau_m_mod.F90` (guarded early branch in `TAU_M` only)
- `Prog/main.F90` (guarded branches at the three `Tau_m` call sites only)
- `Prog/control_mod.F90` (guarded `Control_Precision_tau_vals` only)
- `Prog/gpu_stab_mod.F90` (small: slot access for the τ pass, skip the interim `udvst` export when the τ path is active)
- `Prog/gpu_state_mod.F90`
- `Prog/gpu_wrap_mod.F90` (small additions)
- `Prog/Makefile` (conditional source)
- `docs/gpu/PLAN.md` (transfer table rows only)
- `docs/gpu/scripts/nividia/t207_*`

`Prog/cgr2_2_mod.F90` and `Prog/wrapur_mod.F90` stay untouched (oracles). Anything else: stop and report.

## Acceptance
1. **Kernel tests** (self-test mode, e.g. `ALF_GPU_SELFTEST=4`, fixed fields = the initial configuration). Oracle for the
   default path: the host `TAU_M`/`CGR2_2` in a `-DQR_PREPIVOT` build, as in 206. Oracle for STABLOG: upstream.
   - one full τ pass device vs host from the same G and the same `udvst`: `GT0`, `G0T`, `GTT`, `G00` at **every** time
     slice ≤ 1e-10 relative (report the max per slice block);
   - `CGR2_2` alone from the same UDV states ≤ 1e-10. Both branches (`D1(1) > D2(1)` and else) must be hit; report
     which configs hit which branch, and add a forced-branch test if a branch is not reached naturally;
   - the device τ precision reduction equals `Control_Precision_tau` on the same matrices to ≤ 1e-12 relative;
   - configs: `hub6x6_cb`, `hub6x6_dense`, `hub12x12_cb`, `hub12x12_cb` β = 10, `hub6x6_cb_flux` (complex),
     `hub6x6_cb_symm` (Symm path for the observables); default and LOG.
   - negative controls (faults) must fail.
   - If a config outside this set misses 1e-10, report it with the measured 1-ulp sensitivity, as 206 did. Never change
     the tolerance.
2. **Same seeds, full device path** (updates + stabilization + τ on the device) vs the CPU oracle (`QR_PREPIVOT` build;
   CPU LOG for LOG): `hub6x6_cb` 200 sweeps and `hub12x12_cb` 300 sweeps, one sweep per bin, `Ltau = 1`. Report the
   first bin above 1e-9 and above 1e-4 over **all** datasets, including the τ-resolved ones. Expected: no divergent
   decision.
3. **Precision:** `Precision tau` and `Precision Green` mean and max for `hub12x12_cb` β = 4 and β = 10, default and LOG,
   next to `GPU=none` and the `QR_PREPIVOT` CPU build (CPU LOG for LOG). The GPU tau max must not exceed the larger CPU
   max by more than 1.5×. Choose the sweep count so that each GPU run takes ≤ 15 min, and say what you chose.
4. **Statistics:** `compare2.py` on one config with `Ltau = 1` (e.g. `s8_U4_b2_cb` or smaller), GPU (full device path)
   vs `GPU=none` with different seeds. It must PASS, with τ-resolved observables included in the tested set. If the run
   is too long, reduce bins and say so. This is also the first combined compare2 for 205 + 206.
5. **`GPU=none`:** preprocessed sources of the changed files token-identical to HEAD (`gfortran -E -P`, `diff -B`);
   `data.h5` identical to HEAD on 3 configs (GNU noMPI HDF5, same seeds); ctest 52/52.
6. **`GPU=intel`:** builds on `intel`, and a 6×6 run completes (the τ path is inactive until 302; say so).
- **Report:** diff summary; transfers per τ pass (bytes; interim and permanent separately); time per τ pass GPU vs CPU
  at N = 144, 400, 1024, default and LOG (`ALF_GPU_TAUTIME=1` or similar); s/sweep with `Ltau = 1` before (206 state)
  and after; all tables above.

## nividia / intel rules
As in 204–206:
- work dir `~/alf-gpu/work/t207/` (TMPDIR inside it);
- ≤ 4 cores; ≤ 15 min per run; ≤ 4 GB disk;
- **one GPU process at a time** (check `nvidia-smi` before each GPU run);
- `bash` scripts, guarded `rm -rf` on literal paths;
- intel: load check first, `level_zero:0`, no plain `sycl-ls`, no sudo.

Clean up at the end: keep scripts and logs, remove trees, runs and binaries, and make sure the GPU is released.

## Out of scope
- observables on the device and removing the per-slice copies (208);
- `tau_p` / projector;
- Langevin/HMC;
- `Symm` on the device;
- STAB1/2/3;
- MPI;
- changes to the CPU path.
