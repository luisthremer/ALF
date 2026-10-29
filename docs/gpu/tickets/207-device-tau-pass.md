# 207 – τ pass on the device: PROPR/PROPRM1, WRAPUR, CGR2_2, precision check

- **Phase:** 2
- **Owner agent:** fortran-gpu-dev (review: numerics-reviewer)
- **Depends on:** 204, 205, 206
- **Status:** on hold: awaiting owner decision and verification standard (2026-10-29)

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
   next to `GPU=none` and the `QR_PREPIVOT` CPU build (CPU LOG for LOG). Choose the sweep count so that each GPU run takes
   ≤ 15 min, and say what you chose.
   - **Changed by the owner on 2026-10-29** (it applies to `Precision tau` and to the `Precision Green` max of 206).
     - Original rule: "the GPU tau max must not exceed the larger CPU max by more than 1.5×".
     - Why it was replaced: the max is the single largest event of a heavy-tailed quantity driven by rounding. A CPU run
       whose input is shifted by 1 ulp fails it too (see Resolution).
     - The max is now reported for information only.
   - **The mean** must stay within 1.5× of the larger CPU mean.
   - **Paired same-trajectory test:** the device τ pass and the host τ pass run on identical states (`ALF_GPU_TAUDIAG=2`,
     `t207_diagstat.py`). It must meet all three of these:
     - |mean log10(device/host)| ≤ 0.05 over the events;
     - the spread (sd) of log10(device/host) ≤ 1.25 × the spread between two 1-ulp-shifted host runs;
     - the event counts above 1e-7 and 3e-7 agree within 2·√(device + host).
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

Clean up at the end: make sure the GPU is released. **Replaced on 2026-10-29:** raw outputs, run directories and
binaries are kept until the owner's gate (the earlier rule "remove trees, runs and binaries" destroyed evidence of
acceptance 2–4).

## Out of scope
- observables on the device and removing the per-slice copies (208);
- `tau_p` / projector;
- Langevin/HMC;
- `Symm` on the device;
- STAB1/2/3;
- MPI;
- changes to the CPU path.

## On hold (2026-10-29): state, evidence, open items

207 is **not closed**. Closing needs the new verification standard (`gpu-port/process-lanes`) and the owner's
decision on acceptance 3. Until then no raw output, run directory or binary of 207 is deleted.

### Owner decisions on file
1. **2026-10-29, manager session, approval of the acceptance-3 rewrite** (verbatim): "Yes: for Precision tau (and the
   206 Green max), replace the 1.5x-max rule with the paired same-trajectory test as you proposed; report the max as
   information. Keep the diagnostic in the code. Then review, re-run and commit 207."
2. **2026-10-29, owner instruction after the paired test failed on the default path** (verbatim): "do not drop the sd
   condition, and do not commit 207." The default-path spread stays an open FAIL until the gate. Options for the gate:
   (a) find the cause first (`k_colnorms` pivot order in `la_qrpp`); (b) a new pre-registered criterion, calibrated
   with null and known-bias runs, only after owner approval; (c) keep 207 open.
3. **Caveats recorded with decision 1:** the paired thresholds (0.05, 1.25×, 2√n) were written after the β = 4 LOG
   paired data had been seen; the margin is 6–14 %; the reference spread moves with the chain start
   (`ALF_GPU_TAUDIAG_FROM` 330 vs 331: 0.168 vs 0.181); no paired run at β = 10 exists. Which acceptance-3 rule is in
   force is decided at the gate.
4. **Owner accepted (2026-10-28), pending review:** the `main.F90` hooks outside the three `Tau_m` call sites (setup,
   self-test hook, finalize; pure additions). The device branch is one early return inside `TAU_M`.

### Source states
| State | Commit (local WIP, pushed as `gpu-port/207-wip`) | `gpu_tau_mod.F90` sha256 prefix |
|---|---|---|
| S0 base (206) | 64ecde72 | – |
| S2 Sonnet final | 4d7a66e1 | 07333eba |
| S3 Opus diagnostic | bf25ccc3 | fe86b9f5 |
| S4 review fixes (candidate) | b44c3d6e | d5019652 |

### Evidence status per acceptance item (against S4)
"Stamped" = run from `git archive` of the commit, source sha and binary sha256 on the log line (`logs/mgr207c.txt`);
none of the runs so far used `env -i` or recorded the toolchain and environment, so none meets the new standard.
| # | Item | Status | Evidence |
|---|---|---|---|
| 1 | Kernel self-tests + faults, 6 configs × {default, LOG} | valid for S4 (unstamped: tree hash and binary symbol `upd_to_i` match S4) | `mgr207b.txt`, `selftest_*.txt`, `negctl_*.txt` (11:03); max 1.85e-11, faults 4/4 rc = 1 |
| 2a | Same-seed 6×6, default + LOG | stamped S3 vs S4 bit-identical (`hub6x6_cb_fd`); S3 numbers in `mgr207.txt`, raw lost | `mgr207c.txt` (1) |
| 2b | Same-seed 12×12 default | **stamped S4**: max 3.103e-10, τ 1.742e-10, no bin > 1e-9 | `runs/c_nv_def`, `runs/c_nv_none_pp` |
| 2c | Same-seed 12×12 LOG | stale (S3 binary), raw lost: max 1.214e-9 from bin 184, τ 1.268e-10 | `mgr207.txt` |
| 2d | Sonnet extras (symm/flux/dense fd, NOUPDATE, Ltau = 0) | lost (report only) | – |
| 3a | Precision means (1.5× rule) | stale (S2/S3), raw lost; recomputed from the tables: all pass, worst Green mean 1.36× (β = 4 default) | `gpuprec.txt`, `finprec.txt` |
| 3b | Paired, LOG s1 | **stamped S4**: PASS (ratio 1.062) | `runs/cpair_nv_log` |
| 3c | Paired, default s0 | **stamped S4: FAIL**, recomputed from the raw TAUDIAG lines: n = 560, mean −0.0105, sd(log10 dev/host) 0.1653, reference sd(host+1ulp(b)/host+1ulp(a)) 0.0969, **ratio 1.707 > 1.25** (also 1.79 against host+1ulp/host). Conditions a and c pass | `runs/cpair_nv_def_pp` |
| 3d | Paired at β = 10 | missing | – |
| 3e | Max ratios (information) | stale (S2), raw lost | `gpuprec.txt` |
| 4 | compare2 `s8_U4_b2_cb` | stale (S2), raw lost; JSON self-consistent (PASS, Σt² 10.7, p 0.47) | `compare2_s8_U4_b2_cb.*` |
| 5 | GPU=none (ppcmp, ident, ctest) | valid by source identity (all GPU=none sources of S3 ≡ S4), unstamped | `mgr207.txt` |
| 6 | GPU=intel build + 6×6 run | **blocked**: oneAPI on `intel` is mixed (ifx 2025.3, MKL 2026.2), link fails on SYCL symbols | `make_intel.log` (nividia ifx stand-in, S1) |
| R | Timings, transfers | stale (S2) | `timing_*` |

### Open items
1. **Default paired FAIL** (3c above). Not to be resolved by changing the criterion after the fact.
2. **Default path not run-to-run reproducible.** Same S4 binary twice on `hub12x12_cb_b10_p` (same seeds): `data.h5`
   differs in 20 of 38 datasets by ≤ 2.4e-11 relative per bin (no growth, no divergence); Precision Green mean
   6.184e-10 vs 6.177e-10; the S2 build gave 6.186e-10. LOG repeats bit for bit. This is the likely source of the
   β = 10 drift (6.186/6.191/6.189e-10 across S1/S2/S3), but the cause is **not localised**; prime suspect: the
   `k_colnorms` reduction feeding the `la_qrpp` pivot sort. The default self-test figure `hub12x12_cb_b10` GTT moved
   9.42e-12 (S3) → 1.06e-11 (S4) with a fault-only change.
3. **Stale `nv_def`.** The first manager re-run (`mgr207.txt`) used an `nv_def` built from S3 (10:13). Covered for
   6×6 by the S3/S4 bit identity and for 12×12 default by the stamped S4 re-run (2b); LOG 12×12 (2c) is not re-run.
4. **Lost raw data.** Cannot be recomputed: the precision runs of acceptance 3a/3e, the compare2 run (4), the 12×12 LOG
   and 6×6 same-seed runs of `mgr207.txt`, the Sonnet extras (2d), the block table, the GPU=none run directories, and
   the agents' own self-test, fault, build and ctest logs (overwritten by manager runs with the same file names).
   Surviving raw copies: the paired runs (`raw_207`, `runs/cpair_*`), the Opus `taudiag_*.gz`, and all `mgr207c`
   run directories.
5. **Mismatches M1–M10 of the evidence audit** (2026-10-29):
   - M1: `hub6x6_cb_symm` default GTT reported 3.2e-12, measured 6.258e-12.
   - M2: LOG maxima reported ≤ 1.2e-11 (dense) and ≤ 8.4e-12 (12×12), measured 1.296e-11 and 8.805e-12.
   - M3: default `hub12x12_cb_b10` GTT 9.422e-12 (S3) vs 1.057e-11 (S4), see item 2.
   - M4: β = 4 LOG s3 ratio reported 1.81, recomputed 1.80.
   - M5: "β = 4 s2 at 1.33" listed under the final code came from a pre-NS binary.
   - M6: CPU τ-mean span reported 6.7–7.2e-12, actual 6.7–7.3e-12 (GPU=none 4.4–4.7e-12 not mentioned).
   - M7: "Green mean and max unchanged by 207" contradicted by the β = 10 default drift (item 2).
   - M8: default rows `fin_def`/`fin_nv_def` differ in the 4th digit (τ mean 6.403 vs 6.404e-12).
   - M9: LOG s1 1-ulp reference rows differ between agent and manager runs (chain start 330 vs 331).
   - M10: the ticket referred to a Resolution that did not exist.
6. **Acceptance 6 (Intel)** blocked by the oneAPI version mix (see `docs/gpu/MACHINES.md` update in the process branch).
7. **Sonnet QUESTION 1** (the `main.F90` hooks): accepted by the owner pending review (decision 4); the 2026-10-29 review
   raised no objection to them.
8. **Sonnet QUESTION 2** (judging τ precision on the max at all): part of the gate decision on acceptance 3.
9. **206 code: Green max 3.1× the CPU max** at β = 4 default (2.27e-6 vs 7.2e-7; information only under decision 1,
   no paired Green test exists).
10. **Debug variables reach production runs:** `ALF_GPU_TAUSWAP`, `ALF_GPU_TAU_NS`, `ALF_GPU_NS`, `ALF_GPU_FAULT_STAB`
    are read without a self-test gate (verification-hardening ticket).
11. **Positive controls missing** for the same-seed, paired, precision-mean and compare2 gates (`TAUSWAP=8`,
    `TAU_NS=0` never run through the paired test).
