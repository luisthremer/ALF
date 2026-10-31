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

## Handoff notes

### 2026-10-30, rebase onto 212 (implementer, T0 only; uncommitted working tree on `gpu-port-207-tau` over 93efa861)
- **State:** 207 rebased on 3f8ec8fa (206 + 212). Builds (nvhpc 26.9, default / QR_PREPIVOT / LOG; GPU=none GNU and nvhpc).
  Stopped at acceptance 3: the **default-path paired test still FAILS the sd condition** (dev1 seeds: sd ratio 1.353 > 1.25;
  mean −0.0021, counts 51/47 and 11/14 pass). Not tuned, not re-run. LOG passes (1.062).
- **Code decisions:**
  - every `ALF_GPU_*` read in `gpu_tau_mod` goes through `gpu_env_get`/`gpu_env_is`; all 207 variables are refused in
    production by 212's allow-list (no new code); `ALF_GPU_FAULT_TAU` validated in `gpu_env_check` (invalid → exit 1 in
    every run, armed line); graded fault `eps:<x>` multiplies the largest |GT0| element of the host copy of the device GT0
    in every device-vs-host GT0 comparison of the SELFTEST=4 self-test (fails at 1e-9, passes at 1e-11). Faults now act in
    any debug run where the faulted step runs (as 212's FAULT_STAB), not only under SELFTEST=4.
  - `G00_d` lives in `gpu_tau_mod`, not in `gpu_state_mod`: a fourth Ltau array in `gpu_state_mod` made nvfortran 26.9
    (`-fast`, STABLOG) crash compiling `gpu_stab_mod` (fort1 `malloc(): mismatching next->prev_size`, signal 6;
    reproducible; `-O0` compiles; the unchanged 212 `gpu_stab_mod` crashes too against the modified `gpu_state_mod`).
  - `ALF_GPU_TAUDIAG_SEED` not implemented (meaning unspecified).
- **Harness finding (blocks the paired gate):** 212's line ` GPU debug env: ..., ALF_GPU_TAUDIAG=2, ...` contains
  "TAUDIAG"; `verdict.py:taudiag()` (and `mutate.py`) select lines by `"TAUDIAG" in line` and raise INPUT_ERROR
  ("malformed TAUDIAG line") on it. Dev script `t207_diagstat.py` now selects `TAUDIAG pt ` lines only.
  `verdict.py` replay also compares only the first line of each `Precision` entry (the Max value is on the next line).
- **Evidence (T0, dev seeds, not acceptance):** nividia `~/alf-gpu/work/t207/dev2/` — `logs/t0a_1.txt` (self-tests,
  faults, refusals), `logs/t0b_1_eval.txt` (determinism, Ltau invariance, non-interference, same seed),
  `logs/t0c_1_paired.txt`, `logs/t0c_1_precision.txt`, `logs/t0c_1_compare2.txt`, `gpunone1/`, `logs/t0d_1.txt`;
  run directories `runs/t0a1_*`, `runs/t0b1_*`, `runs/t0c1_*`, `runs/t0c2_*`, `runs/t0d1_*`; binaries `bin/`.
  Scripts `docs/gpu/scripts/nividia/t207_d2_*`.

### 2026-10-31, dev5: accuracy test against a quad-precision reference (decision 19, option E): proposal for the owner
Development on dev seeds (dev1–dev4 = the dev4 seed sets), not a candidate, no acceptance file or harness script touched. Working
tree on `gpu-port-207-tau` over c1cc528e (uncommitted diff sha256 `eba7029c…`).

- **Code (debug aid only):** `Prog/gpu_tau_mod.F90` gets `ALF_GPU_TAUDIAG_QREF=<stride>` (+ `_QVAR=<masks>`, `_QHOST=1`; with
  `_FROM`): at every stride-th stabilization point a stream file `qref/q<k>.bin` with the inputs of the production CGR2_2 (device
  udvr, udvst(nst)), its four outputs, the outputs of `dev_cgr22` with the swp masks of `_QVAR` on the same inputs, the host
  CGR2_2 of the same binary (`_QHOST`), with TAUDIAG=1|2 the host chain's state, and the next block's device propagation B, B⁻¹
  (slice kernels applied to the identity). Diagnostic D→H copies only; production identity: `data.h5` of
  `hub12x12_cb_b10_s` (no debug variable) identical to c1cc528e, 38/38 datasets, default and LOG.
- **Quad reference** (`docs/gpu/scripts/nividia/t207_d5_qref.F90`, offline, gfortran 13.3 `real(16)`; **nvfortran 26.9 has no
  `real(16)`/`complex(16)`**, KIND error): CGR2_2 in binary128 with the same formulae, branch and STABLOG splitting; V1⁻¹ by LU
  with partial pivoting; Householder QR (zlarfg convention) after the column-norm pre-pivot (default) or with exact-norm column
  pivoting (LOG); D3 separated; extended system solved as `solve_extended_System`. Checks: quad vs itself with randomly
  row-permuted QR input and with the other branch: max 3.3e-29 / 3.1e-29 over 144 real points (synthetic easy states 2e-32);
  double host CGR2_2 vs quad on easy synthetic states 8e-15 – 2e-14. Cost 1.8–2.6 s per CGR2_2 at N = 144 (4 threads); a
  280-point run takes 12.5 (default) / 16 (LOG) CPU-minutes. The same program evaluates the double oracle (cgr2_2_mod objects
  of a GNU build, QR_PREPIVOT resp. STABLOG), 9 (default) / 7 (LOG) rounding-order nulls (copy of `nul_cgr22`, compiled with the
  ALF flags so that code 0 reproduces the oracle bit for bit; codes 4, 100004, 1, 100001, 2, 5, 100005, and 1608, 4808 on
  default), and two 1-ulp input perturbations.
- **Finding 1, pointwise accuracy does not see the planted defects.** Max-relative / Frobenius error of the four blocks against
  quad: the device is *more* accurate than the CPU oracle (mean log10 ratio, pooled: default emax −0.157, eF −0.159; LOG −0.039,
  −0.034), and `var8` (device ZUNMQR without the Newton–Schulz step = the `ALF_GPU_TAU_NS=0` / `ALF_GPU_TAUSWAP=8` defect) has
  the identical error (equal to 3 digits on every set). A pointwise accuracy criterion cannot kill either control.
- **Finding 2, the error as the next precision check receives it does.** pTT = max|B (G − Gq)_TT B⁻¹| (B: the next block's
  device propagation; ΔG formed in quad): CGR2_2's ~4e-13 pointwise error becomes ~2.5e-9 (default) / 5e-9 (LOG), the size of
  the Precision-tau numbers. Statistic S_ens = mean over points of log10(pTT_m / g), g = geometric mean of pTT over the correct
  host variants at that point (oracle, nvfortran host CGR2_2, nulls; leave-one-out for a member). Per seed set dev1..dev4
  (210 points each; generated tables `logs/report_agg_{def,log}_pTT.txt`):
  - default: device +0.084, +0.055, +0.034, +0.054 (pooled +0.057 [0.039, 0.074]); `var8` +0.280, +0.294, +0.307, +0.266;
    production control `ALF_GPU_TAU_NS=0` (dev1) +0.280; correct variants max per set +0.074, +0.076, +0.120, +0.073;
  - LOG: device −0.012, −0.030, −0.059, −0.048 (pooled −0.037 [−0.052, −0.021]); `var8` +0.215, +0.194, +0.214, +0.205;
    production control `ALF_GPU_TAUSWAP=8` (dev1) +0.215; correct variants max per set +0.036, +0.058, +0.050, +0.093.
  - The production controls give exactly the stand-in's numbers (the τ pass does not feed back into the CGR2_2 inputs).
  - Against the oracle alone (S_cpu) the device is +0.093…+0.119 (default) / −0.016…+0.022 (LOG); S_cpu separates LOG worse
    (correct variants up to +0.152, control down to +0.220), hence the ensemble reference.
- **206 input part** (TAUDIAG=2 runs, dev1, 140 points; `logs/report_checks.txt`): quad G of the device input vs quad G of the
  host chain's input differ pointwise by eF 3.4e-13 (geometric mean; the size of CGR2_2's rounding error) but carry only
  7.5e-13 / 7.8e-13 into the next check (pTT), against 2.5e-9 / 5.0e-9 for the CGR2_2 error; host CGR2_2 is not harder on the
  device input (mean log10 pTT ratio device input / host input +0.028 default, −0.024 LOG). The Precision-tau numbers are
  dominated by CGR2_2 rounding amplified by the next block, not by the 206 states. (Not measured: the inputs' own error against an
  exact WRAPUR, which would need a quad propagation.)
- **Proposed criterion (replaces the three paired conditions of acceptance 3 if the owner approves; decision 1's sd condition
  stays in force until then):** per path (default with the QR_PREPIVOT-variant-1 oracle, LOG with upstream), `hub12x12_cb_p`
  (β = 4, 400 sweeps), debug run with `ALF_GPU_TAUDIAG_QREF=2 ALF_GPU_TAUDIAG_QHOST=1 ALF_GPU_TAUDIAG_FROM=331` (280 points, 210
  with a next block), offline `t207_d5_qref` with the null codes above; **S_ens(pTT) pooled over two fresh seed sets ≤ 0.15 on
  each path**; positive controls `ALF_GPU_TAU_NS=0` (default) and `ALF_GPU_TAUSWAP=8` (LOG), one fresh seed set each, must give
  S_ens > 0.15; input checks (else INPUT_ERROR): every dump processed, quad self-check (row-permuted, other branch) ≤ 1e-25
  where evaluated, null code 0 identical to the oracle, 210 points with B. Pointwise errors and S_cpu reported as information.
- **Calibration of the bound 0.15** (`logs/report_agg_*`, `logs/report_pairs_015.txt`): margins per set (device max / correct
  max / control min): default +0.084 / +0.120 / +0.266, LOG −0.012 / +0.093 / +0.194; correct variants' bootstrap exceedance of
  0.15 per set: mean 0.005, worst 0.125 (default dev3 n100005) / mean 1e-4, worst 0.003 (LOG); pooled over any two of the four
  sets (6 pairs × 2 paths, 2000 bootstrap draws): P(device > 0.15) = 0, P(control ≤ 0.15) = 0, worst correct variant ≤ 0.001.
  Graded fault: the device's propagated error would have to grow by a factor 10^(0.15 − S) = 1.24 (default, pooled) / 1.54
  (LOG) to fail; the controls are factors 1.31 / 1.11 above the bound (per-set minimum). Thresholds chosen after seeing all
  dev data (same caveat as decision 1's thresholds): judge only on fresh seeds.
- **Old conditions:** sd condition: measures rounding-pattern similarity (dev4), not accuracy: replace. Mean condition: its
  defect separation was +0.050…+0.064 (default) against correct ≤ 0.032, a margin of 0.018 at the bound 0.05; the proposed
  statistic separates controls from correct variants by 0.146 (default) / 0.101 (LOG). I see no added detection power in
  keeping it as a gate; report-only is possible (it needs a TAUDIAG=2 run, 11–12 GPU-min per path). Count condition: tail
  information only.
- **Run length and GPU minutes:** QREF runs 8.0–9.7 min (default) / 9.4–10.4 min (LOG) on the A30, plus 12.5 / 16 CPU-min
  (4 threads, GPU idle) and 1.6 GB transient dumps per run. T2 as proposed (2 sets × 2 paths + 2 controls): ≈ 55 GPU-min,
  ≈ 85 CPU-min; T1 (one set per path + 2 controls) ≈ 37 GPU-min. Needs harness support (build of the offline program and the GNU
  oracle objects from the archive, dump processing, verdict from `qref.txt`): manager/owner.
- **Open:** the default-path device is at the top of the correct band (S_ens +0.034…+0.084; S_cpu ≈ +0.10, i.e. ~25 % larger
  propagated error than the oracle, like the blocked-QR null n4808 at +0.095); not attributed to a device step (QVAR host-twin
  masks 16/32/1/2/4/256 in one TAUDIAG run would do it). nvfortran 26.9 crashed once (fort2 signal 11) on `gpu_wrap_mod` /
  `gpu_stab_mod` in a parallel make; the rebuild passed.
- **Evidence:** nividia `~/alf-gpu/work/t207/dev5/`: runs `runs/{s,m,c,d}_*/hub12x12_cb_p/` (run.log, info, data.h5,
  provenance.json, stamp.txt, qref.txt, qref_proc.json, qref_manifest.txt with the sha256 of every dump; every 40th dump kept,
  smoke runs all), identity runs `runs/id_*`; generated reports `logs/report_*.txt|json`, chains `logs/chain_*.txt`; binaries
  `bin/` (GPU `nv_def_pp2` 0637587e…, `nv_log2` fd6022df…; base c1cc528e `base_nv_*`; offline `old_q/v2/qref_{pp,log}`
  d0c1f5bd…/9e17e710… for the main and control runs, `qref_{pp,log}` e989f36d…/b27c2643… for the din runs). GPU total 140.3
  min. Scripts `docs/gpu/scripts/nividia/t207_d5_*`.
