# 302b – 204–206 device paths on the Arc B580 (GPU=intel)

- **Phase / lane:** 3 (Intel lane)
- **Owner agent:** fortran-gpu-dev (Sonnet)
- **Review:** numerics-reviewer (code, acceptance file, scripts) · **Verification:** acceptance-verifier
- **Depends on:** 302a merged into `gpu-port` (ONEMKL backend of `alf_la`); 212 (env reader, deterministic sums; in
  `gpu-port`); prerequisites P1–P8 below
- **Base:** `TBD` = `gpu-port` after the merge of the harness branch `gpu-port-harness-302b` (H1/H2); 302a is merged
  (`890185c0`), 207 is not in the base
- **Branch:** `gpu-port-302b` (the 2026-10-30 draft branch `gpu-port-302b-b580` is superseded)
- **Acceptance files:** `testsuite/gpu_reference/acceptance/302b.json` (Intel part, `intel`) and `302b.nvidia.json`
  (`nividia` part, `acc_driver.py --part nvidia`). Both are drafted together in `302b.draft.json`; the nvidia part sits under
  the key `_part_nvidia` and is split into its own file at pre-registration.
- **Status:** DRAFT (2026-10-30, updated 2026-10-31 with decisions 20–23 and the manager's answers); not pre-registered.
  Remaining before pre-registration: review (c) of the harness branch, the P8 timing probe, review (b) of this file

## Goal
With `GPU=intel` (ifx 2026.2 OpenMP offload, spir64 JIT, plus the 302a oneMKL backend), the 204 wraps (checkerboard and
dense), the 205 WRAPGR/updates and the 206 stabilization (default and STABLOG) run on the B580. They pass the 206/212
gates against the ifx CPU oracle. The 207 τ pass stays off under `GPU=intel` until 302c. `GPU=none` and `GPU=nvidia` stay
unchanged.

## Context
- **Guards in the code at 3e23a682** (unchanged since d8cc74f7; grep `ALF_LA_ONEMKL` in `Prog/`):
  - `Prog/gpu_stab_mod.F90:139–142`, `gpu_stab_setup`: `#if defined(ALF_LA_ONEMKL)` → `GPU stabilization (206): INACTIVE …
    (ticket 302)` + `return` / `#else` … `#endif` (end of the subroutine). PLAN cites 137–140; 212 shifted the lines.
  - `Prog/gpu_state_mod.F90:113–115` (`use alf_la_mod, only: la_init, la_finalize`), `577–579` (`la_init` in
    `gpu_state_setup`) and `923–925` (`la_finalize` in `gpu_state_finalize`). Line 351 only picks the backend name for
    the 212 build line; it is not a guard.
  - `Prog/gpu_wrap_mod.F90:98–102`: `la_available = .false.` under ONEMKL. Because of it, `dense_left`/`dense_right` stop
    with "dense hopping needs alf_la (ticket 302 for GPU=intel)" (lines 427–430, 473–476).
  - Today under `GPU=intel` the 205 updates and the checkerboard wraps run on the device; stabilization, τ and
    measurements run on the host. Measured on `intel`, 2026-10-28 (ifx 2025.3, old driver): 6×6 β = 4 cb took 0.475 s per
    sweep with host stabilization (`~/alf-gpu/work/t206/runs/intel_hub6x6_cb/info`).
- **τ path.** The current base has no device τ code: 207 lives on `gpu-port-207-tau` (WIP, ported onto 212). With 206
  active and `Ltau = 1`, the 206 INTERIM path copies `udvst(NST)` D→H in `gpu_stab_exchange`, and the host `TAU_M` runs.
  That is the NVIDIA state that 206/212 accepted. On the 207 branch, `gpu_tau_setup` checks only `gpu_stab_on`, `Ltau`,
  Projector and reconstruction, and has **no ONEMKL guard**. So once 207 is merged, removing the 302 guard would also turn
  on the device τ pass on the B580 (PLAN, "τ path on Intel").
  **The τ switch:** a compile-time guard at the top of `gpu_tau_setup`:
  `#if defined(ALF_LA_ONEMKL)` → info line ` GPU tau pass (207): HOST under GPU=intel (device tau pass refused until ticket 302c)` → `return` / `#endif`.
  `gpu_stab_tau_dev` stays `.false.`, so the host τ path runs as in 206.
  - The line must not contain the word `INACTIVE`: `verdict.py` treats that word as an inactive GPU path, and every run
    would become an INPUT_ERROR. If 207 is in the base, this line is added to the `active_lines`.
  - If 207 lands **after** 302b, the same guard must come with the 207 merge (the manager's merge checklist). Otherwise
    the device τ pass turns on under `GPU=intel` untested.
  - **Decided (owner, decision 22, 2026-10-31): "CPU τ pass".** Device τ off, host `TAU_M` on, the run continues.
    The same-seed τ datasets and Precision tau therefore check the host τ pass fed by the device stabilization.
  - **Manager (2026-10-31):** 207 is not in the base, so 302b touches no `gpu_tau_mod.F90`. The guard above is an item
    of the 207 merge checklist (in the 207 ticket); the NVIDIA lane adds it when 207 merges.
- **ALF_LA_TEST_\* refusal.** 302a's ticket says that 302b extends the 212 production refusal to `ALF_LA_TEST_*`. The 212
  reader (`gpu_env_read`, `gpu_state_mod.F90` ~200–235) collects only `ALF_GPU_*` entries. 302b extends it to the prefix
  `ALF_LA_TEST_`:
  - in a production run: ` GPU debug env refused: ALF_LA_TEST_<X> (production run; set ALF_GPU_DEBUG=1 for a debug run)`, exit 1;
  - in a debug run: listed in the sorted ` GPU debug env:` line.
  - This touches the shared part of `gpu_state_mod.F90` for every backend (CUDA would refuse these names too), so the
    manager sequences it.
  - **Manager (2026-10-31): all GPU backends.** The reader sits inside `#ifdef ALF_GPU` (`gpu_state_mod.F90:100`), so
    `GPU=none` is untouched; one rule for every GPU build is simpler to check than a backend-specific one, and
    `cuda_identity` runs without `ALF_LA_TEST_*`, so it stays bit-identical.
- **302a findings that matter here.** oneMKL `zgeqrf`/`zungqr`/`zunmqr` are hybrid, with a fixed cost of ~0.45 ms per call.
  STABLOG keeps host ZGEQP3 (`la_geqp3_host`). The BLAS calls run on the GPU. The interop queue is the OpenMP runtime's
  SYCL queue. `la_init` prints ` alf_la: ONEMKL device: Intel(R) Arc(TM) B580 Graphics` (seen in the 302a T2 logs).
- **Toolchain** (decision 13): `. ~/alf-gpu/env/intel.sh`, ifx 2026.2.0, MKL 2026.2, driver 26.35. No package hold: before
  every T1/T2 the manager records `dpkg-query` into the run directory and stops on any difference from
  `intel:~/alf-gpu/backup/intel_versions_2026-10-30.txt`. `intel.sh` does **not** set `ONEAPI_DEVICE_SELECTOR`.

## Files allowed
- `Prog/gpu_stab_mod.F90`: remove the ONEMKL guard in `gpu_stab_setup` (4 preprocessor/statement lines). Nothing else.
- `Prog/gpu_state_mod.F90`: the three ONEMKL guards (113–115, 577–579, 923–925). With the manager's sequencing, also the
  `ALF_LA_TEST_` prefix in the env reader and the allow-list (shared part).
- `Prog/gpu_wrap_mod.F90`: `la_available` (98–102). The stale "(ticket 302 for GPU=intel)" message text may stay.
- `Prog/gpu_tau_mod.F90`: does not exist in the base (207 not merged); not touched by 302b.
- `docs/gpu/scripts/intel/t302b_*` (development only).
- Not allowed:
  - `Libraries/Modules/alf_la_mod.F90`, `alf_la_sycl_shim.cpp`, `configure.sh`, Makefiles, `testsuite/**`;
  - `Prog/main.F90`, `Prog/tau_m_mod.F90`, `Prog/Wrapgr_mod.F90`; any CPU-path file;
  - the CUDA code paths.
- How the file list is enforced: GPU=none identity and CUDA identity cover CPU and CUDA behaviour on 6×6 cb only;
  `acc_gpunone.sh` skips every `Prog/gpu_*` file. Edits elsewhere in `Prog/gpu_*` (e.g. `gpu_update_mod.F90`, the CUDA
  dense path, shared self-test code) would not be caught by a gate. So the manager checks the candidate's
  `changed_files.txt` (written by the driver) against this list before T2, and review (a) checks the diff line by line
  (review (b) 2026-10-31, finding 5).- A defect found in `alf_la` (ONEMKL) is reported, not fixed here: its fix would change this file's expected counts
  (re-registration). **Manager (2026-10-31): no `alf_la` edits in 302b**; a defect found there gets its own ticket.

## Acceptance
The draft is in `302b.draft.json`; `verdict.py` lint passes for both parts in DRAFT mode. Tolerances are as in 212.
**Machine `intel`:**
- 7 builds: `in_def`, `in_def_pp` (`-DQR_PREPIVOT`), `in_log` (GPU=intel); `in_none`, `in_none_pp`, `in_none_log`,
  `base_in_none` (ifx CPU).
- 10 configs, identical to 212's, with their parameter hashes. The driver refuses if pyALF on `intel` prepares different
  parameters.
- 67 GPU runs, 10 CPU runs, 4 ctests, 15 gates, 18 controls.

Every `GPU=intel` run gets `ONEAPI_DEVICE_SELECTOR=level_zero:0` from `acc_run.sh` itself (H2; recorded in the
provenance env), and every production run must show the B580 device line as an active line (the refusal and
invalid-value runs stop in `gpu_env_check` before `la_init` and declare their own lines). The UHD 730 has no FP64.
The device τ pass must stay off (decision 22): `forbidden_lines` holds the 207 ACTIVE line, so a run that prints it is
an INPUT_ERROR.

| Gate | Oracle / statistic | Threshold | Control (expect) |
|---|---|---|---|
| `selftest_wrap` | 203 mirror + 204 wrap self-tests, `SELFTEST=2`, 6×6 cb and dense (first dense `la_gemm`/`la_gemm_cr` in ALF) | PASSED lines; numbers ≤ 1e-10 | `FAULT=hop`, `FAULT=vexp` (FAIL) |
| `selftest_update` | 205 `UPDTEST=1`, cb and flux | G ≤ 1e-10, phase ≤ 1e-12 | `FAULT=upd`, `FAULT=nsig` (FAIL) |
| `selftest_stab` | 206 `SELFTEST=3`, 5 configs × (default `_pp`, LOG) | U, G ≤ 1e-10; precision reduction ≤ 1e-12 | perm / eps:1e-9 (FAIL); eps:1e-11 (PASS) |
| `faults` | 8 stab faults (perm, phase, nan, eps:1e-9 × 2 builds) + hop, vexp, upd, nsig | rc 1, named FAILED line, armed line | eps:1e-11 armed but undetectable (FAIL) |
| `faults_mirror` | `SELFTEST=fault\|fault_v\|nan` | rc 1, FAILED naming hopping / Op_V / non-finite | clean st2 run (FAIL) |
| `refusal` | 15 `ALF_GPU_*` names of 212, production run | exact 212 line, rc 1 | `DEBUG=1` with NS, TL accepted (FAIL) |
| `refusal_la_test` | `ALF_LA_TEST_BREAK_ORDER`, `_POISON_INFO`, `_WRONG_LD` | same line format, rc 1 | `DEBUG=1` + POISON_INFO accepted (FAIL); harness change H1 |
| `invalid_value` | `FAULT_STAB=xyz`, `FAULT=xyz` | exact line, rc 1 | valid eps:1e-11 (FAIL) |
| `same_seed` | B580 vs ifx CPU, same seeds: 6×6 LOG / CPU LOG, 6×6 default / CPU `QR_PREPIVOT`, 12×12 default / CPU `QR_PREPIVOT` | no bin > 1e-4; overall and τ ≤ 1e-9; 38 datasets, 5 τ | `scale 1e-3` (FAIL), `nan` (INPUT_ERROR) |
| `determinism` | same binary twice: 6×6 default, 6×6 LOG, 12×12 β = 10 default, 12×12 β = 10 LOG | bit-identical `data.h5` + Precision lines | `scale 1e-12` (FAIL) |
| `precision` | 12×12 β = 4, 20 × 20: GPU mean / larger ifx CPU mean, Green and τ | ≤ 1.5 | `ALF_GPU_NS=0` (FAIL; see Q-NS) |
| `compare2` | `s8_U4_b2_cb_short` 120 × 30, B580 default vs ifx GPU=none, fresh seeds | m = 11, α = 0.01 Šidák, `min_bins` 8 | shift 8 (FAIL; decision 23); CPU-vs-CPU null (PASS) |
| `ifx_none_identity` | ifx 2026.2 GPU=none, base vs candidate, 6×6 | bit-identical | `scale 1e-12` (FAIL) |
| `alf_la_tests_onemkl` | ctest `^t_` under ONEMKL on the merged base | 100 % of 14; backend and B580 lines | `ALF_LA_TEST_WRONG_LD` (FAIL), replaces all 3 checks |
| `t_determinism_onemkl` | 212's micro-test under ONEMKL (also run in 302a's re-acceptance; regression check here) | 1/1 pass | `ALF_LA_TEST_PLANT_ULP=1` (FAIL) |

**Machine `nividia`** (`_part_nvidia`):

| Gate | Statistic | Threshold | Control |
|---|---|---|---|
| `gpu_none_ppcmp` | `acc_gpunone.sh` counts (only `Prog/gpu_*` changes; excluded from the comparison) | `0 0 0 0 0 0` | plant on `Prog/main.F90` (FAIL) |
| `gpu_none_ctest` | GNU devel ctest | 100 % of 53 (confirmed, decision 18) | one planted failure (FAIL) |
| `cuda_identity` | CUDA base vs candidate builds, same seeds, default and LOG, 6×6 | bit-identical | `scale 1e-12` (FAIL) |
| `none_identity` | nvfortran GPU=none base vs candidate | bit-identical | `scale 1e-12` (FAIL) |

CUDA "unchanged" is checked at run time, because `acc_gpunone.sh` does not preprocess `Prog/gpu_*` under the CUDA macros.
"HOST unchanged" is read as GPU=none (the alf_la HOST backend is used only by the tests).

### ifx/oneMKL-specific risks and the gate that covers each
- **Reduction order and determinism:**
  - Since 212 the floating-point sums are in fixed order. Only max reductions and integer counts remain, plus the setup
    checksum compared with a tolerance.
  - The host part of oneMKL's hybrid QR runs with `MKL_NUM_THREADS=1` (`acc_run.sh`, and `acc_ctest.sh` since the
    harness branch `gpu-port-harness-302b`).
  - Covered by `determinism` (both paths, β = 10) and `t_determinism_onemkl`.
  - Not covered: MKL with more than one thread on PVC (conditional numerical reproducibility, CNR) → 304.
- **Level Zero memory and ordering:**
  - The OpenMP target regions and the shim share one queue. The per-site D→H/H→D `target update`s are covered by
    `selftest_update` + nsig, `same_seed` and the self-tests.
  - 302a's `t_order` runs once inside `alf_la_tests_onemkl` (one repetition, not 302a's 100-repetition canary).
- **FP64 on Xe2 and fast math:**
  - `configure.sh` IntelLLVM uses `-fp-model=fast=2 -no-prec-div -xHost`. These flags reach `ALF_FLAGS_PROG`, so device
    code gets them too.
  - The CPU oracle has the same flags (common mode). The self-tests (1e-10 against the host routines), `same_seed` (1e-9)
    and `precision` cover the GPU-vs-CPU part. A bias common to ifx's flags is not covered here: GNU vs ifx was checked
    statistically in 106.
- **Device selection:** the selector in every GPU run, the B580 active line, and la_init's check that the OpenMP device
  equals the queue device (302a).
- **JIT vs AOT:** 302b stays JIT (`spir64`, what `configure.sh` builds). AOT is a build-path question for PVC (301/303) and
  saves ~1 s per process at most (302a profile).
  - JIT output depends on the IGC/driver version: covered by the dpkg record (decision 13).
  - `intel.sh` sets `NEO_CACHE_PERSISTENT=0`, so every process JIT-compiles its kernels afresh; replay pairs therefore
    also check that two JIT compilations give the same result.
- **Host fallback:** `OMP_TARGET_OFFLOAD=MANDATORY` (harness).

### Run lengths, time and budget
Run lengths are those of 212 (statistics, same hashes). No time limit applies (decision 16). `timeout_s` 3600 is a hang
guard only; the driver has a single global timeout.

**How the time was estimated:**
- A30 wall times come from 212's T2 driver log (`nividia:~/alf-gpu/work/t212/cand/mgr/8022c4b0de64/logs/`): 6×6 fd 70 s,
  12×12 fd 473 s, 12×12 β = 10 160 sweeps ~600 s, 12×12 β = 4 400 sweeps ~600 s, c2 1073 s, self-tests ~1 s.
- The B580/A30 factor f is ~1.35 at 6×6: 0.475 s/sweep (B580, host stabilization, 2026-10-28) vs 0.35 s/sweep (A30).
- Both GPUs look latency-bound at these sizes. 302a's per-call costs add ~0.1–0.4 s per sweep for the device
  stabilization (~1–5 ms per decompose/CGR, ~30 calls per sweep). FP64 GEMMs are 5–8× slower on the B580.
- Assumed: f = 1.5 central, 3 pessimistic.

| Block (GPU runs) | A30 (212) | B580 central | B580 pessimistic |
|---|---|---|---|
| 51 short debug/refusal/fault runs (JIT per process) | ~1 min | ~9 min | ~15 min |
| acc_ns, acc_tl, acc_la_poison | 3.5 min | 5 min | 11 min |
| same-seed + determinism 6×6 (4 runs) | 4.7 min | 7 min | 14 min |
| ss_def_12 | 7.9 min | 12 min | 24 min |
| determinism 12×12 β = 10 (4 runs) | 40 min | 60 min | 120 min |
| precision (3 runs) | 30 min | 45 min | 90 min |
| c2_gpu | 18 min | 27 min | 54 min |
| ctests (4) | – | ~5 min | ~8 min |
| **per tier (T1 or T2)** | **~104 min** | **~170 GPU-min** | **~335 GPU-min** |

- The nividia part takes ~5 GPU-min per tier (4 CUDA runs) plus CPU builds and ctest.
- The CPU oracle runs (10) run alongside the GPU runs on `intel` (`cpu_parallel` 2) and add no GPU time.
- **Budget: decision 20** (2026-10-31): no fixed GPU budget on the B580. The P8 timing probe still runs before
  pre-registration, to set `timeout_s` (hang guard) and to plan the T1/T2 slots; no cuts are planned.

### compare2
- Same config and length as 212's T2: 120 × 30. With decision 16 the run length needs no exception any more.
- The positive control plants 8 SEs: **decision 23** (2026-10-31, "Approve shift 8"). At m = 11, α = 0.01 Šidák
  (z ≈ 3.32, t ≈ k/√2), a shift of 5 is caught ~59 % of the time, 8 ~99 %.

### Other open points
- **Q-NS (owner decision pending, review (b) 2026-10-31, finding 1):** `ALF_GPU_NS=0` exists because of cuSOLVER's Q
  orthogonality (A30: CPU mean 4.1e-11, with NS 5.4e-11, without NS 9.7e-11, 2.4× the CPU). oneMKL's hybrid QR may give a
  LAPACK-quality Q, so the control may not reach 1.5× on the B580; then the precision gate has no killed control and T2
  cannot close. T0 runs after pre-registration, so "stop before pre-registration" cannot work. The code allowed in 302b
  has no other knob that degrades precision in production (`eps:` acts only under `SELFTEST=3`), and `mutate.py` cannot
  edit the Precision lines; a fallback control needs a harness or code change and the owner's approval.
- 12×12 LOG same-seed on the B580: left out; it belongs to ticket 213 (decision 3b). Adding it would need an owner
  decision naming 302b.
- Counts: `alf_la_tests_onemkl` 14 (302a's 13 + `t_determinism`) and GNU ctest 53 are confirmed by 302a's
  re-acceptance on the 212 base (decision 18); `gpu-port` added no test since.

## Prerequisites (manager) — state 2026-10-31
- **P1 done.** venv `intel:~/alf-gpu/venv` (decision 21): pyALF 2.6.post42, numpy, scipy, h5py as on `nividia`. The
  driver starts with the venv python.
- **P2 done.** Prebuilt HDF5 for ifx 2026.2 in `~/alf-gpu/hdf5/`.
- **P3 H1 implemented** on branch `gpu-port-harness-302b` (`05e677c4`): `acc_run.sh` accepts `ALF_LA_TEST_*` only with
  `--allow-debug` and records them as debug env; `verdict.py` `DEBUG_PREFIXES`, the gate lint and the driver's ctest
  check include the prefix; VERIFICATION §2 lists it. Review (c): changes required (selector evidence check, selfcheck
  cases); fixes in progress, then a fresh review (c).
- **P4 H2 implemented** on the same branch: `acc_run.sh` sets `ONEAPI_DEVICE_SELECTOR=level_zero:0` for `GPU=intel`
  builds and refuses overrides; the draft no longer declares it. **H3 done** outside the harness: `intel.sh` sets
  `NEO_CACHE_PERSISTENT=0` (2026-10-30); the old `~/.cache/neo_compiler_cache` waits for the owner.
- **P5 done.** 302a merged (`890185c0`); counts confirmed. The base is fixed when the harness branch merges.
- **P6** dpkg record before T1/T2 (decision 13), against `intel:~/alf-gpu/backup/intel_versions_2026-10-30.txt`.
- **P7 done.** Decisions 20 (budget), 22 (τ), 23 (shift 8).
- **P8 done** (2026-10-31, test-runner, `intel:~/alf-gpu/work/t302b/probe/`): base `GPU=intel` build of `3e23a682`
  (206 INACTIVE, stabilization on the host; 205 updates and cb wraps on the B580) vs the ifx CPU build, Ltau = 1, versions
  equal to the 2026-10-30 record. Not a 302b result.

  | config (probe length) | B580 build s/sweep | ifx CPU s/sweep | A30 212 s/sweep (device stab.) |
  |---|---|---|---|
  | hub6x6_cb β = 4 (3 × 20) | 0.141 | 0.027 | ~0.35 |
  | hub12x12_cb β = 4 (2 × 20) | 0.688 | 0.340 | ~1.5 |
  | hub12x12_cb β = 10 (2 × 10) | 1.687 | 0.852 | – |
  | s8_U4_b2_cb (3 × 30) | 0.127 | 0.035 | – |

  - Wall per sweep, from `date` around the run; info "CPU Time" agrees within ~10 %. GPU mirror setup 0.54 s per run.
  - The 2026-10-28 point (0.475 s/sweep at 6×6, ifx 2025.3, old driver) is obsolete: the same path now takes 0.141 s.
  - The time estimate above (f = 1.5–3 vs the A30) is pessimistic for the 204/205 part (measured f ≈ 0.4–0.5 with host
    stabilization). The device stabilization adds an unknown cost (302a: ~0.45 ms fixed per oneMKL QR call), so the
    estimate table stays as the upper bound; `timeout_s` 3600 stays a generous hang guard.
  - On the B580 the GPU build is 2–5× slower than one ifx CPU core at these sizes (information; no performance gate).
  - `alf_la: ONEMKL device:` is not printed on the base: `la_init` is skipped under the 302 guard
    (`gpu_state_mod.F90:577–579`). 302b removes that guard, so the active line in the draft is expected from 302b on;
    T0 must show it.

## Approvals
Empty (draft).

## Evidence
- `intel:~/alf-gpu/work/t302b/`: T0 in `dev/`; candidates in `cand/agent/<sha12>/` (T1) and `cand/mgr/<sha12>/` (T2).
- `nividia:~/alf-gpu/work/t302b/`: the nvidia part, same layout.
- Never delete raw outputs or binaries before the owner's gate; ≤ 1.5 GB of live evidence per machine.
- Every log carries the provenance stamp. Development scripts are `docs/gpu/scripts/intel/t302b_*`.
- Report: diff summary, transfers (the PLAN table is unchanged apart from the backend), s/sweep on the B580 from the run
  walls (information only).

## Handoff notes
(implementer)

## Out of scope
- The device τ pass on the B580 (302c).
- PVC, AOT, the SuperMUC profile and Slurm (301, 303, 304).
- Host vs device thresholds for small linear algebra (306).
- Performance work.
- Fixes to `alf_la`/the shim (see Files allowed).
- 213 (STABLOG same-seed excess).
- MKL threading/CNR (conditional numerical reproducibility) on PVC.
