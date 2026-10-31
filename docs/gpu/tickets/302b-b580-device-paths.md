# 302b – 204–206 device paths on the Arc B580 (GPU=intel)

- **Phase / lane:** 3 (Intel lane)
- **Owner agent:** fortran-gpu-dev (Sonnet)
- **Review:** numerics-reviewer (code, acceptance file, scripts) · **Verification:** acceptance-verifier
- **Depends on:** 302a merged into `gpu-port` (ONEMKL backend of `alf_la`); 212 (env reader, deterministic sums; in
  `gpu-port`); prerequisites P1–P8 below
- **Base:** `TBD-after-302a-merge` (gpu-port after the 302a merge; after the 207 merge too if 207 lands first)
- **Branch:** `gpu-port-302b-b580`
- **Acceptance files:** `testsuite/gpu_reference/acceptance/302b.json` (Intel part, `intel`) and `302b.nvidia.json`
  (`nividia` part, `acc_driver.py --part nvidia`). Both are drafted together in `302b.draft.json`; the nvidia part sits under
  the key `_part_nvidia` and is split into its own file at pre-registration.
- **Status:** DRAFT (2026-10-30); not pre-registered; open questions are marked `QUESTION:`

## Goal
With `GPU=intel` (ifx 2026.2 OpenMP offload, spir64 JIT, plus the 302a oneMKL backend), the 204 wraps (checkerboard and
dense), the 205 WRAPGR/updates and the 206 stabilization (default and STABLOG) run on the B580. They pass the 206/212
gates against the ifx CPU oracle. The 207 τ pass stays off under `GPU=intel` until 302c. `GPU=none` and `GPU=nvidia` stay
unchanged.

## Context
- **Guards in the code at d8cc74f7** (grep `ALF_LA_ONEMKL` in `Prog/`):
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
  - `QUESTION:` (owner) does "τ path refused at setup" mean this (device τ off, host τ on, the run continues), or should
    a `GPU=intel` run with `Ltau = 1` stop? Stopping would force `Ltau = 0` configs: new parameter hashes, and no τ
    datasets in the same-seed gate. The draft assumes the host τ fallback.
  - `QUESTION:` (manager) `gpu_tau_mod.F90` is NVIDIA-lane-owned (PLAN ownership table). A one-line guard there needs
    your sequencing.
- **ALF_LA_TEST_\* refusal.** 302a's ticket says that 302b extends the 212 production refusal to `ALF_LA_TEST_*`. The 212
  reader (`gpu_env_read`, `gpu_state_mod.F90` ~200–235) collects only `ALF_GPU_*` entries. 302b extends it to the prefix
  `ALF_LA_TEST_`:
  - in a production run: ` GPU debug env refused: ALF_LA_TEST_<X> (production run; set ALF_GPU_DEBUG=1 for a debug run)`, exit 1;
  - in a debug run: listed in the sorted ` GPU debug env:` line.
  - This touches the shared part of `gpu_state_mod.F90` for every backend (CUDA would refuse these names too), so the
    manager sequences it.
  - `QUESTION:` (manager) accept this design, or should the refusal stay ONEMKL-only?
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
- `Prog/gpu_tau_mod.F90`: the ONEMKL guard in `gpu_tau_setup`, only if 207 is in the base, and with the manager's
  sequencing.
- `docs/gpu/scripts/intel/t302b_*` (development only).
- Not allowed:
  - `Libraries/Modules/alf_la_mod.F90`, `alf_la_sycl_shim.cpp`, `configure.sh`, Makefiles, `testsuite/**`;
  - `Prog/main.F90`, `Prog/tau_m_mod.F90`, `Prog/Wrapgr_mod.F90`; any CPU-path file;
  - the CUDA code paths.
- A defect found in `alf_la` (ONEMKL) is reported, not fixed here: its fix would change this file's expected counts
  (re-registration). `QUESTION:` (manager) should 302b instead be allowed ONEMKL-only edits of `alf_la_mod.F90` and the
  shim, with the counts of 302a's nvidia part?

## Acceptance
The draft is in `302b.draft.json`; `verdict.py` lint passes for both parts in DRAFT mode. Tolerances are as in 212.
**Machine `intel`:**
- 7 builds: `in_def`, `in_def_pp` (`-DQR_PREPIVOT`), `in_log` (GPU=intel); `in_none`, `in_none_pp`, `in_none_log`,
  `base_in_none` (ifx CPU).
- 10 configs, identical to 212's, with their parameter hashes. The driver refuses if pyALF on `intel` prepares different
  parameters.
- 67 GPU runs, 10 CPU runs, 4 ctests, 15 gates, 18 controls.

Every GPU run declares `ONEAPI_DEVICE_SELECTOR=level_zero:0`, and every run must show the B580 device line as an active
line. The UHD 730 has no FP64.

| Gate | Oracle / statistic | Threshold | Control (expect) |
|---|---|---|---|
| `selftest_wrap` | 203 mirror + 204 wrap self-tests, `SELFTEST=2`, 6×6 cb and dense (first dense `la_gemm`/`la_gemm_cr` in ALF) | PASSED lines; numbers ≤ 1e-10 | `FAULT=hop`, `FAULT=vexp` (FAIL) |
| `selftest_update` | 205 `UPDTEST=1`, cb and flux | G ≤ 1e-10, phase ≤ 1e-12 | `FAULT=upd`, `FAULT=nsig` (FAIL) |
| `selftest_stab` | 206 `SELFTEST=3`, 5 configs × (default `_pp`, LOG) | U, G ≤ 1e-10; precision reduction ≤ 1e-12 | perm / eps:1e-9 (FAIL); eps:1e-11 (PASS) |
| `faults` | 8 stab faults (perm, phase, nan, eps:1e-9 × 2 builds) + hop, vexp, upd, nsig | rc 1, named FAILED line, armed line | eps:1e-11 armed but undetectable (FAIL) |
| `faults_mirror` | `SELFTEST=fault\|fault_v\|nan` | rc 1, FAILED naming hopping / Op_V / non-finite | clean st2 run (FAIL) |
| `refusal` | 15 `ALF_GPU_*` names of 212, production run | exact 212 line, rc 1 | `DEBUG=1` with NS, TL accepted (FAIL) |
| `refusal_la_test` | `ALF_LA_TEST_BREAK_ORDER`, `_POISON_INFO`, `_WRONG_LD` | same line format, rc 1 | `DEBUG=1` + POISON_INFO accepted (FAIL); **needs harness change H1** |
| `invalid_value` | `FAULT_STAB=xyz`, `FAULT=xyz` | exact line, rc 1 | valid eps:1e-11 (FAIL) |
| `same_seed` | B580 vs ifx CPU, same seeds: 6×6 LOG / CPU LOG, 6×6 default / CPU `QR_PREPIVOT`, 12×12 default / CPU `QR_PREPIVOT` | no bin > 1e-4; overall and τ ≤ 1e-9; 38 datasets, 5 τ | `scale 1e-3` (FAIL), `nan` (INPUT_ERROR) |
| `determinism` | same binary twice: 6×6 default, 6×6 LOG, 12×12 β = 10 default, 12×12 β = 10 LOG | bit-identical `data.h5` + Precision lines | `scale 1e-12` (FAIL) |
| `precision` | 12×12 β = 4, 20 × 20: GPU mean / larger ifx CPU mean, Green and τ | ≤ 1.5 | `ALF_GPU_NS=0` (FAIL; see Q-NS) |
| `compare2` | `s8_U4_b2_cb_short` 120 × 30, B580 default vs ifx GPU=none, fresh seeds | m = 11, α = 0.01 Šidák, `min_bins` 8 | shift 8 (FAIL; **owner approval needed**); CPU-vs-CPU null (PASS) |
| `ifx_none_identity` | ifx 2026.2 GPU=none, base vs candidate, 6×6 | bit-identical | `scale 1e-12` (FAIL) |
| `alf_la_tests_onemkl` | ctest `^t_` under ONEMKL on the merged base | 100 % of 14; backend and B580 lines | `ALF_LA_TEST_WRONG_LD` (FAIL), replaces all 3 checks |
| `t_determinism_onemkl` | 212's micro-test under ONEMKL (never run there) | 1/1 pass | `ALF_LA_TEST_PLANT_ULP=1` (FAIL) |

**Machine `nividia`** (`_part_nvidia`):

| Gate | Statistic | Threshold | Control |
|---|---|---|---|
| `gpu_none_ppcmp` | `acc_gpunone.sh` counts (only `Prog/gpu_*` changes; excluded from the comparison) | `0 0 0 0 0 0` | plant on `Prog/main.F90` (FAIL) |
| `gpu_none_ctest` | GNU devel ctest | 100 % of 53 (count to be confirmed on the merged base) | one planted failure (FAIL) |
| `cuda_identity` | CUDA base vs candidate builds, same seeds, default and LOG, 6×6 | bit-identical | `scale 1e-12` (FAIL) |
| `none_identity` | nvfortran GPU=none base vs candidate | bit-identical | `scale 1e-12` (FAIL) |

CUDA "unchanged" is checked at run time, because `acc_gpunone.sh` does not preprocess `Prog/gpu_*` under the CUDA macros.
"HOST unchanged" is read as GPU=none (the alf_la HOST backend is used only by the tests).

### ifx/oneMKL-specific risks and the gate that covers each
- **Reduction order and determinism:**
  - Since 212 the floating-point sums are in fixed order. Only max reductions and integer counts remain, plus the setup
    checksum compared with a tolerance.
  - The host part of oneMKL's hybrid QR runs with `MKL_NUM_THREADS=1` (`acc_run.sh`).
  - Covered by `determinism` (both paths, β = 10) and `t_determinism_onemkl`.
  - Not covered: MKL with more than one thread on PVC (conditional numerical reproducibility, CNR) → 304.
- **Level Zero memory and ordering:**
  - The OpenMP target regions and the shim share one queue. The per-site D→H/H→D `target update`s are covered by
    `selftest_update` + nsig, `same_seed` and the self-tests.
  - 302a's `t_order` is kept as a regression in `alf_la_tests_onemkl`.
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
  - The NEO compiler cache (below) means that replay pairs run on cached kernels.
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
- **This exceeds decision 11** (T1 ≈ 90, T2 ≈ 120). Proposal for the owner: T1 ≈ 200 and T2 ≈ 220 GPU-min on the B580.
  The basis is a timing probe the manager runs **before pre-registration** (P8). If the probe gives f > 2, the manager
  returns with cuts. The recommended first cut: the β = 10 determinism config at 4 bins × 20 sweeps (saves ~30–60 min).
  Shorter precision runs are not recommended (their false-alarm rate is not calibrated on oneMKL).
- `QUESTION:` (owner) approve this budget, or prefer the β = 10 cut now?

### compare2
- Same config and length as 212's T2: 120 × 30. With decision 16 the run length needs no exception any more.
- The positive control plants 8 SEs. Decision 12 approved that for 212 only. At m = 11, α = 0.01 Šidák (z ≈ 3.32,
  t ≈ k/√2), a shift of 5 is caught ~59 % of the time, so 5 does not make the gate's control reliable.
- `QUESTION:` (owner) approve shift 8 for 302b (draft) — or keep 5 and accept ~41 % of T2s failing the control?

### Other open points
- **Q-NS:** `ALF_GPU_NS=0` exists because of cuSOLVER's Q orthogonality. oneMKL's Q may not need it, so the control may
  not reach 1.5×. T0 shows the ratio on dev seeds. If it stays ≤ 1.5, stop: a different precision control needs the
  owner's approval before pre-registration.
- `QUESTION:` (owner) gate 12×12 LOG same-seed on the B580? It is ticket 213's (decision 3b). A `report_only` gate needs a
  decision naming 302b. The draft leaves it out.
- `QUESTION:` (manager) the `alf_la_tests_onemkl` count of 14 (302a's 13 + `t_determinism`) depends on how the 302a merge
  resolves `testsuite/alf_la.tests/CMakeLists.txt`. Confirm it on the merged base. The same holds for GNU ctest 53.

## Prerequisites (manager)
- **P1** A venv on `intel` (`~/alf-gpu/venv`). It does not exist today, and system python has neither `h5py` nor
  `py_alf`. It needs the same pyALF version as nividia (otherwise the parameter hashes refuse), plus numpy, scipy and
  h5py. The driver must start with the venv python (212's `env -i` PATH issue). It is a user-space install:
  `QUESTION:` owner OK.
- **P2** Prebuilt HDF5 for ifx 2026.2 in `~/alf-gpu/hdf5/`. It is being built; assumed present.
- **P3** Harness change **H1**, with review (c) before pre-registration:
  - `acc_run.sh` accepts `ALF_LA_TEST_*` only with `--allow-debug`;
  - `verdict.py` `DEBUG_PREFIXES` adds `ALF_LA_TEST_`;
  - VERIFICATION §2 lists it.
  Without H1, drop `refusal_la_test` and its control.
- **P4** Optional harness changes, `QUESTION:` (manager):
  - **H2:** `acc_run.sh` sets `ONEAPI_DEVICE_SELECTOR=level_zero:0` for GPU=intel builds, as `acc_ctest.sh` already does.
    The draft instead declares it in every GPU run.
  - **H3:** keep the NEO compiler cache inside the run's tmp (`NEO_CACHE_DIR`) or switch it off (`NEO_CACHE_PERSISTENT=0`).
    Today every GPU process on `intel` writes `~/.cache/neo_compiler_cache` (16 MB on 2026-10-30, 302a included), which is
    outside `~/alf-gpu/` (confinement).
- **P5** 302a merged; the 302b base fixed; counts confirmed (above).
- **P6** dpkg record before T1/T2 (decision 13).
- **P7** Owner approvals: shift 8, budget, the τ-switch reading.
- **P8** Timing probe of the base `GPU=intel` build: 1 bin each of 6×6 fd, 12×12 β = 4, 12×12 β = 10, s8. This is not a
  302b result.

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
