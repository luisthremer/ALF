# Build options and runtime switches of the GPU port

Reference page. Every entry was taken from the code (state of this worktree); "not determined" means the code does not say.
Line numbers are anchors and drift; search for the name. Categories:

- **production**: meant for production runs.
- **tuning**: changes speed (and possibly rounding), not the physics; needs an acceptance check before production.
- **diagnostic**: prints or measures, does not change the result; slows the run.
- **test-only**: self-tests, fault injection, negative controls. Never in a production run.

## 0. The production guard (read this first)

A GPU build calls `gpu_env_check` at setup (`Prog/gpu_state_mod.F90:360`, called from `gpu_state_setup`, line 635).

- It reads `/proc/self/environ` (Linux only; if unreadable the run stops, line 217) and collects every variable named
  `ALF_GPU_<x>` or `ALF_LA_TEST_<x>` with a non-empty value. The value `0` counts as set.
- It prints one line `GPU debug env: ...` (or `none`) and one line `GPU build: git <hash> backend <CUDA|ONEMKL|HIP|HOST>` to stdout and `info`.
- **A run is a debug run** only if `ALF_GPU_SELFTEST`, `ALF_GPU_UPDTEST` or `ALF_GPU_DEBUG` is set to a value other than `0` (line 409).
- **In a production run every `ALF_GPU_*` and `ALF_LA_TEST_*` variable except `ALF_GPU_DISABLE` and `ALF_GPU_DEBUG` stops the run** with
  `GPU debug env refused: <name>` (line 415). So all switches in the tables below named `ALF_GPU_*` or `ALF_LA_TEST_*` need
  `ALF_GPU_DEBUG=1` (or SELFTEST/UPDTEST). `ALF_GPU_DEBUG=1` itself does nothing else.
- Consequence: a leaked `ALF_GPU_CC` (a build variable, section 1) in the job environment also stops a production run. `configure.sh`
  exports it for `HELMAGPU`. Not tested; follows from the code. `ALF_GPU` (no suffix) is not matched.
- `ALFSH_*`, `ALF_HIP_*`, `ALFDEVQR_*`, `ALF_UNSTABLE_THRESHOLD` are outside the guard. They are never refused, so they must be checked by hand.
- Fault variables with an unknown value stop every run, debug or not (lines 423-472).

## 1. Build options (`configure.sh`)

Usage: `. ./configure.sh MACHINE MODE STAB [HDF5] [Devel] [GPU=none|nvidia|intel|amd]` (`configure.sh:5`, `:267`). An unknown `GPU=` value or a wrong
MACHINE is an error (no fallback, `:546`). `GPU=none` or no argument leaves the CPU build unchanged.

| Name | Values / default | Effect | Category | Where |
|---|---|---|---|---|
| `GPU=none` | default | Host build, upstream behaviour. No GPU macro. `alf_la` uses plain BLAS/LAPACK (`ALF_LA_HOST` semantics). | production | `configure.sh:644` |
| `GPU=nvidia` | needs MACHINE `PGI` or `HELMAGPU` | nvfortran, `-mp=gpu -gpu=<cc>`, `-cudalib=cublas,cusolver`; defines `ALF_GPU`, `ALF_LA_CUDA`; exports `ALF_GPU=nvidia`. | production | `configure.sh:612-620` |
| `GPU=intel` | needs MACHINE `IntelLLVM` or `IntelX` | ifx, OpenMP offload `spir64` (JIT), oneMKL + SYCL shim; defines `ALF_GPU`, `ALF_LA_ONEMKL`; removes `-static`; exports `ALF_GPU=intel`. | production | `configure.sh:635-641` |
| `GPU=amd` | needs MACHINE `GNU` (GCC >= 14 with amdgcn offload) | gfortran `-foffload=amdgcn-amdhsa`, hipBLAS/hipSOLVER, own HIP kernels; defines `ALF_GPU`, `ALF_LA_HIP`; removes `-static` and `-ffast-math`; exports `ALF_GPU=amd`. | production | `configure.sh:622-632` |
| `ALF_GPU_CC` | env; default `cc80`; `HELMAGPU`: `cc90,cuda12.9` (nvhpc 25.7) or `cc90,cuda12.6` | Target compute capability for `GPU=nvidia`, e.g. `cc80,cc90`. Build time only; do not leave it exported at run time (section 0). | production | `configure.sh:612`, `:506-508` |
| `ALF_HELMA_NVHPC` | env; `25.7` (default) or `24.11` | MACHINE `HELMAGPU`: picks the nvhpc and OpenMPI modules. 24.11 builds but the L=8 run aborts with `la_qrpp: NaN or Inf` (measured 2026-11-11); do not use. | production | `configure.sh:502-506` |
| `ALF_AMD_ARCH` | env/make; default `gfx90a` (MI200) | AMD GPU architecture for gfortran offload and `hipcc --offload-arch`. | production | `configure.sh:625`, `Libraries/Modules/Makefile:59` |
| `ROCM_PATH` | env; default `/opt/rocm` | ROCm prefix (include/, lib/, bin/hipcc) for the HIP shim and link line. | production | `configure.sh:629`, `Makefile:57,60` |
| `HIPCC` | make variable; default `$ROCM_PATH/bin/hipcc` | Compiler of `alf_hip_kernels.hip.cpp` and `alf_hip_qr.hip.cpp`. | production | `Libraries/Modules/Makefile:60` |
| `ALF_SYCL_BLK_FLAGS` | make/env; default empty (JIT) | Extra `icpx` flags for `libalfsb.so` (async block kernels) and `libalfdevqr.so`. For PVC ahead-of-time: `-fsycl-targets=spir64_gen -Xs "-device pvc"`. | production | `Libraries/Modules/Makefile:48,52` |
| `MKLROOT` | set by the oneAPI environment | oneMKL include/lib path for `GPU=intel`. | production | `configure.sh:636`, `Makefile:43` |
| `ALF_FLAGS_EXT` | env; default empty | Extra compiler flags appended by `configure.sh`. This is how the optional `-D` macros of section 2 are set. | production | `configure.sh:38,592` |
| `ALF_HDF5_DIR` | env; default self-built `HDF5/<compiler>` | Use an existing HDF5 install instead of building one. | production | `configure.sh:69` |
| `libalfdevqr.so` (make target `lib`) | built only for `GPU=intel` | Device QR library, loaded at run time by `dlopen` only if `ALFSH_DEVQR`/`ALFSH_DEVUNMQR` is set. The default build does not need it. | tuning | `Libraries/Modules/Makefile:50-52` |

## 2. Preprocessor macros

| Name | Values / default | Effect | Category | Where |
|---|---|---|---|---|
| `ALF_GPU` | set by `configure.sh` for any GPU target | Compiles the GPU modules (`gpu_state/wrap/lowrank/update/stab/tau/obs_mod`) and the hooks in `main.F90`, `Wrapgr_mod.F90`, `tau_m_mod.F90`, `Hamiltonian_Hubbard_smod.F90` (there is no GPU code in `upgrade_mod.F90`; it only carries the `ALF_TRACE_UPDATES` hook). Unset (GPU=none): none of this is compiled. | production | `Prog/Makefile:17`, `Prog/main.F90:710` |
| `ALF_LA_CUDA` | `GPU=nvidia` | `alf_la` backend cuBLAS/cuSOLVER on device-resident arrays. Library handles use the stream of the OpenMP synchronous regions. cuSOLVER deterministic mode on and cuBLAS atomics off are fixed (no switch). | production | `alf_la_mod.F90:1195-1217` |
| `ALF_LA_ONEMKL` | `GPU=intel` | `alf_la` backend oneMKL (OpenMP dispatch for BLAS, SYCL shim for QR/LU). Also enables the async block path code in `gpu_lowrank_mod.F90`. | production | `alf_la_mod.F90:105`, `gpu_lowrank_mod.F90:128` |
| `ALF_LA_HIP` | `GPU=amd` | `alf_la` backend hipBLAS/hipSOLVER via C shim, plus own HIP kernels (`alf_hipk_mod`). | production | `alf_la_hip_mod.F90:8` |
| `ALF_LA_OM` | derived, do not set | Defined inside `alf_la_mod.F90` when `ALF_LA_ONEMKL` or `ALF_LA_HIP` is set. `LA_GUARD` heads the host branch of 16 public routines; under `ALF_LA_OM` that branch is preprocessed away because each of them calls its `om_*` twin (oneMKL: `alf_la_om_mod`; HIP: `alf_la_hip_mod`), and a routine without device code would stop the run instead of calling host BLAS on device data. No public `la_*` routine lacks a device implementation today, so `LA_GUARD` never fires in this tree (`NOT_IMPLEMENTED.md` 5.1). | production | `alf_la_mod.F90:105-120` |
| `ALF_LA_HOST` | implicit when no backend macro is set | Plain BLAS/LAPACK. Only tested for a clash with `ALF_LA_CUDA` (error). | production | `alf_la_mod.F90:29,111` |
| `ALF_LA_SYNC_ALWAYS` | off | CUDA only. Brackets every library call with `la_sync`. Safe fallback if stream ordering is in doubt; slow. Not used by the OM backends. | diagnostic | `alf_la_mod.F90:1323,1330` |
| `ALF_LA_DEBUG` | off | CUDA only. Prints the library stream handle once in `la_init` (`0x2` expected). | diagnostic | `alf_la_mod.F90:1208` |
| `ALF_LA_TEST_BADSTREAM` | off | CUDA only. `la_init` uses a `cudaStreamCreate` stream, not ordered with OpenMP regions. Negative control for `t_order`; it must fail. | test-only | `alf_la_mod.F90:1201,1238` |
| `ALF_TEAM_TIMERS` | off | Wall time per sweep phase (sweep, wrapgrup, wrapgrdo, wrapur, wrapul, cgr, obser) and block timers of the lowrank path; one `TEAM_TIMERS` line per bin on stdout. Small overhead from `system_clock`. | diagnostic | `Prog/main.F90:159,1399`, `gpu_lowrank_mod.F90:79` |
| `ALF_TIMERS` | off | Exclusive stack timers of the sweep (update, flush, vwrap, hop, stab, meas, ...). Appends to `timers.txt` per bin. Without the macro the module is empty. | diagnostic | `Prog/Timers_mod.F90:40,79` |
| `ALF_TRACE_UPDATES` | off | Appends one line per Metropolis decision to `trace_updates.txt` (for trajectory comparison of update schemes). Large file, slow. | diagnostic | `Prog/Trace_updates_mod.F90:63`, `upgrade_mod.F90:116` |
| `ALF_LOWRANK_DEBUG` | off | Shadow Green function for delayed/submatrix updates; checks every accept and flush; terminates on a violation; writes `lowrank_debug_summary.txt`. Slow. | diagnostic | `Prog/Lowrank_update_mod.F90:117`, `Wrapgr_mod.F90:68` |
| `QR_PREPIVOT` | off; set with `ALF_FLAGS_EXT=-DQR_PREPIVOT` | Pre-pivoted QR in `cgr2_2_mod` (the CPU variant of the device algorithm). Needed as host oracle for `ALF_GPU_SELFTEST=3` and `=4` with the default stabilization; those self-tests refuse to run without it. The meaning of the value (`==2`, `==3`) is in `cgr2_2_mod.F90:405-472`, not described here. | test-only | `gpu_tau_mod.F90:2118`, `gpu_stab_mod.F90:1709` |

`STABLOG`, `STAB1/2/3`, `MPI`, `TEMPERING`, `PARALLEL_PARAMS` are upstream macros set through the `configure.sh` mode arguments. The GPU path
supports default and `STABLOG` stabilization; `STAB1/2/3` makes it inactive (`gpu_state_mod.F90:495-540`).

## 3. Runtime switches common to all backends

Rows with `ALF_GPU_*` are subject to the guard of section 0: they need `ALF_GPU_DEBUG=1` unless stated. Only `ALF_GPU_DISABLE`, `ALF_GPU_DEBUG` and
`ALF_UNSTABLE_THRESHOLD` are allowed in production. "Hook" switches (self-tests, timings) run once after setup, before the first sweep.

| Name | Values / default | Effect | Category | Where |
|---|---|---|---|---|
| `ALF_GPU_DISABLE` | `1`; default unset (GPU path on if supported) | The GPU build runs the host path (reason `disabled by ALF_GPU_DISABLE=1` in `info`). Allowed in production without DEBUG. | production | `gpu_state_mod.F90:513` |
| `ALF_GPU_DEBUG` | any value other than `0` (`1` by convention) declares a debug run; `0`, empty or unset: production | Only effect: lifts the guard for the other `ALF_GPU_*` variables. | test-only | `gpu_state_mod.F90:409,415` |
| `ALF_UNSTABLE_THRESHOLD` | positive real; default 10 (upstream) | Abort threshold of the G precision check (`XMAX > threshold` aborts). Smaller value = stricter. Non-positive or unparsable value: exit 2. Not GPU-specific, not guarded. Production on SuperMUC uses 1e-4 (DECISIONS.md, decision 46). | production | `Prog/control_mod.F90:98-112` |
| `ALF_SHM_CHUNK_SIZE_GB` | positive integer; default unset = normal allocation | Upstream: chunk size of MPI shared-memory blocks. Not a GPU switch. | production | `Prog/main.F90:226,305` |
| `ALF_GPU_SELFTEST` | any non-empty value except `0` runs the mirror self-test and the run continues, **except `3` and `4`**, which skip it (`gpu_state_mod.F90:686-695`: the code excludes exactly these two); `1` is the usual value. `2` runs the mirror self-test **and** the wrap hook (`gpu_wrap_mod.F90:1254`); `3` stabilization hook; `4` tau pass hook; unset or `0`: none. Any value other than `0` also makes the run a debug run | Compares device with host on fixed fields. Values `fault`, `fault_c`, `fault_v`, `nan`, `nan_v` plant a 1e-9 or NaN defect in the host copy (negative controls, must fail with a nonzero exit). If the GPU path is inactive the run stops. `3` and `4` need `QR_PREPIVOT` for the default stabilization. | test-only | `gpu_state_mod.F90:642,686-700` |
| `ALF_GPU_UPDTEST` | `1` | One-slice test of the device update at the first `WRAPGRUP` call; state restored afterwards (Control acceptance counters are not; combine with `ALF_GPU_EXIT=1`). Declares a debug run. | test-only | `gpu_update_mod.F90:167`, `Wrapgr_mod.F90:547` |
| `ALF_GPU_EXIT` | `1` | Ends the run with exit 0 after a hook (`SELFTEST`=2/3/4, `*TIME`, `TIMING`, `UPDTEST`) has run. No effect alone. | test-only | `gpu_wrap_mod.F90:1262`, `gpu_stab_mod.F90:1589`, `gpu_tau_mod.F90:2014`, `Wrapgr_mod.F90:663` |
| `ALF_GPU_TIMING` | `1` | Hook: timing of the device wraps vs host. | diagnostic | `gpu_wrap_mod.F90:1258` |
| `ALF_GPU_STABTIME` | `1` | Hook: time per stabilization step, device vs host. | diagnostic | `gpu_stab_mod.F90:1576` |
| `ALF_GPU_TAUTIME` | `1` | Hook: one tau pass without observables, device vs host. Needs the device tau pass (`Ltau = 1`), else the run stops. | diagnostic | `gpu_tau_mod.F90:1999` |
| `ALF_GPU_TAUTEST_V` | `1` | Verbose per-slice differences of the `SELFTEST=4` self-test. | diagnostic | `gpu_tau_mod.F90:2138` |
| `ALF_GPU_NOUPDATE` | `1`; default unset (device update path) | Keeps the host `WRAPGRUP/DO` in a GPU build (A/B tests). Also changes which copies the tau/stabilization code makes. | test-only | `gpu_update_mod.F90:105` |
| `ALF_GPU_NOOBS` | `1`; default unset (device observables where supported) | Keeps the host `Obser` (A/B tests). Read with `get_environment_variable`, still covered by the guard. | test-only | `gpu_obs_mod.F90:93` |
| `ALF_GPU_OBS_SKIP` | comma list of kinds from 1..4 (1 Green_eq, 2 SpinZ_eq for N_FL=1, 3 Den_eq, 4 SpinMz); default empty | Those equal-time observables are not measured on the device and **stay zero in the output**. Never use for results. | test-only | `gpu_obs_mod.F90:95-100,349` |
| `ALF_GPU_TRACE` | file name; default unset | Per `WRAPGRUP/DO` call writes the flipped sites of the slice and the phase to that file (first-divergence tests). | diagnostic | `gpu_update_mod.F90:184` |
| `ALF_GPU_TL`, `ALF_GPU_TR` | integer >= 1; defaults 1 and 2 | Columns (left) / rows (right) per thread of the per-group wrap kernels. Unparsable or < 1: ignored. Tuned on the A30. Whether the result is bitwise independent of the value: not determined. | tuning | `gpu_wrap_mod.F90:98,1250-1251` |
| `ALF_GPU_VB` | integer >= 1; default 8 | Columns/rows per thread of the fused Op_V slice kernels. Bitwise independence: not determined. | tuning | `gpu_wrap_mod.F90:99,1252` |
| `ALF_GPU_NS` | `0` switches off; default on | Newton-Schulz step on Q in the device QR (stabilization). Off lowers the orthogonality of Q (measured Precision Green 9.7e-11 vs 5.4e-11 mean, ticket 206). Debug only. | test-only | `gpu_stab_mod.F90:158,520` |
| `ALF_GPU_TAU_NS` | `0` switches off; default on | Same for the explicit Q in the tau pass (`CGR2_2`); off = `ZUNMQR` applied directly to X, as upstream. Needs the device tau pass. | test-only | `gpu_tau_mod.F90:254` |
| `ALF_GPU_CGR_LU` | `1` or `2`; `0`, unset or any other value: off | Determinant phase in the device `CGR` (`gpu_stab_mod.F90:176-177,887-898`). Default: the tracked det(U) of the two U factors where that is enabled (oneMKL only), else the LU of the RHS (`la_getrf`). `1`: always the LU (the path before the LU-free phase). `2`: compute both, use the tracked phase, print `GPU CGR phase check (ALF_GPU_CGR_LU=2): CGR calls compared ..., largest |phase_tau - phase_LU|` at the end. Needs `ALF_GPU_DEBUG=1` (guard). **Only has an effect with `GPU=intel`**: with CUDA and HIP `decompose_slot` sets `dw_ok = .false.` (`#if !defined(ALF_LA_ONEMKL)`), so the LU is always used and `2` compares nothing. | diagnostic | `gpu_stab_mod.F90:176,215,471,887` |
| `ALF_GPU_STABCHECK` | `1`: every CGR repeated on the host; `2`: every WRAPUR/WRAPUL block; default off | Host check of the device states, deviations summarised in `info`. Slow. | diagnostic | `gpu_stab_mod.F90:1542` |
| `ALF_GPU_STABMIX` | `1`: use the host CGR result; `2`: the host WRAPUR/UL result; default off | Isolation experiment (the host result replaces the device result). Implies a check level. | test-only | `gpu_stab_mod.F90:1544-1550` |
| `ALF_GPU_FAULT` | `hop`, `hopc`, `vexp`, `nan`, `order` (with `SELFTEST=2`); `upd`, `nsig`, `hostonly` (update path); other value: error | Plants a defect in the device copy only (negative control). `hostonly`: the "device" run of the test is a second host run. | test-only | `gpu_state_mod.F90:423`, `gpu_wrap_mod.F90:1253,1271`, `gpu_update_mod.F90:143`, `Wrapgr_mod.F90:578` |
| `ALF_GPU_FAULT_STAB` | `perm`, `phase`, `nan`, `eps:<x>` with 0 < x < 1e300 | Negative controls for the stabilization. `eps:<x>` multiplies the largest element of the device G by (1+x), in the `SELFTEST=3` comparison, or (decision 26) in every stabilization of a declared control sweep. | test-only | `gpu_state_mod.F90:431-450`, `gpu_stab_mod.F90:1557-1570` |
| `ALF_GPU_FAULT_TAU` | `perm`, `prop`, `prec`, `blocks`, `eps:<x>` | Negative controls for the tau pass. `eps:<x>` perturbs the four device outputs of `CGR2_2` (self-test, or every pass of a control sweep). Needs the device tau pass. | test-only | `gpu_state_mod.F90:453-472`, `gpu_tau_mod.F90:1955,1966` |
| `ALF_GPU_FAULT_TAU_AT` | `<nst>:<nfe>` (`nst` 0 = every stabilization point; `nfe` >= 1 flavour) | Restricts the `eps` fault of the sweep to one stabilization point and flavour. Out-of-range value stops the run. | test-only | `gpu_tau_mod.F90:1974-1991` |
| `ALF_GPU_TAUSWAP` | integer mask 0..4095 | Replaces steps of the device `CGR2_2` by host twins in the production calls (bits listed at `gpu_tau_mod.F90:170-174`). Isolation experiment. | test-only | `gpu_tau_mod.F90:830` |
| `ALF_GPU_TAUDIAG` | `1` (all variants), `2` (no step twins), `0`/unset off | Per stabilization point and flavour one `TAUDIAG` line: device precision numbers against a lockstep host tau pass, 1-ulp perturbed chains. Very slow. | diagnostic | `gpu_tau_mod.F90:835-839` |
| `ALF_GPU_TAUDIAG_FROM` | integer >= 0; default 0 | First pass for the perturbation chains (and `QREF` dumps). | diagnostic | `gpu_tau_mod.F90:840` |
| `ALF_GPU_TAUDIAG_SEED` | integer >= 0; default derived from the run RNG state | Seed of the 1-ulp perturbations. | diagnostic | `gpu_tau_mod.F90:846` |
| `ALF_GPU_TAUDIAG_VAR` | comma list (max 10) of masks 0..4095 | Extra chains of `dev_cgr22` with those `TAUSWAP` masks. | diagnostic | `gpu_tau_mod.F90:855` |
| `ALF_GPU_TAUDIAG_AMP` | comma list (max 10) of reals in (0, 1e10) | Calibration chains with perturbation amplitude a * 4e-16. | diagnostic | `gpu_tau_mod.F90:869` |
| `ALF_GPU_TAUDIAG_NUL` | comma list (max 12) of integer codes (`flags + 100*nb + 100000*stream`) | Host `CGR2_2` chains with rounding-order "null variants"; code table at `gpu_tau_mod.F90:208-212`. | diagnostic | `gpu_tau_mod.F90:883` |
| `ALF_GPU_TAUDIAG_QREF` | stride >= 0; default off | Dumps `CGR2_2` inputs/outputs to `qref/q<point>.bin` in the run directory (directory must exist) for an offline quad-precision reference. | diagnostic | `gpu_tau_mod.F90:691` |
| `ALF_GPU_TAUDIAG_QVAR`, `_QHOST` | comma list of masks 0..4095; `0` or `1` | With `QREF`: also dump `dev_cgr22` outputs for those masks / the host `CGR2_2` outputs. | diagnostic | `gpu_tau_mod.F90:697,710` |

If the device tau pass is off (for example `Ltau /= 1`), any of `SELFTEST=4`, `TAUTIME`, `FAULT_TAU[_AT]`, `TAUSWAP`, `TAU_NS`, `TAUTEST_V`, `TAUDIAG*`
stops the run, so a set switch is never silently ignored (`tau_off_refuse`, `gpu_tau_mod.F90:286-300`).

## 4. Intel (`GPU=intel`, oneMKL + SYCL shim)

These variables are read by the C++ shim and libraries, not by Fortran, so the guard of section 0 does **not** cover them.
Only the first character of the value is tested unless stated. Status of the production decision: `docs/gpu/sessions/supermuc_large_prompt/ISSUE_INTEL_VS_A30.md`, section 9.

| Name | Values / default | Effect | Category | Where |
|---|---|---|---|---|
| `ALFSH_ASYNCBLK` | `1` (exactly); default unset = off | Submatrix block sequence runs as SYCL kernels plus oneMKL GEMM on the OpenMP interop queue, without a wait per op (library `libalfsb.so`). Prints `ALFSH_ASYNCBLK=1, submatrix blocks run asynchronously`. **Status: bitwise identical to off (measured, V2), measured gain, not yet cleared for production by the owner** (release step "V2 next" in the Intel doc, after a same-seed check of the production binary; `SPEEDUPS.md`, `NOT_IMPLEMENTED.md` 6). The measured gain differs with the test: +15 to 20 % per card in the one-chain comparison on one card (L=24: 7.05 to 6.55 s per sweep, `ISSUE_INTEL_VS_A30.md` section 2, V1 vs V2), but only +3 % per node with 32 ranks on 8 tiles (2.57 to 2.65 sweeps/s per node; `sessions/supermuc_large_prompt/phase2/team/p3/STATUS.md`, phase 3b step 3). With several ranks per tile the device is already shared, so less waiting is left to hide. | tuning | `Prog/gpu_lowrank_mod.F90:1334-1346` |
| `ALFSH_DEVQR` | `1`: own device ZGEQRF (panel + ZLARFT, GEMM update); `2`: T-free panel + apply kernel; other/unset: oneMKL (default) | Loads `libalfdevqr.so` by `dlopen`; stops the run if the library or a symbol is missing. Falls back to oneMKL for large M. Rounding-level differences from oneMKL; the Intel doc says "not yet" for production (waits for numerics review; `=2` also for the L6T rule and tau-slot review). | tuning | `alf_la_sycl_shim.cpp:119,126-148` |
| `ALFSH_DEVUNMQR` | `1`; default off | Device ZUNMQR (apply kernel, side L/R, trans N/C), same library and status as `DEVQR=2`. Works without `DEVQR`. | tuning | `alf_la_sycl_shim.cpp:121` |
| `ALFSH_QRNB` | integer > 0; default 32 (16 if `DEVQR=2` and unset) | Block size of the device QR. | tuning | `alf_la_sycl_shim.cpp:123-125` |
| `ALFSH_DEVQR_LIB` | path; default `<exe dir>/../Libraries/Modules/libalfdevqr.so`, else `libalfdevqr.so` on the loader path | Which device-QR library to load (used to plant defects in tests). | tuning | `alf_la_sycl_shim.cpp:127-134` |
| `ALFSH_TIMING` | `1`; default off | Counts calls and wall time per shim routine; prints a table at `alfsh_finalize` (`ALFSH_TIMING <routine> calls n total s`). No effect on results. With the async path it includes queued work. | diagnostic | `alf_la_sycl_shim.cpp:71-80` |
| `ALFSH_LAUNCHBENCH` | any value (presence); default unset | At the first device-QR call prints the launch cost of trivial kernels/GEMMs and panel timing (`ALFSH_LAUNCHBENCH`, `ALFSH_DEVINFO`, `ALFSH_PROBE` lines; those three are output labels, not variables). Only with `DEVQR`. | diagnostic | `alf_la_sycl_shim.cpp:228`, `alf_la_sycl_devqr.cpp:70` |
| `ALFDEVQR_G` | 1, 2 or 3; default by panel height (1 for M <= 256, else 2) | Forces the number of sub-groups per panel column in the T-free panel kernel (`DEVQR=2`). Experiments. | tuning | `alf_la_sycl_qr2.hpp:212` |
| `ALFDEVQR_ANV` | 2 or 4; default 2 | Columns per sub-group of the apply kernel for 289 <= nq <= 1152 (`4` selects the other instantiation). | tuning | `alf_la_sycl_qr2.hpp:384` |
| `ALFDEVQR_NBS`, `ALFDEVQR_NBSM` | integers; defaults 32 and 256 | Panel width used for short panels (height <= `NBSM`) when `nb == 16` (`DEVQR=2`). | tuning | `alf_la_sycl_qr2.hpp:393-394,423` |
| `ALFDEVQR_PROF` | any value (presence) | Per-stage host timing with a wait after every launch; prints every 50 calls to stderr (`ALFDEVQR_PROF`, `..._PROF2`). Changes timing, not results. | diagnostic | `alf_la_sycl_devqr.hpp:319`, `alf_la_sycl_qr2.hpp:416` |
| `ALFDEVQR_SYNC` | any value (presence) | Host wait after every launch of the `DEVQR=1` kernels. | diagnostic | `alf_la_sycl_devqr.hpp:323` |
| `ALF_LA_TEST_BREAK_ORDER` | `1`; test hooks of the oneMKL backend | GEMM goes to a second, unordered queue; nobody waits. Must make the order test fail. | test-only | `alf_la_mod.F90:396,564` |
| `ALF_LA_TEST_POISON_INFO` | `1` | Solver info slots are never read back; failures go unnoticed. | test-only | `alf_la_mod.F90:397,509` |
| `ALF_LA_TEST_WRONG_LD` | `1` | `la_gemm` uses a valid but wrong leading dimension (`lda = m`). | test-only | `alf_la_mod.F90:398,563` |

Each armed hook prints `alf_la test hook armed: <name>`. They are also refused in production by the guard (prefix `ALF_LA_TEST_`).
`ALFSH_BATCH`, mentioned in old team scripts, no longer exists in the source.

## 5. AMD (`GPU=amd`, hipBLAS/hipSOLVER + own HIP kernels)

| Name | Values / default | Effect | Category | Where |
|---|---|---|---|---|
| `ALF_HIP_QR` | first character `0` = off; default on | Own blocked Householder QR and Q-apply (`alf_hip_qr.hip.cpp`) for geqrf/unmqr with m <= 2304; `0` (or larger m) uses hipSOLVER. About 2x faster geqrf and 5x unmqr at n=576. Whether the result is bitwise equal to hipSOLVER: not determined. | tuning | `alf_la_hip_mod.F90:38,251,518` |
| `ALF_HIP_ASYNC` | first character `0` = off; default on **only if** the patched libgomp GCN plugin registered its presync hook, else off | Deferred synchronisation: HIP calls are not followed by `hipDeviceSynchronize`; the patched plugin calls `alf_hip_presync()` before every OpenMP kernel launch and device copy. Output bit-identical to the synchronous mode (README of the AMD platform). Prints `HIP deferred synchronisation: ON/OFF (...)` once. Use `0` to debug or for A/B timing. | production | `alf_hip_kernels.hip.cpp:25-46` |
| `ALF_LA_TEST_*` | not active | `alf_la_hip_mod.F90` declares the three hook flags and a `hook_armed` function, but nothing calls it; the hooks cannot be armed on HIP. The guard still refuses the names. | test-only | `alf_la_hip_mod.F90:25,218,345` |

All section 3 switches also apply on AMD. Build options: section 1. Vendor variable `GCN_NUM_THREADS`: section 7.

## 6. NVIDIA (`GPU=nvidia`, cuBLAS/cuSOLVER)

No runtime variable is read only by the CUDA backend. All section 3 switches apply. Build/macro items that are CUDA only:
`ALF_GPU_CC`, `ALF_HELMA_NVHPC` (section 1); `ALF_LA_SYNC_ALWAYS`, `ALF_LA_DEBUG`, `ALF_LA_TEST_BADSTREAM` (section 2).

| Setting | Value | Effect | Category | Where |
|---|---|---|---|---|
| cuSOLVER deterministic mode | fixed on | Bitwise reproducible solver results. No switch. Needs CUDA libraries with `cusolverDnSetDeterministicMode` (CUDA 11.8 default libs lack it, hence `cuda12.6` on Helma). | production | `alf_la_mod.F90:1215`, `configure.sh:508` |
| cuBLAS atomics | fixed off | Reproducible library results. No switch. | production | `alf_la_mod.F90:1216` |
| One GPU per process | by device selection | The OpenACC "current device" used for the stream query must equal the OpenMP default device. Bind ranks with `CUDA_VISIBLE_DEVICES` (see section 7). In the code, rank `r` uses OpenMP device `r mod number of devices`, with `r` the **global** `MPI_COMM_WORLD` rank, not the node-local one (`MPI_GPU.md`, "Rank-to-device binding"). | production | `alf_la_mod.F90:44-47`, `gpu_state_mod.F90:552` |

## 7. Vendor and runtime variables used in job scripts

These are set by the user or the job script, not read by ALF. None starts with `ALF_GPU_`, so the guard ignores them.

| Name | Typical value | Effect | Where used |
|---|---|---|---|
| `OMP_TARGET_OFFLOAD` | `MANDATORY` | Fail loudly if no device is found instead of running on the host. Used in all production and test scripts. | `docs/gpu/MPI_GPU.md:12`, `scripts/supermuc_ng2/throughput.sbatch:14` |
| `OMP_NUM_THREADS` | `1` (Helma), `14` (SuperMUC, one rank per tile) | Host threads per rank. | same |
| `LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE` | `1048576` (1 MB) | Intel. Staging buffer of the Level Zero plugin; the 32 KB D->H copy per block took 380 us with the default and 22 us with 1 MB. Largest single lever (L=24: 10.6 -> 7.1 s/sweep); bitwise identical. Recommended "V1" for production. | `ISSUE_INTEL_VS_A30.md` sec. 2, 9 |
| `ZE_FLAT_DEVICE_HIERARCHY` | `FLAT` | Intel. Each tile is one device (8 per node); rank `r` uses tile `r mod 8`. Flat is the default; set it explicitly. | `scripts/supermuc_ng2/env.sh:3`, `platforms/intel_pvc_supermuc_ng2.md:9,26` |
| `ZEX_NUMBER_OF_CCS` | `0:4,...,7:4` | Intel. 4 compute command streamers per tile; with 4 ranks per tile it is the throughput optimum (32 ranks per node, `sessions/supermuc_large_prompt/COMPARISON.md`, later sessions with the other switches in place). No gain for a single process, and no gain in the first session (`platforms/intel_pvc_supermuc_ng2.md`, 8 ranks, without the later switches: `ZEX_NUMBER_OF_CCS=0:4` listed there under "no gain"). The two statements are about different conditions (rank count per tile, switches), not a contradiction. | `sessions/supermuc_large_prompt/COMPARISON.md:21`, `SYNC.md:16` |
| `ZE_AFFINITY_MASK` | `0` | Intel. Restricts a run to one tile; used in single-tile timing and test scripts. | `sessions/.../team/*.sh` |
| `ONEAPI_DEVICE_SELECTOR` | `level_zero:0` | Intel test machine (Arc B580): select the GPU. Do not run plain `sycl-ls` there (CLAUDE.md). | `docs/gpu/MACHINES.md`, `CLAUDE.md` |
| `I_MPI_PIN_DOMAIN`, `KMP_AFFINITY` | `auto`, `granularity=fine,compact,1,0` | Intel MPI/OpenMP pinning in the throughput job. | `scripts/supermuc_ng2/throughput.sbatch:18` |
| `LIBOMPTARGET_PLUGIN_PROFILE`, `LIBOMPTARGET_DEBUG` | `T`, `1` | Intel offload runtime profiling and call counts (diagnostic; slow). | `sessions/.../team/p7/NOTES.md` |
| `GCN_NUM_THREADS` | `4` (default 16) | AMD. Threads of the GCN OpenMP runtime; parallel-region startup scales with it; 4 was best at L=16. | `platforms/amd_mi210/README.md:14` |
| `GCN_PROF` | `1` | AMD. Per-kernel wall time at exit; only with the patched libgomp GCN plugin. | `platforms/amd_mi210/README.md:10` |
| `CUDA_VISIBLE_DEVICES` | GPU index per rank | NVIDIA. Bind each rank to one GPU (`rank / ranks_per_GPU`); several ranks share one GPU via MPS. | `scripts/helma/rank_wrap.sh` |
| `CUDA_MPS_PIPE_DIRECTORY`, `CUDA_MPS_LOG_DIRECTORY` | per-GPU directories | NVIDIA. Select the MPS daemon of the GPU; one daemon per GPU. | `scripts/helma/run_pp.sh:12-33`, `rank_wrap.sh` |
| `CUDA_MPS_ACTIVE_THREAD_PERCENTAGE` | `25` (production) | NVIDIA. SM share per MPS client; 25 % and 50 % equal within 2 % at 8-16 ranks. MPS is mandatory for several ranks per GPU (2 ranks without it: 2.6-6.5x slower per chain). | `platforms/helma_h100.md:58,112-115` |

The restart of a production chain is by `confin` files, not by a variable (`docs/gpu/platforms/helma_h100.md`).
