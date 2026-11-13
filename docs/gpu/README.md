# ALF GPU port

> **Branch `gpu-port`: the general GPU code for all ALF users.** Upstream ALF (`master` 521ddfca)
> plus opt-in GPU backends; the Hubbard model and the measurement window are upstream's, and a
> `GPU=none` build is bitwise identical to master. The space-time point-defect model lives on
> `gpu-port-unified`. Where this page or the pages it links mention the defect, `ham_h0` on
> time slice 1, `SpinZ_localtau` or `LOBS_ST = 0`, they describe that branch. See
> [BRANCHES.md](BRANCHES.md).

The GPU port runs the ALF auxiliary-field QMC sweep on a GPU. The Green
function, the delayed and submatrix updates, the wraps, the stabilisation and the equal-time
measurements stay on the device for the whole sweep. One source tree builds for four targets:

| `GPU=` | Hardware tested | Compiler | Libraries | Backend macro |
|---|---|---|---|---|
| `none` (default) | any CPU | any ALF compiler | BLAS/LAPACK | none: the CPU code is unchanged |
| `nvidia` | A30, H100 | nvfortran (nvhpc 25.7 / 26.9) | cuBLAS, cuSOLVER | `ALF_LA_CUDA` |
| `intel` | Max 1550 (SuperMUC-NG Phase 2) | ifx (oneAPI 2025.x) | oneMKL, SYCL shim | `ALF_LA_ONEMKL` |
| `amd` | MI210 | gfortran 14 with amdgcn offload | hipBLAS, hipSOLVER, HIP kernels | `ALF_LA_HIP` |

The port was developed for the Hubbard model with a space-time point defect (time-dependent
hopping, `N_FL = 2`), with MPI and parallel tempering (`PARALLEL_PARAMS`). What is and is not
supported is listed in [FEATURES.md](FEATURES.md) and [NOT_IMPLEMENTED.md](NOT_IMPLEMENTED.md).

## Build

```bash
source configure.sh PGI MPI GPU=nvidia HDF5        # NVIDIA (HELMAGPU on Helma)
source configure.sh IntelLLVM MPI GPU=intel HDF5   # Intel
source configure.sh GNU MPI GPU=amd HDF5           # AMD
make lib program
```

Details per machine, module versions and known compiler problems: [BUILD.md](BUILD.md),
[MACHINES.md](MACHINES.md) and the platform pages below.

## How correctness is checked

A GPU build is compared with a CPU build made by the same compiler (`GPU=none`), same
parameters and seeds. The trajectories must agree to rounding (about 1e-13 with `Mz = .F.`,
about 1e-10 with `Mz = .T.`). Method, tolerances and the scripts are in
[VERIFICATION.md](VERIFICATION.md); the driver is `scripts/defect/run_cmp.py`.
The CI (`.github/workflows/tests.yml`) builds and tests the CPU path only.

## Speed (seconds per sweep, Beta = L, submatrix update, Nwrap 10)

| Platform, layout | L=24 | L=32 | Source |
|---|---|---|---|
| A30, one chain | 5.25 | 18.2 | [TUNING_DEFECT.md](TUNING_DEFECT.md) |
| A30, 4 ranks with MPS | 10.5 per chain | | [MPI_GPU.md](MPI_GPU.md) |
| H100, one chain | 3.51 | 9.82 | [platforms/helma_h100.md](platforms/helma_h100.md) |
| H100, 8 ranks with MPS | 6.80 per chain | 34.5 per chain | same |
| Max 1550, one chain per tile | 4.95 | 16.4 | [ISSUE_INTEL_VS_A30.md](sessions/supermuc_large_prompt/ISSUE_INTEL_VS_A30.md) |
| Max 1550, 4 chains per tile | 9.3 per chain | | same |
| MI210, one process | 4.41 | 17.5 | [platforms/amd_mi210/README.md](platforms/amd_mi210/README.md) |

The Intel single-chain numbers need the opt-in switches listed in
[RUNTIME_SWITCHES.md](RUNTIME_SWITCHES.md). Open speedup candidates: [SPEEDUPS.md](SPEEDUPS.md).

## Where things are

| Topic | Document |
|---|---|
| Design: programming model, `alf_la` interface, device data | [ADR-001](ADR-001-programming-model.md), [ADR-002](ADR-002-alf-la-and-data-model.md) |
| Delayed and submatrix update on the device | [DELAYED_SUBMATRIX.md](DELAYED_SUBMATRIX.md) |
| Time-dependent hopping (point defect) | [TIMEDEP_HOPPING.md](TIMEDEP_HOPPING.md) |
| MPI and parallel tempering with the GPU path | [MPI_GPU.md](MPI_GPU.md), [PP_GPU.md](PP_GPU.md) |
| Tuning (`N_delay`, `Nwrap`) for the defect runs | [TUNING_DEFECT.md](TUNING_DEFECT.md) |
| Build macros and environment switches | [RUNTIME_SWITCHES.md](RUNTIME_SWITCHES.md) |
| Benchmarks, CPU baseline and profile | [BENCHMARKS.md](BENCHMARKS.md), [BASELINE.md](BASELINE.md), [PROFILE.md](PROFILE.md) |
| What makes the port hard; plan; owner decisions | [HARD.md](HARD.md), [PLAN.md](PLAN.md), [DECISIONS.md](DECISIONS.md) |
| Claims in the literature or docs that measurement contradicted | [DISCREPANCIES.md](DISCREPANCIES.md) |
| How to switch the GPU path off | [ROLLBACK.md](ROLLBACK.md) |
| How this branch was put together | [CONSOLIDATION.md](CONSOLIDATION.md) |

Platforms:

| Platform | Page | Scripts |
|---|---|---|
| NVIDIA A30 (`nividia`) | [TUNING_DEFECT.md](TUNING_DEFECT.md), [MPI_GPU.md](MPI_GPU.md) | `scripts/defect/`, `scripts/pp/`, `scripts/nividia/` |
| NVIDIA H100 (Helma, NHR@FAU) | [platforms/helma_h100.md](platforms/helma_h100.md) | `scripts/helma/` |
| Intel Max 1550 (SuperMUC-NG Phase 2) | [platforms/intel_pvc_supermuc_ng2.md](platforms/intel_pvc_supermuc_ng2.md) (first build and layout), [sessions/supermuc_large_prompt/](sessions/supermuc_large_prompt/) (device QR, async blocks, production) | `scripts/supermuc_ng2/`, `sessions/supermuc_large_prompt/{jobs,prod,phase2}` |
| AMD MI210 | [platforms/amd_mi210/README.md](platforms/amd_mi210/README.md) | `platforms/amd_mi210/` |

Working material of the development sessions (job scripts, logs, measurement notes, WBS files)
is kept under [sessions/](sessions/), [tickets/](tickets/) and [reports/](reports/). It records
how each result was obtained and is not needed to build or run the code.

## Source layout

| Files | Content |
|---|---|
| `Libraries/Modules/alf_la_mod.F90` | `la_*` interface: one API over host LAPACK, cuBLAS/cuSOLVER, oneMKL and HIP |
| `Libraries/Modules/alf_la_sycl_*.{cpp,hpp}` | Intel: SYCL shim, async block path (`libalfsb.so`), device Householder QR (`libalfdevqr.so`) |
| `Libraries/Modules/alf_la_hip_mod.F90`, `alf_hipk_mod.F90`, `alf_hip_*.cpp`, `alf_la_hip_shim.cpp` | AMD: HIP backend, kernels and device QR |
| `Prog/gpu_state_mod.F90` | device state, environment switches, activation rules |
| `Prog/gpu_wrap_mod.F90`, `gpu_stab_mod.F90` | wraps and stabilisation (QR, CGR) on the device |
| `Prog/gpu_update_mod.F90`, `gpu_lowrank_mod.F90` | local updates: delayed and submatrix scheme, block mode |
| `Prog/gpu_obs_mod.F90`, `gpu_tau_mod.F90` | equal-time and time-displaced measurements |
| `testsuite/alf_la.tests/`, `testsuite/gpu_reference/` | unit tests of `alf_la`, statistical and acceptance harness |
