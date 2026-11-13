# Building ALF for GPU targets

`configure.sh` takes an optional argument `GPU=none|nvidia|intel|amd` (case-insensitive, any position):

```sh
. ./configure.sh PGI noMPI HDF5 GPU=nvidia NO-INTERACTIVE      # nvfortran + OpenMP offload + cuBLAS/cuSOLVER
. ./configure.sh IntelLLVM noMPI GPU=intel NO-INTERACTIVE      # ifx + OpenMP offload (spir64 JIT) + oneMKL
. ./configure.sh GNU noMPI GPU=amd NO-INTERACTIVE              # gfortran >= 14 amdgcn offload + hipBLAS/hipSOLVER + HIP kernels
make -j8 all
```

`GPU=none` (or no `GPU=`) leaves every exported `ALF_*` variable identical to the upstream script
(checked for GNU/IntelLLVM/PGI x noMPI/MPI x default/LOG x HDF5 x Devel with `docs/gpu/scripts/nividia/env_diff.sh`).

| `GPU=` | MACHINE | compile flags (modules, prog, ana) | added to `ALF_LIB` | defines | exported |
|---|---|---|---|---|---|
| none | any | - | - | - | - |
| nvidia | `PGI` or `HELMAGPU` (`HELMAGPU` loads `gpu-env/2025`, **nvhpc 25.7** + `openmpi/5.0.8-nvhpc25.7-cuda`, **CUDA 12.9**: `ALF_GPU_CC` default `cc90,cuda12.9`; nvhpc 24.11 is broken, see below) | `-mp` replaced by `-mp=gpu -gpu=${ALF_GPU_CC:-cc80}` | `-mp=gpu -gpu=${ALF_GPU_CC:-cc80} -cudalib=cublas,cusolver` | `-DALF_GPU -DALF_LA_CUDA` | `ALF_GPU=nvidia` |
| intel | `IntelLLVM`/`IntelX` only | `-fiopenmp -fopenmp-targets=spir64 -I${MKLROOT}/include`, `-static` removed | `-fiopenmp -fopenmp-targets=spir64 -qmkl -lmkl_sycl_blas -lmkl_sycl_lapack -lsycl -lstdc++ -ldl -L<ALF>/Libraries/Modules -lalfsb -Wl,-rpath,<ALF>/Libraries/Modules` (`configure.sh:641`; `libalfsb.so` is built by `Libraries/Modules/Makefile`, target `lib`, part of `make all`) | `-DALF_GPU -DALF_LA_ONEMKL` | `ALF_GPU=intel` |
| amd | `GNU` only (gfortran >= 14 with GCC amdgcn offload, ROCm with hipBLAS/hipSOLVER/hipcc) | `-foffload=amdgcn-amdhsa -foffload-options=amdgcn-amdhsa=-march=${ALF_AMD_ARCH:-gfx90a}`, `-static` and `-ffast-math` removed (runs went unstable with it) | `-foffload=amdgcn-amdhsa -foffload-options=amdgcn-amdhsa=-march=${ALF_AMD_ARCH:-gfx90a} -L${ROCM_PATH:-/opt/rocm}/lib -Wl,-rpath,${ROCM_PATH:-/opt/rocm}/lib -lhipsolver -lhipblas -lamdhip64 -lstdc++` (`configure.sh:631`) | `-DALF_GPU -DALF_LA_HIP` | `ALF_GPU=amd` |

- Any other MACHINE/GPU combination, or an unknown `GPU=` value, prints an error and `return 1` (no fallback).
- `ALF_GPU_CC` (environment, e.g. `cc80,cc90`) overrides the NVIDIA compute capability. Default `cc80` (A30/A100); `HELMAGPU` defaults to `cc90,cuda12.9`.
- **`GPU=amd`** (recipe: `platforms/amd_mi210/README.md`, `configure.sh:622-632`, `Libraries/Modules/Makefile`): needs gfortran >= 14 built with amdgcn offload and a ROCm
  installation (`ROCM_PATH`, default `/opt/rocm`; `HIPCC`, default `$ROCM_PATH/bin/hipcc`); `ALF_AMD_ARCH` (default `gfx90a`, MI200/MI210) sets the target. `make` compiles the
  hipBLAS/hipSOLVER shim `alf_la_hip_shim.cpp` (g++) and the HIP kernels `alf_hip_kernels.hip.cpp` and HIP QR `alf_hip_qr.hip.cpp` (`hipcc`) into the module library. For the speed measured
  on the MI210 the optional patched libgomp GCN plugin is needed (`platforms/amd_mi210/plugin-gcn-pool-prof.patch`, `plugin-gcn-presync-hook.patch`; without the presync hook every
  call synchronises, `ALF_HIP_ASYNC`). `Checkerboard = .T.` with real groups hangs on this backend (known issue in the AMD README); use `Checkerboard = .F.`.
- **Helma (H100)**: nvhpc **25.7** with CUDA **12.9** (`cc90,cuda12.9`). nvhpc 24.11 compiles and links but the L = 8 GPU run aborts with `la_qrpp: NaN or Inf` (not usable); without
  `cuda12.9` the default CUDA 13 device compile fails. Details: `platforms/helma_h100.md`. `ALF_HELMA_NVHPC=24.11` is kept only for reproducing that failure.
- **Ahead-of-time compilation for Intel exists**: the default `GPU=intel` build is JIT (`spir64`), startup of one process about 65 s on SuperMUC-NG Phase 2. AOT for PVC (`spir64_gen`, `ocloc`) gives
  about 14 s: replace `-fopenmp-targets=spir64` in `ALF_FLAGS_MODULES`, `ALF_FLAGS_ANA`, `ALF_FLAGS_PROG` and `ALF_LIB` after `configure.sh` by
  `-fopenmp-targets=spir64_gen -Xopenmp-target-backend -device -Xopenmp-target-backend pvc` (no quotes: the Makefiles re-split the flags; recipe `aottt`/`aot` in
  `sessions/supermuc_large_prompt/phase2/team/build_intel.sh`; the wrapper `scripts/supermuc_ng2/ifx_aot` does the same for the first-session build). The two SYCL libraries are
  built AOT with `ALF_SYCL_BLK_FLAGS='-fsycl-targets=spir64_gen -Xs "-device pvc"'` (environment or make variable; default empty = JIT; used for `libalfsb.so` and `libalfdevqr.so`).
  Results are identical to JIT, about 3 % faster than a warm JIT run (`sessions/supermuc_large_prompt/FINDINGS.md` section 4).
- `ALF_FLAGS_QRREF` (libqrref) is left unchanged. No flag contains a space (CMake testsuite splits on spaces).
- The link line gets no compile flags, hence the offload flags are also in `ALF_LIB`.

## Known limits
- (Historical, ticket 201: the GPU flags alone did not change numerics; since ticket 202 `-DALF_GPU` and `-DALF_LA_*` compile the GPU modules and the `alf_la` backends.)
- Helma (H100): `. ./configure.sh HELMAGPU PP GPU=nvidia HDF5 NO-INTERACTIVE`, then `make lib program` (login shell needed for `module`).
- MPI + GPU: `. ./configure.sh PGI MPI GPU=nvidia HDF5 NO-INTERACTIVE` builds and runs (see MPI_GPU.md); `GPU=intel` with `MPI` and `PP` is the SuperMUC-NG Phase 2 production build
  (`platforms/intel_pvc_supermuc_ng2.md`). The device of rank `r` is `r mod number of devices` with the global rank (MPI_GPU.md, "Rank-to-device binding").
- `GPU=intel` runs on a GPU only with a level-zero device; select it with `ONEAPI_DEVICE_SELECTOR=level_zero:N`
  (on the `intel` test machine `level_zero:0` is the B580, `level_zero:1` the iGPU). `OMP_TARGET_OFFLOAD=MANDATORY`
  makes a missing device an error instead of a host fallback.
- Helper scripts: `docs/gpu/scripts/nividia/env_diff.sh` (identity check), `build_gpu.sh` (build + ldd + ctest).
