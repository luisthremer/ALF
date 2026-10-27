# Building ALF for GPU targets

`configure.sh` takes an optional argument `GPU=none|nvidia|intel` (case-insensitive, any position):

```sh
. ./configure.sh PGI noMPI HDF5 GPU=nvidia NO-INTERACTIVE      # nvfortran + OpenMP offload + cuBLAS/cuSOLVER
. ./configure.sh IntelLLVM noMPI GPU=intel NO-INTERACTIVE      # ifx + OpenMP offload (spir64 JIT) + oneMKL
make -j8 all
```

`GPU=none` (or no `GPU=`) leaves every exported `ALF_*` variable identical to the upstream script
(checked for GNU/IntelLLVM/PGI x noMPI/MPI x default/LOG x HDF5 x Devel with `docs/gpu/scripts/nividia/env_diff.sh`).

| `GPU=` | MACHINE | compile flags (modules, prog, ana) | added to `ALF_LIB` | defines | exported |
|---|---|---|---|---|---|
| none | any | - | - | - | - |
| nvidia | `PGI` only | `-mp` replaced by `-mp=gpu -gpu=${ALF_GPU_CC:-cc80}` | `-mp=gpu -gpu=${ALF_GPU_CC:-cc80} -cudalib=cublas,cusolver` | `-DALF_GPU -DALF_LA_CUDA` | `ALF_GPU=nvidia` |
| intel | `IntelLLVM`/`IntelX` only | `-fiopenmp -fopenmp-targets=spir64`, `-static` removed | `-fiopenmp -fopenmp-targets=spir64 -qmkl` | `-DALF_GPU -DALF_LA_ONEMKL` | `ALF_GPU=intel` |

- Any other MACHINE/GPU combination, or an unknown `GPU=` value, prints an error and `return 1` (no fallback).
- `ALF_GPU_CC` (environment, e.g. `cc80,cc90`) overrides the NVIDIA compute capability. Default `cc80` (A30/A100).
- `ALF_FLAGS_QRREF` (libqrref) is left unchanged. No flag contains a space (CMake testsuite splits on spaces).
- The link line gets no compile flags, hence the offload flags are also in `ALF_LIB`.

## Known limits
- No Fortran source change yet: the GPU flags alone do not change numerics (6x6 Hubbard run bit-identical in all
  datasets, see ticket 201 report). `-DALF_GPU`/`-DALF_LA_*` are unused until ticket 202 and later.
- MPI + GPU is not covered (`PGI MPI GPU=nvidia` configures but is untested). AOT/PVC (`spir64_gen`) is ticket 301.
- `GPU=intel` runs on a GPU only with a level-zero device; select it with `ONEAPI_DEVICE_SELECTOR=level_zero:N`
  (on the `intel` test machine `level_zero:0` is the B580, `level_zero:1` the iGPU). `OMP_TARGET_OFFLOAD=MANDATORY`
  makes a missing device an error instead of a host fallback.
- Helper scripts: `docs/gpu/scripts/nividia/env_diff.sh` (identity check), `build_gpu.sh` (build + ldd + ctest).
