# Machines

## `nividia` (host `gpuNV-test`) – probed 2026-10-26, read-only + 10 s smoke tests

| Item | Value |
|---|---|
| OS | Ubuntu 24.04.5 LTS, kernel 6.8.0-147 |
| CPU | AMD EPYC 7763 (Zen 3), **16 cores** visible (VM slice, 1 socket, 1 NUMA node, SMT off), AVX2, no AVX-512 |
| RAM / disk | 31 GB / `/` 227 GB, **25 GB free (90 % used)** |
| GPU | 1× **NVIDIA A30**, 24 GB HBM2, compute capability 8.0 (sm_80), 165 W, ECC on, MIG disabled |
| Driver | 615.71.10, CUDA UMD 13.4 |
| CUDA toolkit | 13.2 at `/usr/local/cuda` (not on PATH); NVHPC ships CUDA 13.3 |
| NVHPC | **26.9** at `/opt/nvidia/hpc_sdk/Linux_x86_64/26.9` (not on PATH): nvfortran/nvc/nvc++, cuBLAS, cuSOLVER, NVHPC BLAS/LAPACK, HPC-X OpenMPI, NCCL, NVSHMEM |
| GNU | gfortran/gcc 13.3.0; OpenBLAS (pthread) as system `libblas.so.3`, reference `liblapack` |
| MPI on PATH | none. Available: NVHPC HPC-X (`comm_libs/mpi/bin/mpifort`), Intel MPI 2021.19 |
| Intel oneAPI | `/opt/intel/oneapi` (system-wide, many versions). **ifx only in `compiler/2025.3`** (ifx 2025.3.3); `compiler/latest` = 2026.2 has icx/icpx but no ifx. MKL 2026.2 (incl. SYCL BLAS/LAPACK), Intel MPI 2021.19 (`mpiifx`). **No `ocloc`** |
| Build tools | cmake 3.28.3, make, git (GitHub over https works) |
| Profilers | perf (`perf_event_paranoid=4`, so unprivileged perf is restricted [verify]), gprof, nsys at `/usr/local/bin/nsys`, ncu not on PATH [verify in NVHPC `profilers/`] |
| Python | 3.12.3, venv works, **no pip, numpy, scipy, h5py or pyALF** |

### Smoke tests (all cleaned up afterwards)

| Test | Result |
|---|---|
| nvfortran 26.9 `-mp=gpu -gpu=cc80`, complex(8) `target teams distribute parallel do`, `OMP_TARGET_OFFLOAD=MANDATORY` | PASS, ran on device, correct |
| nvfortran + `use cublas`, `cublasZgemm` on OpenMP-mapped arrays via `use_device_addr` | PASS, correct (0,512) |
| Same with plain `zgemm` (host BLAS symbol) on device addresses | link fails / crashes, as expected: device calls must go through an explicit backend, never the host BLAS name |
| ifx 2025.3 `-fiopenmp -fopenmp-targets=spir64` (JIT), `OMP_TARGET_OFFLOAD=DISABLED` | PASS, host fallback correct |
| Same, `OMP_TARGET_OFFLOAD=MANDATORY`, no Intel GPU present | aborts (core dump), so a host fallback fails loudly, as wanted on SuperMUC |
| ifx `-fopenmp-targets=spir64_gen -Xopenmp-target-backend "-device pvc"` (AOT) | **FAIL: `ocloc` missing** (part of Intel's GPU driver stack, system package) |
| ifx `-qmkl` host `zgemm` in an offload build | PASS |

### Environment snippets

```sh
# NVIDIA
export NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.9
export PATH=$NVHPC/compilers/bin:$NVHPC/comm_libs/mpi/bin:$PATH
# Intel (ifx lives only in 2025.3; setvars picks 'latest' = 2026.2 for the rest)
source /opt/intel/oneapi/setvars.sh --force >/dev/null
export PATH=/opt/intel/oneapi/compiler/2025.3/bin:$PATH
```

### Consequences
- User-space oneAPI install is **not needed**: ifx, MKL and Intel MPI are already there. A 2025.3 ifx with
  a 2026.2 MKL/runtime is a version mix; SuperMUC's module versions must be checked [verify, owner].
- AOT `spir64_gen` for PVC can't be tested here without `ocloc` (system install, needs the owner's OK).
  JIT `spir64` compiles and catches front-end/offload-lowering bugs.
- Upstream ALF has no nvfortran profile in `configure.sh` (R1 §7: only `PGI`) [verify in step 3].
- pyALF analysis and `test_vs_ed` need numpy/scipy/h5py + pyALF in a venv (pip bootstrap via `python3 -m venv`).
- The CPU reference is 16 Zen 3 cores, not a full HPC node; BENCHMARKS.md must say so.
- Disk is tight: keep builds and run outputs small, delete after use.

## `intel` (host `IntelGPU`) – probed 2026-10-27, read-only

Reached via `ssh intel` (jump host). **Test sandbox on a shared workstation** (other users and a desktop session).
Rules: CLAUDE.md "Machine rules". Our dir `~/alf-gpu/` (mode 700).

| Item | Value |
|---|---|
| OS | Ubuntu 25.04 |
| CPU | Intel Core i5-14400 (6P + 4E cores, 16 threads), 1 socket, 61 GB RAM |
| GPU | **Intel Arc B580** (Battlemage G21, Xe2, discrete, `/dev/dri/renderD128`); iGPU UHD 730 (`renderD129`) |
| GPU driver | `intel-opencl-icd` 24.52, `libze-intel-gpu1` 24.52, Level Zero loader 1.20.3, `libigc2` 2.5 |
| oneAPI | 2025.1 … 2026.1 in `/opt/intel/oneapi`. **ifx only up to 2025.3** (2025.3.3). MKL per release, Intel MPI 2021.15–2021.18 |
| `ocloc` | system-wide → AOT `spir64_gen` builds possible here (not on nividia) |
| Other | gfortran, cmake, python3, clinfo; ample disk |

Environment: `. ~/alf-gpu/env/intel.sh` (= `/opt/intel/oneapi/2025.3/oneapi-vars.sh --force`; ifx 2025.3.3,
`MKLROOT=/opt/intel/oneapi/2025.3`, `mpiifx`).

**Blocker (2026-10-27): no GPU access for our account.** The render nodes are `root:render 0660` and the account is
not in `render`, so `sycl-ls` and `clinfo` show only the OpenCL CPU device. Fix by the owner (root): add the account
to `render`, then log in again. Agents never do this themselves.

To verify once access works: `sycl-ls` lists `level_zero:gpu` for the B580; **FP64 support and rate** (`clinfo`
"Double-precision Floating-point support", and a measured DGEMM/ZGEMM rate). ALF needs complex(8) throughout, so a
weak FP64 rate limits what the B580 can say about performance; correctness and the ifx offload/oneMKL/SYCL paths are
what this machine is for.

## SuperMUC-NG Phase 2 (LRZ) – not reachable by agents
2× Xeon Platinum 8480+, 4× Intel Data Center GPU Max 1550 (2 stacks each) per node, 240 nodes.
oneAPI module versions, Slurm partitions and limits: [verify, owner].
