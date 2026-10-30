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

**State 2026-10-29 (read-only checks):** compute mode **Default**, MIG disabled; `/` 19 GB free (93 %); NVIDIA MPS
binaries present (`/usr/bin/nvidia-cuda-mps-control`, `-server`), no daemon running; other users log in at times
(no GPU use seen). Two CUDA processes at once are time-sliced (13× slower each, DISCREPANCIES.md); MPS trial results:
gate report of 2026-10-29 and CLAUDE.md "nividia GPU sharing".

## `intel` (host `IntelGPU`) – probed 2026-10-27, read-only

Reached via `ssh intel` (jump host). **Test sandbox on a shared workstation** (other users and a desktop session).
Rules: CLAUDE.md "Machine rules". Our dir `~/alf-gpu/` (mode 700).

| Item | Value |
|---|---|
| OS | Ubuntu 26.04.1 LTS (was 25.04 until the admin's update of 2026-10-29) |
| CPU | Intel Core i5-14400 (6P + 4E cores, 16 threads), 1 socket, 61 GB RAM |
| GPU | **Intel Arc B580** (Battlemage G21, Xe2, discrete, `/dev/dri/renderD128`); iGPU UHD 730 (`renderD129`) |
| GPU driver | `intel-opencl-icd` 26.35.39758.17, `libze-intel-gpu1` 26.35.39758.17, Level Zero loader `libze1` 1.32.0, `libigc2` 2.41.15 (checked 2026-10-30; 24.52 / 1.20.3 / 2.5 until 2026-10-29) |
| oneAPI | `/opt/intel/oneapi`: compiler 2025.1, 2025.3 (ifx 2025.3.3), **2026.2 (ifx 2026.2.0)**; MKL **2026.2** only (`mkl/2025.1` is an empty dir); Intel MPI 2021.15–2021.19; TBB 2023.2, UMF 1.2, TCM 1.6. **No `apt-mark hold`** on any intel package (checked 2026-10-30) |
| `ocloc` | system-wide → AOT `spir64_gen` builds possible here (not on nividia) |
| Other | gfortran, cmake, python3, clinfo; ample disk |

Environment: `. ~/alf-gpu/env/intel.sh`, pinned per component since 2026-10-30: `tcm/1.6 umf/1.2 tbb/2023.2
mpi/2021.19 compiler/2026.2 mkl/2026.2` (each `env/vars.sh`; never `latest`). Gives ifx 2026.2.0 20261001,
`MKLROOT=/opt/intel/oneapi/mkl/2026.2`, `mpiifx` 2021.19. The old script is kept as `intel.sh.2025.3`.
oneAPI's `vars.sh` parses and resets the positional parameters of whoever sources it, so until 2026-10-30 sourcing
`intel.sh` inside a script emptied the caller's `"$@"` (found by the 302a implementer: `-V`/`-R` silently dropped; the
harness's `acc_ctest.sh` sources the env file inside `bash -c` and then reads `$1`–`$3`). Since 2026-10-30 the script sources
each `vars.sh` inside a function with an empty argument list; the resulting environment is identical (`env | sort`
diff empty), the previous file is kept as `intel.sh.argreset_2026-10-30`. The nividia env files do not have this problem
(checked 2026-10-30). No acceptance run had used `intel` before the fix.
GPU runs on `intel` wrote the Intel GPU compiler cache to `~/.cache/neo_compiler_cache` (16 MB by 2026-10-30, from the
302a runs), outside `~/alf-gpu/` (found by the 302b draft). Since 2026-10-30 `intel.sh` sets `NEO_CACHE_PERSISTENT=0`
(kernels are JIT-compiled at each start; previous file `intel.sh.neocache_2026-10-30`). The existing cache directory is
left in place until the owner decides.
Checked 2026-10-30: `ifx -fiopenmp -fopenmp-targets=spir64 -qmkl ... -lmkl_sycl_blas -lmkl_sycl_lapack -lsycl` links
(the link that failed on 2026-10-29); ZGEMM plus a target region runs on the B580 under `OMP_TARGET_OFFLOAD=MANDATORY`
with the correct result; `ldd`: all MKL/SYCL/iomp libraries from 2026.2. Sourcing prints a harmless UMF warning on
stderr.

**No package hold (owner, 2026-10-30, DECISIONS.md row 13):** before every Intel-lane T1/T2 run, record
`dpkg-query -W 'intel-oneapi-compiler-fortran*' 'intel-oneapi-mkl*' intel-opencl-icd libze-intel-gpu1 libze1 libigc2`
into the run directory and compare it with the ticket's recorded set; on any difference stop and report.

**Fixed 2026-10-29 evening by the admin** (a matching 2026.2 set instead of decision 6's MKL 2025.3 + hold; no hold, so a
later system update can mix versions again). Earlier state, kept for the record:
**State 2026-10-29: the oneAPI installation is mixed and `GPU=intel` does not link.** A system update (owner's
maintenance, 2026-10-28) removed MKL 2025.3 (`dpkg` state `rc`) and installed MKL 2026.2, without a Fortran compiler
2026.2: `compiler/` has 2025.1, 2025.3 (ifx 2025.3.3), 2026.2 (icx/icpx, no ifx); `mkl/` has 2025.1 (empty) and 2026.2.
`intel.sh` no longer finds MKL; pairing ifx 2025.3 with MKL 2026.2 fails at link time (`sycl::_V1` symbols from
`libmkl_sycl_*.so` 2026.2 need the 2026.2 SYCL runtime). Intel supports MKL only with the matching compiler release.
Snapshot: `~/alf-gpu/backup/oneapi_state_2026-10-29.txt`. Fix options: gate report of 2026-10-29 (step 5).

**GPU access (fixed 2026-10-27 by the owner: account added to `render`).** Checked the same day:
- Level Zero sees `level_zero:0` = Arc B580 (`[0xe20b]`, driver 1.6.32224) and `level_zero:1` = UHD 730. Select
  the B580 with `ONEAPI_DEVICE_SELECTOR=level_zero:0`.
- B580 aspects include `fp64` and `atomic64`. An ifx 2025.3 OpenMP target region (complex(8) multiply, `exp`,
  real(8) `sqrt`) runs on the device under `OMP_TARGET_OFFLOAD=MANDATORY`, max rel. diff vs host 4.4e-16, both JIT
  (`-fopenmp-targets=spir64`) and AOT (`spir64_gen -Xopenmp-target-backend "-device bmg-g21"`, i.e. ocloc works).
- **FP64 rate not measured yet** (a ZGEMM/DGEMM benchmark; will come with the first oneMKL work).
- **Don't run plain `sycl-ls`:** it lists the devices and then aborts in the driver (`drm_neo.cpp:410`), and Ubuntu's
  apport writes a ≈ 110 MB report to `/var/crash/_opt_intel_oneapi_compiler_2025.3_bin_sycl-ls.<uid>.crash` (one
  file, overwritten on repeat). One exists from this probe; it is outside `~/alf-gpu/`, so agents leave it to the
  owner. `sycl-ls --verbose` with the selector did not abort.

**Machine split (owner, 2026-10-27):** Intel targets (`GPU=intel`, ifx offload, oneMKL/SYCL) are tested on
`intel`; NVIDIA targets on `nividia`. CPU-only checks (unit tests, `GPU=none` identity, references) stay on nividia.

## SuperMUC-NG Phase 2 (LRZ) – not reachable by agents
2× Xeon Platinum 8480+, 4× Intel Data Center GPU Max 1550 (2 stacks each) per node, 240 nodes.
oneAPI module versions, Slurm partitions and limits: [verify, owner].
