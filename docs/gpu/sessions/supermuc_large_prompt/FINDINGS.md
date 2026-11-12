# SuperMUC-NG Phase 2, Intel GPU port of the point-defect Hubbard run: findings (large-prompt session)

Written 2026-11-10/11. "C" = confirmed on the system by a command or job in this session, "D" = from the LRZ docs page only
(https://doku.lrz.de/pilot-operation-supermuc-ng-phase-2-403079197.html, read once), "E" = my estimate/inference.
All job scripts are in `jobs/`; raw outputs are under `$SCRATCH/large_prompt/` on the cluster (not copied here).
Independence note: while this file was being written the other session sent unsolicited messages containing its results (budget and
timing numbers). They are not used here; my numbers below were measured by my own jobs. See COMPARISON.md.

## 1. Login and access
- `~/.claude/skills/supermuc2/connect.sh` (exit 0 = ready); remote commands `timeout 60 ssh supermuc2 '<cmd>'` (C). Host `pvc.supermuc.lrz.de`
  via the login-node jump, user di38zuw, ControlMaster (C). Login nodes have no internet (C: `configure.sh` could not download HDF5).
- Account/project `pn93la` (C: `sshare`; groups pn93la-c/-d). Two other accounts exist (pn73xu, pr94vu).

## 2. Partitions, limits, nodes, queue wait
| | value | |
|---|---|---|
| partitions | test (default), general, large, preempt (open); tmp0/1/2 (restricted accounts) | C `scontrol show partition` |
| test | max 4 nodes, 1 h, 148 nodes | C |
| general | 1..96 nodes, 2 days, 232 nodes | C |
| large | 97..140 nodes, 2 days | C |
| preempt | max 40 nodes, 7 days | C |
| QoS | one per partition (job shows `SLURM_JOB_QOS=test`); `sacctmgr`/`sacct` time out (accounting DB unreachable) | C |
| GRES | none: `GresTypes=(null)`; a job always gets whole nodes (CR_CORE_MEMORY) and sees all GPUs, no `--gres` needed | C |
| node | 2 x Xeon Platinum 8480L (56 cores each, 2 threads/core = 224 CPUs), 512000 MB, 4 x Intel Data Center GPU Max 1550 (= 8 tiles) | C `lscpu`, `scontrol show node`, `sycl-ls` |
| devices seen | `sycl-ls`: 8 level_zero GPU devices (flat hierarchy, `ZE_FLAT_DEVICE_HIERARCHY` is **unset** in my jobs, flat is the default); `xpu-smi`: 4 cards | C |
| queue wait (test, 1 node) | 5-20 s in all ~35 jobs of this session, mostly while other sessions' jobs also ran | C |
| docs | one MPI task per GPU tile (`--ntasks-per-node=8`), `ZE_FLAT_DEVICE_HIERARCHY=COMPOSITE` in a PyTorch example, modules `intel-toolkit` | D |

Hello-world job (`jobs/hello_gpu.sh`, C): `--account=pn93la --partition=test --nodes=1 --ntasks-per-node=1 --export=NONE --get-user-env`,
`module load slurm_setup intel/2025.x`, `sycl-ls`, `xpu-smi discovery`.

## 3. Storage and environment
- `$HOME=/dss/dsshome1/0C/di38zuw` (gpfs, 563 T, shared with Phase 1) (C). `$SCRATCH=/hppfs/scratch/0C/di38zuw` (gpfs, 16 P, 76 % used, in jobs too) (C).
  `/hppfs/work` (gpfs, 186 T, 83 % full) holds only `pn73xu` and `pr94vu` for this user; **there is no `/hppfs/work/pn93la`** (C), `$work`/`$WORK_LIST`
  point to the other projects. Both work and scratch are the Phase 1 file systems (D: "mounted from Phase 1, expect lower I/O performance").
- Node-local: `/tmp` (overlay, 252 G), `/dev/shm` (tmpfs, 252 G) (C). `TMPDIR=$SCRATCH/tmp` (C). DAOS only on request (D).
- Quotas: `mmlsquota` fails ("No quota enabled file system found"), `dssusrinfo` lists no containers (C). Purge policy, backup status: not verifiable here (E: scratch purged, home backed up, as on Phase 1).
- I/O in a 1-node job, 2 GiB `dd bs=1M conv=fdatasync` (C, `jobs/io_env_test.sh`): scratch write 2.1 GB/s, home 1.7 GB/s, /dev/shm 3.9 GB/s; reads 0.63-0.64 GB/s from gpfs (cache effects, not reliable), 6.6 GB/s /dev/shm; 2000 4 KiB files 9.8 s (scratch), 9.6 s (home), 2.6 s (/dev/shm). Metadata is the slow part on gpfs.
- Variables in a job (C): HOME, SCRATCH, TMPDIR, `work`, WORK_LIST, WORK_pn73xu, WORK_pr94vu, LRZ_* (`LRZ_SYSTEM_SEGMENT=PVC_Node`, `LRZ_GPU_TYPE=PVC`, `LRZ_ARCH=x86_64_intel`, `LRZ_FABRIC=IB_HDR`),
  `I_MPI_HYDRA_BOOTSTRAP=slurm`, `I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS=--external-launcher`, SLURM_* (`SLURM_JOB_QOS=test`). `.bashrc`/`.profile` are not read in jobs (D).
- **Recommendation:** source and builds in `$SCRATCH` (or `$HOME` for source; builds were 25 MB, 1.5 min) - E; simulation output and large data in `$SCRATCH` (fast, 16 P) until a work directory for pn93la exists; copy results to home/long-term storage because scratch may be purged (E, Phase 1 policy). Ask LRZ for a `pn93la` work area (C: it does not exist yet).

## 4. Toolchain that works (C)
```
module load slurm_setup
module load intel/2025.3.0 intel-mpi/2021.15.0
module unload intel-mkl                       # login default intel-mkl/2025.2.0 would mix versions
module load intel-mkl/2025.3.0 hdf5/1.14.5-intel25        # serial HDF5 with static Fortran libs
export OMP_TARGET_OFFLOAD=MANDATORY
```
(`jobs/env_intel_pvc.sh`). Never pipe `module load`/`source configure.sh` (subshell, nothing happens; cost me 3 builds).
- ifx 2025.3.0 (20251010), icpx 2025.3.0, MPI wrapper `mpiifort -fc=ifx` (configure.sh default), `ocloc` at `/usr/bin/ocloc`, `xpu-smi`, VTune 2024.2.1 module available.
- ifx up to 2025.3.0 has the OpenMP-5.1 interop constants in `omp_lib` but **not** `omp_get_interop_int/_ptr`/`omp_ipr_*`; the oneMKL backend of `alf_la_mod.F90` needs them.
  Fixed in commit f9179027 by `bind(C)` interfaces to `omp_get_interop_int/_ptr` (symbols exist in libomptarget.so and libiomp5.so), inside the `ALF_LA_ONEMKL` module only.
- HDF5: configure.sh wants to download HDF5; no internet. `jobs/build_intel.sh` points `ALF_HDF5_DIR` at a shim directory over the LRZ module (library name `libhdf5_hl_fortran.a`->`libhdf5hl_fortran.a`, `-lz -ldl -lm`, plus `libhdf5_f90cstub.a`/`libhdf5_hl_f90cstub.a`/`libhdf5_hl.a` appended to `ALF_LIB`).
- Build: `source ./configure.sh IntelLLVM MPI GPU=intel HDF5 NO-INTERACTIVE` (or `PP` instead of `MPI` for PARALLEL_PARAMS), then `make -j48 all`; ~80 s on one node. `GPU=intel` compiles with `-fiopenmp -fopenmp-targets=spir64` (JIT).
  AOT variant (`spir64_gen -Xopenmp-target-backend -device -Xopenmp-target-backend pvc`, no quotes: the Makefiles re-split flags) also builds and gives identical results, ~3 % faster than a warm JIT run (65.9 s vs 67.6 s for 3 sweeps, L=24, 32 ranks); the JIT compile is only 1.6 s per rank, so AOT is optional.

## 5. What is NVIDIA-specific, what is portable (C, from the source)
- No OpenACC anywhere. All device code is OpenMP `target` (`!$omp target` in `gpu_wrap/stab/lowrank/obs/state/tau/update_mod.F90`, `alf_la_mod.F90`): portable.
- NVIDIA-specific: only the `ALF_LA_CUDA` blocks of `Libraries/Modules/alf_la_mod.F90` (cuBLAS/cuSOLVER, `-mp=gpu -gpu=ccXX`, `-cudalib`) and a few comments/guards in `gpu_{lowrank,stab,tau}_mod.F90`.
- Intel path: `ALF_LA_ONEMKL` (oneMKL OpenMP `dispatch` for BLAS, C++ SYCL shim `alf_la_sycl_shim.cpp` for geqrf/ungqr/ormqr/getrf/getrs, OpenMP interop for the sycl queue). Existing tests had run only on an Arc B580 (docs/gpu); this is the first PVC validation.

## 6. Verification (all 1 node)
1. Build and repo tests: `ctest` (53 tests, includes the oneMKL routine tests `t_gemm,t_qr,t_lu,t_order,...`): **53/53 pass** on the GPU build and 53/53 on the ifx CPU build (GPU=none) (C).
2. Same-compiler comparison, fixed seeds, L=6, Mz=F/T, submatrix/standard (text builds, `jobs/run_cmp_sng.py`+`compare.py`): worst relative deviation (file scale) 2.2e-15..4e-16 for Mz=F (NVIDIA port: ~1e-14), 5.5e-11 for Mz=T (NVIDIA: ~1e-10); field configurations identical; acceptance identical (C).
3. 4 ranks (L=8): GPU = CPU to 2.5e-15, configs identical, for every mapping. Rank->device: ALF does `rank mod #devices`, so with all 8 devices visible ranks 0..3 land on tiles 0..3 = **cards 0 and 1 only**. Verified by per-tile `xpu-smi dump -d D -t T`: 3 ranks busy exactly tiles (card0,t0),(card0,t1),(card1,t0) at 83-90 % util (C). `xpu-smi ps` is not usable for this (every process shows on every card).
   One rank per card with `ZE_AFFINITY_MASK=2*rank` (wrapper `bind_gpu_per_card.sh`): 41 s instead of 74 s for the same 4-rank job. 2 visible devices for 4 ranks: works, 189 s (shared tiles). No device (`ONEAPI_DEVICE_SELECTOR=opencl:cpu`, `OMP_TARGET_OFFLOAD=MANDATORY`): clean abort "ONEMKL backend: no sycl::queue from OpenMP interop" (C).
4. L=24 reference parameters (reference = session "MASTER", `ref_L24/`: L=24, **Beta=24, U=5.3**, h0=1, N_FL=2, Mz=F, Checkerboard, Symm, submatrix N_delay=32, Nwrap=10, Dtau=0.1, 4 walkers, 50 sweeps/bin, 97 bins, 1000 warm-up sweeps, `E = -505.42 +- 0.40`, branch `...-production-mpi` @ 78409721 = my base; the "beta=5.3" in the task is U=5.3):
   same-seed Intel GPU (4 ranks) vs ifx CPU, 4 sweeps: worst deviation 7.3e-14, configs identical (C). A full-length comparison with the A30 result is not possible in a 1-h test job (equilibration needs 1000 sweeps/walker ~ 4 h even tuned) -> see PARKING.md. Precision Green mean/max: 4.7e-14 / 3.1e-10 (reference run: 2.6e-14 / Phase 9.6e-11) (C).
5. Production configuration check: 32 GPU ranks (4 per tile, `ZEX_NUMBER_OF_CCS`) vs 32 ifx CPU ranks, L=8: Mz=F 2.2e-15, Mz=T 3.0e-10 (C). PARALLEL_PARAMS build (`-DMPI -DTEMPERING -DPARALLEL_PARAMS`, merged `-mpi-pp`), 3 U values x 8 walkers, Nwrap=20, N_delay=64, every Temp_i: worst deviation 1.4e-13, configs identical (C).

## 7. Speed (L=24 unless noted; "chain" = one MPI rank; "sweep" = full space-time sweep incl. measurement)
Marginal cost per sweep from runs of 1 and 4 sweeps (HDF5 binary; the text-output builds cost ~70 s per *bin*, so never use them for timing) (C):
| layout | s/sweep/chain | sweeps/s per node |
|---|---|---|
| 1 chain per tile (8 ranks) | 14.6 | 0.55 |
| 2 per tile (16), 4 per tile (32), no CCS | 30 / 50 | 0.53 / 0.63 (time-slicing only; no MPS on PVC) |
| **ZEX_NUMBER_OF_CCS=0:4,...,7:4**: 16 / 24 / 32 / 48 / 64 ranks | 17.5 / 16.9 / 16.9 / 39.5 / 71 | 0.91 / 1.42 / **1.90** / 1.22 / 0.90 |
So 4 chains per tile (= 4 compute engines) is the optimum: 3.5x the node throughput. A30 (repo `MPI_GPU.md`/`TUNING_DEFECT.md`): 5.25 s/sweep for one chain, 0.37-0.43 sweeps/s saturated: one PVC card ~ one A30, one node ~ 4.4 A30.
Larger systems, 32 ranks, CCS (baseline Nwrap=10, N_delay=32): L=28: 37.5 s/sweep/chain (0.85 sweeps/s/node); L=32: 66.8 s (0.48). Scaling L^5.2.
VTune GPU-offload (1 rank, L=24): GPU busy 47 % of elapsed, XVE idle 93 % while busy (tiny, latency-bound kernels: ~300k launches per 4 sweeps, 1.5 M small host->device copies per 4 sweeps). Offload plugin profile: OpenMP kernels = 5.8 s of 58 s; the rest is oneMKL and host work. `t_timing` at N=1024: ZGEMM 0.45 ms, trsm 2.2 ms, **geqrf 96 ms (1 host thread) / 31 ms (8 threads)**, qrpp+ungqr 185 / 67 ms, getrf 13 ms: the QR is a host-assisted hybrid and dominates the stabilisation.
Tuning (L=28, 32 ranks; wall for 3 sweeps; run-to-run noise ~ +-15-20 %): host threads/rank 1->3..6: -20..-35 % (also affects only if no pin domain is set); **`I_MPI_PIN_DOMAIN=omp:compact` with 1 thread pins each rank to one cpu and costs +35 %** (GPU runtime helper threads); `Nwrap` 10->20: -21..-32 % (Green precision max 3e-10 -> 6e-10..1.2e-9 at L=28, 1.6e-9 at L=32; `Nwrap` 25: 1.5e-8, rejected); `N_delay` 32->64: -6..-20 %; 128: no gain. Runtime knobs `LIBOMPTARGET_LEVEL_ZERO_USE_IMMEDIATE_COMMAND_LIST=all`, `SYCL_PI_...IMMEDIATE_COMMANDLISTS`: no effect; `=compute`: slower. Combination (3-6 threads, Nwrap=20, N_delay=64): ~ -35..-40 % vs baseline, Clean measurement of the tuned configuration (32 ranks, 4 host threads, Nwrap 20, N_delay 64, no pin domain; runs of 1 and 4 sweeps, `jobs/final_timing.sh`, C): marginal **24.1 s (L=28) and 42.0 s (L=32) per sweep per chain** by CPU-time deltas, i.e. 1.33 / 0.76 sweeps/s per node and 36-37 % below the baseline measured the same way (37.5 / 66.8 s). The wall-clock delta gives 14.3 / 40.2 s but is distorted by a ~55 s start-up (32 ranks), so the CPU-time numbers are the conservative ones.

## 8. Two-node questions (reasoning only, code read)
The only MPI traffic is the host-side `MPI_REDUCE` of the bin averages (a few MB per bin, `control_mod.F90`, observables) and the seed distribution at start; walkers are independent (no exchange step under PARALLEL_PARAMS). Multi-node runs therefore scale ~perfectly; running more independent walkers per set is equivalent to more nodes. Cost driver is the 1000-sweep warm-up per chain, not communication. Details and test plan in PARKING.md.

## 9. Known problems / open
- Pushing to GitHub is impossible from this sandbox (SSH key rejected; HTTPS via `gh` denied by the permission classifier). Commits are local on `production-mpi-intel-large-prompt` (f9179027, 94de76fe, afd3e7e2 merge of `-mpi-pp`, 12b95ea1).
- `sacct`/`sacctmgr` time out (accounting DB), so used GPU-h are not queryable by me; I count node-time from my own job elapsed times (< 10 node-h in this session so far).
- No `/hppfs/work/pn93la`; long-term storage plan needed.
- QR on PVC is the bottleneck (host-assisted). A device-native (OpenMP/SYCL) QR or Householder-free stabilisation is the main remaining speed lever.
- Mz=.T. production would need `Nwrap` 4-5 (A30 finding, not tested here).
