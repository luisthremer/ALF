# GPU=intel on SuperMUC-NG Phase 2 (Intel Data Center GPU Max 1550)

Status 2026-11-10, branch `intel-supermuc-ng2` (production-mpi + `-mpi-pp`). Everything below was measured on the system
(project `pn93la`, partition `test`), helper scripts in `docs/gpu/scripts/supermuc_ng2/`.

## System (from `sinfo`, `scontrol`, `xpu-smi`)
- Login `pvc.supermuc.lrz.de` (no internet on the login nodes). Partitions: `test` (<= 4 nodes, 1 h), `general` (<= 96 nodes, 2 d),
  `large` (97..140 nodes, 2 d), `preempt` (<= 40 nodes, 7 d). No micro queue, no `--gres`: whole nodes only.
- Node: 2 x Xeon 8480+ (112 cores, 224 threads), 512 GB, 4 x Max 1550 = 8 tiles (`ZE_FLAT_DEVICE_HIERARCHY=FLAT`: 8 devices, driver NEO 23.43, Level Zero 1.3).
- Budget is in GPU-h: 1 node-hour = 4 GPU-h; `pn93la` has 900 GPU-h (`budget_and_quota` after `module load lrztools`), no core-h.
- Storage: `$HOME` = `/dss/dsshome1/0C/di38zuw` (100 GB quota), `$SCRATCH` = `/hppfs/scratch/0C/di38zuw`. WORK is mounted from Phase 1 (`/hppfs/work/<project>`),
  `/hppfs/work/pn93la` did not exist yet (provisioning of a new entitlement takes up to one working day); `$WORK` is not set on the login node.

## Toolchain that works
`intel-toolkit/2025.1.0` (ifx 2025.1.1, Intel MPI 2021.15, MKL 2025.1) + `hdf5/1.14.5-intel25-impi` (`scripts/supermuc_ng2/env.sh`).
- oneAPI 2025.3 compiles but the OpenMP-to-SYCL interop crashes at run time with the Level Zero V2 adapter, and without `SYCL_UR_USE_LEVEL_ZERO_V2=1`
  Level Zero does not initialise (`sycl-ls` fails). 2025.1 works without any such variable. ifx < 2026 has no `omp_get_interop_int/ptr` in `omp_lib`
  (libiomp5 exports them): `Libraries/Modules/alf_la_mod.F90` declares them with `bind(C)` (the only source change needed).
- HDF5: ALF downloads and builds its own HDF5, impossible here. `scripts/supermuc_ng2/build.sh` points `ALF_HDF5_DIR` at the LRZ module through a shim directory
  (CMake-style library names, `h5fc` wrapper, `-lhdf5_f90cstub -lhdf5_hl_f90cstub`).
- **AOT for PVC**: `ocloc` is at `/usr/bin/ocloc`. configure.sh cannot carry a flag with a space, so `scripts/supermuc_ng2/ifx_aot` (copy to a directory first in `PATH`
  as `ifx`, `REAL_IFX` = the real ifx; `build.sh` does this for `AOT=1`) rewrites `-fopenmp-targets=spir64` into
  `-fopenmp-targets=spir64_gen -Xopenmp-target-backend "-device pvc"`. JIT start-up of one process is about 65 s, AOT 14 s.
- Builds: `AOT=1 bash -l build.sh prodaot IntelLLVM MPI GPU=intel HDF5 NO-INTERACTIVE` (production), `... ppaot IntelLLVM PP GPU=intel HDF5 ...` (parameter scans),
  `... gpuaot IntelLLVM noMPI GPU=intel ...` and `cpuref2 IntelLLVM noMPI ...` (same-compiler CPU reference, text output).
- Run: `OMP_TARGET_OFFLOAD=MANDATORY`, `ZE_FLAT_DEVICE_HIERARCHY=FLAT`; the device of rank r is `r mod 8` = tile r, so 8 ranks fill all 4 cards (`xpu-smi dump`: 67-86 % utilisation on every card).
  With fewer than 8 ranks per node the ranks land on the first cards: bind one rank per card with `ZE_AFFINITY_MASK=2*local_rank`.

## Numerical verification (same-compiler CPU reference, `run_cmp.py`, fixed seeds, `verify_vs_cpu.sbatch`)
Worst relative deviation of all observable files, acceptance identical in every case (criterion 1e-8):
Mz=F: standard 3.0e-15 (L=6), submatrix N_delay 8 / 32 / 128: 3.0e-15 (L=6), 8.0e-15 (L=8), 2.1e-14 (L=12, 32), 2.7e-14 (L=12, 128).
Mz=T (intrinsic precision of G): 2.5e-11 (L=6, both schemes), 3.4e-9 (L=12, N_delay 128); Mz=T is not used in production.
The MPI + HDF5 production binary (1 rank) reproduces the text output of the same seeds (Ener_scal bins -33.9246, -33.1407).
Precision Green at L=24 / 28 / 32: mean 1.7e-14 / 4.6e-14 / 1.4e-13, max 1.2e-11 / 1.1e-10 / 5.3e-10. Acceptance 0.324.

## Performance (L = 24, Beta = 24, N_FL = 2, one node, seconds per sweep and chain unless noted)
| change | per chain | node aggregate |
|---|---|---|
| first working build (JIT, `OMP_NUM_THREADS=1`, `N_delay` 32), 8 ranks | 14.6 | 0.55 sweeps/s |
| AOT instead of JIT | same (start-up 65 s -> 14 s) | |
| 14 host threads per rank (host MKL zgemm of the submatrix flush, 112 cores / 8 ranks) | 9.7 | 0.83 |
| `N_delay` 32 -> 128 (64: -24 %, 128 and 192: -29 %, repeatable to 0.1 s) | 8.2 | **0.97** |
| oversubscription (ranks per tile): 8 / 16 / 32 ranks, 1 thread each | 14.6 / 29 / 51 | 0.55 / 0.54 / 0.63 |
| 16 ranks x 7 threads | 27 | 0.59 |
| `LIBOMPTARGET_LEVEL_ZERO_USE_IMMEDIATE_COMMAND_LIST=all` / `compute`, `UR_L0_USE_IMMEDIATE_COMMANDLISTS=2`, `ZEX_NUMBER_OF_CCS=0:4` | no gain (first hint of -16 % did not repeat) | |

No MPS on PVC: ranks sharing a tile are time-sliced, the aggregate stays flat, so one rank per tile. The A30 reference (4 walkers, MPS) is 0.37 sweeps/s per GPU,
a PVC card (2 tiles) now reaches 0.24 sweeps/s, a node (4 cards) 0.97. After the fixes a stack sample of one rank shows ~60 % waiting for the GPU in
`blk_begin`/synchronisation, ~28 % in oneMKL-on-GPU calls, 9 % in host MKL; the device time of the OpenMP kernels is only ~10 % of the wall time (2.9 s of 29.5 s for 2 sweeps, about 76 000 small launches per sweep).
Further gains need fewer, larger kernels (code work), not tuning.

Larger lattices (8 ranks per node, `N_delay` 128): L = 28: 16.0 s/sweep, 0.50 sweeps/s per node, 1.30 GB device memory per rank;
L = 32: 26.3 s/sweep, 0.30 sweeps/s per node, 2.49 GB per rank (a tile has 64 GB).

## Parameter scans (PARALLEL_PARAMS)
`Scripts_and_Parameters_files/make_PP_Uscan.py` writes `parameters`, `Temp_<i>/parameters`, `seeds` (and a filled `job.sh` from
`JobfileSuperMUC-NG-Phase2-GPU-PP.sh`, never submitted). Tested on 2 nodes: 4 U values x 4 ranks, L = 28 (every group GPU ACTIVE, output in `Temp_<i>/`).

### Sizing of the U scan (L = 28 and 32, 10 U values, 16 ranks per set)
10 sets x 16 ranks = 160 ranks = **20 nodes per lattice size** (40 nodes if both run together, limit of `general` is 96). One rank per tile keeps all tiles busy;
fewer nodes (more ranks per tile) only time-slice. A job of T hours on N nodes costs 4 N T GPU-h, so the 900 GPU-h of `pn93la` (225 node-h) are the limit,
not the node count: 20 nodes for 4 h = 320 GPU-h.
| L | s/sweep | 20 nodes, T h | GPU-h | sweeps per chain | sweeps per U (16 chains) |
|---|---|---|---|---|---|
| 28 | 16.0 | 4 | 320 | 900 | 14 400 |
| 32 | 26.3 | 6.5 | 520 | 890 | 14 200 |
The 840 GPU-h of both rows leave 60 GPU-h for tests. The L = 24 reference run has about 15 000 usable sweeps per U in total, but its chains ran about 4 800 sweeps each (1000 sweeps of
warm-up, MASTER report of 2026-11-10); with 900 sweeps per chain the warm-up would eat the whole run. More budget (the projects `pn73xu` and `pr94vu` also show 900 GPU-h),
fewer U values or fewer chains with longer runs, or warm starts from `confout` of a shorter run are the options.
