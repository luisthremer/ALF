# ALF GPU on Helma (NHR@FAU, NVIDIA H100): build, benchmarks, L = 32 U-scan production

Status 2026-11-11 (written 15:20 CET). Everything below marked "measured" comes from the runs listed in section 5; "projection" is arithmetic on
those numbers, not a measurement. The production chains were submitted on 2026-11-11 at about 14:50; **at the time of writing no production
link had started** (section 6 gives the current state, section 8 the open points).

## 1. Machine and project

- Helma, login `ssh helma` (ProxyJump `csnhr`), Slurm 25.11, project/account `b133ae` (test/porting contingent 3000 GPU-h).
- GPU node `h100` partition: 4 x H100 (95.8 GB, cc 9.0), 128 cores, driver 615.71.09 (CUDA UMD 13.4). Per GPU 32 cores, ~175 GB host RAM.
- Partitions used: `preempt` (default, 1 node, max 48 h, `PreemptMode=CANCEL`, `GraceTime=900`: a preempted job is cancelled, **not** requeued)
  and `h100` (24 h). `h200` is not allowed.
- Workspace (GPU nodes mount only `$HOME`, `/hnvme`, `$TMPDIR`): `/hnvme/workspace/b133ae29-ALF_GPU` (90 days).

## 2. Build recipe (measured)

MACHINE `HELMAGPU` in `configure.sh` (login shell, `module` only exists there):

```bash
export ALF_HDF5_DIR=/hnvme/workspace/b133ae29-ALF_GPU/hdf5_257      # ALF builds its own HDF5 1.14.6 here
. ./configure.sh HELMAGPU PP GPU=nvidia HDF5 NO-INTERACTIVE
make lib program                                                    # serial
```

`HELMAGPU` loads `gpu-env/2025`, `nvhpc/25.7`, `openmpi/5.0.8-nvhpc25.7-cuda` (pulls `cuda/13.0.2`), and compiles with `-mp=gpu -gpu=cc90,cuda12.9`
(`ALF_GPU_CC` overrides). Other pairing: `ALF_HELMA_NVHPC=24.11` (`nvhpc/24.11` + `openmpi/5.0.5-nvhpc24.11-cuda`, flags `cc90,cuda12.6`).

Findings while building (all measured 2026-11-11):

| Item | Result |
|---|---|
| `nvhpc/24.11` GPU build | compiles and links (needs `cuda12.6`: the default CUDA 11.8 libs lack `cusolverDnSetDeterministicMode`), but the **L = 8 GPU run aborts**: `alf_la: la_qrpp: NaN or Inf in the matrix`. The same-compiler CPU run is fine. Not used. |
| `nvhpc/25.7` GPU build | L = 8 run correct (section 3). Without `cuda12.9` the default (CUDA 13) device compile fails (`nvnvvmd: open ... nvhpc_cuda_real2_runtime.bc failed`). |
| HDF5 1.14.6 with `nvc` | link error `__truncdfhf2` (Float16): `HDF5/install_hdf5.sh` now passes `--disable-nonstandard-feature-float16` for `nvc`. No HDF5 module exists for nvhpc 25.7, so ALF's own installer is used. |
| Open MPI under Slurm | `mpirun` sees one slot per task: use `--oversubscribe` (the jobs have 32 cores per GPU; `--bind-to core --map-by core` is used and `--report-bindings` shows one core per rank). The `/dev/knem` warning is silenced with `--mca smsc ^knem`. |
| GPU nodes | no `/usr/bin/time`, system python has no numpy/h5py (comparison scripts use `/hnvme/workspace/b133ae29-CLAUDE/venv_h5/bin/python`, read-only). |

## 3. Correctness (measured)

L = 8, `N_FL = 2`, `h0 = 1`, `Mz = .F.`, Checkerboard, U = 5.0, `N_delay` 64, 4 bins x 10 sweeps, 1 rank, same seed; GPU `ALF-gpu-nvhpc25.7-cc90.out`
vs CPU `ALF-cpu-nvhpc25.7.out` (`GPU=none`, same nvfortran 25.7; `docs/gpu/scripts/helma/smoke.sh`, job 1302917):

- worst relative deviation over all datasets of `data.h5` and the final `confout`: **5.2e-13** (norm: max|a-b| / max|a| per dataset, `cmp_h5.py`)
- acceptance 0.3377294921875 identical; GPU Precision Green mean 8.4e-15, max 2.4e-12; Precision Phase max 5.2e-13
- `info`: `GPU path (gpu_state): ACTIVE`, `GPU observables: device equal-time measurement ACTIVE`.

Restart (production = chain of <= 24 h links; `restart_test.sh`, job 1302960): ALF reads `confin_*`, not `confout_*`, so `run_pp.sh` renames
`confout_* -> confin_*` in every `Temp_i` before each link (what `out_to_in.sh` does) and keeps the previous `data.h5`/`confin` as `.bak`/`.prev`.
`NBin` counts the bins of the *current* invocation, `data.h5` is appended. L = 8, 2 groups x 2 walkers: 3 bins + restart + 6 bins = 9 bins; the
first 6 bins of the restarted run agree with an uninterrupted 6-bin run (max |diff| 7e-15 for `Ener_scal`, 0 for `Green_eq`, `SpinZ_local`),
so the chain (including the RNG state) is continued, not restarted. Without the rename ALF starts a fresh configuration.

## 4. Benchmarks (measured)

Settings: parameters of `/jakob_params` (Hubbard, `h0 = 1`, `N_FL = 2`, `Mz = .F.`, Beta = L, Dtau 0.1, Nwrap 10, submatrix), only L, `N_Phi = L^2/2`,
N_delay and U differ; `NSweep = 2` x `NBin = 3` = **6 sweeps** per test, one GPU (H100, node in `preempt`), all ranks of one test on that GPU
(PARALLEL_PARAMS binary, one group, `mpi_per_parameter_set` = ranks per GPU), `OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1`, `mpirun --bind-to core`,
one CUDA MPS daemon, "MPS" = `CUDA_MPS_ACTIVE_THREAD_PERCENTAGE`. s/sweep per chain = (FIN TIME - START TIME in `info`) / 6, taken from rank 0 of
the group (excludes the ~0.3-0.6 s mirror setup and the MPI start; includes the first-sweep warm-up of the device kernels, so it is slightly
pessimistic). Sweeps/s per GPU = ranks / (s/sweep). GPU memory = max `fb` of `nvidia-smi dmon -s um` (all ranks). SM % = mean `sm` over dmon samples
with sm > 0 (10 s interval, so coarse). Raw: `docs/gpu/scripts/helma/results/bench1_results.csv`.

### 4.1 L = 24 (like-for-like with the old A30 numbers), `N_delay` 32

| ranks / GPU | MPS | H100 s/sweep per chain | H100 sweeps/s per GPU | GPU mem MiB | Precision Green max |
|---|---|---|---|---|---|
| 1 | off | 3.51 | 0.285 | 1656 | 2.1e-11 |
| 2 | off | 29.5 | 0.068 | 3308 | 4.5e-11 |
| 2 | 25 % | 4.54 | 0.440 | 2628 | 5.0e-11 |
| 4 | 25 % | 4.98 | 0.804 | 5187 | 5.0e-11 |
| 8 | 25 % | 6.80 | 1.176 | 10305 | 1.3e-10 |

Old A30 numbers (read from `nividia:/home/luis/alf-gpu/work/prodmpi/runs/{A_mpi1,B_n1,B_n2,B_n4,B_n8,PROD_L24_U5.3_h1_4w}`, nothing was run there; all are
L = 24, Beta = 24, `N_FL = 2`, `h0 = 1`, `Mz = .F.`, `N_delay` 32, 1 x A30, MPI ranks on one GPU; the benchmark runs have 20-30 sweeps, `NSweep` 10; the
production has `NSweep` 50). They are `CPU Time` of rank 0 in `info` divided by the number of sweeps, i.e. **including start-up and a different
timer than the H100 numbers**; whether MPS was active in the B runs is not recorded in the run directories (the PP_GPU.md measurements used MPS 25 %):

| A30 ranks | CPU time (s) / sweeps | s/sweep per chain | sweeps/s per GPU |
|---|---|---|---|
| 1 (A_mpi1) | 110.4 / 20 | 5.5 | 0.18 |
| 1 (B_n1) | 166.2 / 30 | 5.5 | 0.18 |
| 2 (B_n2) | 210.2 / 30 | 7.0 | 0.29 |
| 4 (B_n4) | 317.5 / 30 | 10.6 | 0.38 |
| 8 (B_n8) | 547.0 / 30 | 18.2 | 0.44 |
| 4 (PROD_L24_U5.3_h1_4w, 14.17 h CPU time, 97 bins x 50 sweeps) | 51012.8 / 4850 | 10.5 | 0.38 |

At 8 ranks per GPU the H100 does 1.18 sweeps/s against 0.44 (A30): **2.7x** (single rank 1.5-1.6x). The A30 "single chain" figure of TUNING_DEFECT.md for L = 24
(5.25 s/sweep, `N_delay` 32, nvfortran 26.9, measured separately on 2026-11-07) agrees with the 5.5 s above.

### 4.2 L = 32, U = 5.3 (1 GPU = one group), `N_delay` 64 unless stated

| ranks / GPU | MPS | s/sweep per chain | sweeps/s per GPU | GPU mem MiB | SM % | Precision Green max |
|---|---|---|---|---|---|---|
| 1 | off | 9.82 | 0.102 | 3584 | 55 | 1.8e-10 |
| 1 | 25 % | 14.81 | 0.068 | 3270 | 80 | 1.4e-10 |
| 2 | off | 42.7 | 0.047 | 7164 | 93 | 1.8e-10 |
| 2 | 25 % | 16.26 | 0.123 | 6472 | 81 | 1.4e-10 |
| 4 | off | > 420 s for 6 sweeps (timeout, not finished) | - | 14325 | 97 | - |
| 4 | 25 % | 19.98 | 0.200 | 12875 | 91 | 6.0e-10 |
| 4 | 50 % | 18.32 | 0.218 | 13395 | 98 | 7.2e-10 |
| 8 | off | > 420 s for 6 sweeps (timeout, not finished) | - | 28647 | 99 | - |
| 8 | 25 % | 34.49 | 0.232 | 25683 | 94 | 6.0e-10 |
| 8 | 50 % | 33.91 | 0.236 | 26723 | 100 | 7.2e-10 |
| 16 | 25 % | 63.99 | 0.250 | 51294 | 100 | 6.0e-10 |
| 16 | 50 % | 62.96 | 0.254 | 53374 | 100 | 7.2e-10 |

`N_delay` at 8 ranks / GPU, MPS 25 %: 32 -> 66.5 s/sweep (0.120 sweeps/s), 64 -> 34.5 (0.232), 128 -> 34.3 (0.233). So 64 (as in TUNING_DEFECT.md for L >= 32)
is the choice; 128 gives nothing. Old A30 L = 32 number (TUNING_DEFECT.md, 2026-11-07, `N_FL = 2`, `Mz = .F.`, `N_delay` 64, one rank on the GPU): 18.2 s/sweep,
i.e. 3.0 min per 10-sweep bin, 3077 MiB. H100, one rank without MPS: 9.8 s/sweep (**1.85x**), 3584 MiB. No multi-rank L = 32 numbers exist for the A30.

Reading of the numbers:
- **MPS is mandatory.** Several ranks on one GPU without MPS time-slice: 2 ranks are 6.5x (L = 24: 29.5 vs 4.54 s/sweep) and 2.6x (L = 32: 42.7 vs 16.3 s/sweep) slower per chain than with MPS,
  and 4 or 8 ranks at L = 32 did not finish 6 sweeps in 420 s (with MPS 6 sweeps take 122 s and 210 s of wall time).
- A single rank cannot fill an H100 (SM 55 %); aggregate throughput rises from 0.10 (1 rank, no MPS) to 0.23 (8 ranks) and 0.25 sweeps/s (16 ranks), +6 % from 8 to 16.
- MPS thread share 25 % and 50 % are indistinguishable at 8 and 16 ranks (within 2 %); 25 % was used for production (as in `PP_GPU.md`).
- Memory per rank ~3.2 GB (L = 32): 16 ranks use 51-53 GB of the 95.8 GB.
- Precision Green max 1.4e-10 ... 7.2e-10 for all L = 32 runs (requirement < 1e-8).
- Not measured: the full 4-GPU node (4 U on 4 GPUs of one node), 5-6 sweeps are short (the first-sweep warm-up and one outlier would show),
  GPU utilisation at finer time resolution. The 4-GPU test was dropped from the first-day set because production jobs are single-GPU (section 7); the rank-to-GPU
  binding (`rank_wrap.sh`: rank r sees only GPU r / RPG) is implemented and was used on 2 GPUs only in `restart_test.sh` (2 groups, correct restart on both).

## 5. Runs behind the numbers

| Job | What |
|---|---|
| 1302917 | L = 8 GPU vs same-compiler CPU (nvhpc 25.7) |
| 1302913 | L = 8 GPU (25.7) vs CPU built with 24.11: 5.2e-13 |
| 1302857 | L = 8 with the nvhpc 24.11 GPU build: aborts (NaN in `la_qrpp`) |
| 1302960 | restart test (2 GPUs, 2 x 2 ranks) |
| 1302959 | 1-GPU benchmark set (L = 24 and 32 tables) |
| 1302937 / 1302916 | earlier benchmark attempts (cancelled: `mpirun` slot error, then no per-config timeout) |

All in `preempt`. GPU-h used for build tests and benchmarks: well below 10 GPU-h (benchmark job 1302959 ran 1 h on 1 GPU; other jobs minutes).

## 6. Queue findings (2026-11-11, `sbatch --test-only`, 24 h unless noted)

| Request | Estimated start |
|---|---|
| `h100`, 1 GPU | 2026-11-13 14:58 (later 13:27) |
| `h100`, 1 node (4 GPU) | 2026-11-13 15:06 |
| `h100`, 2 nodes | 2026-11-13 15:18 |
| `h100`, 4 nodes | 2026-11-13 16:59 |
| `h100`, 8 nodes | 2026-11-13 19:15 |
| `preempt`, 1 GPU / 2 GPU / 1 node, 24 h or 48 h | 2026-11-11 17:13-17:21 (later estimates: 19:31) |

On 2026-11-11 `h100` had 88 of 93 nodes allocated, `preempt` 178-181 allocated. Early `preempt` jobs (1 GPU, 5-30 min) started within seconds
(idle nodes at that moment); the 20 production jobs submitted at about 14:50 were all still pending ("Priority") at 15:20.
The estimates are for jobs that are alone in the queue; `preempt` has the earlier start because the partition spans more nodes, but its jobs can be
cancelled after 2 h when `h100`/higher-priority jobs need the node.

Chosen layout: **20 independent single-GPU job chains in `preempt`** (smallest allocation backfills best; a preempted link only loses the bin in progress of its 8 walkers; every link
restarts from `confout`). A multi-node PP job was rejected: `preempt` allows 1 node, and `h100` multi-node jobs start days later.

## 7. Production layout (as submitted)

Parameters: Jakob's reference file (`docs/gpu/scripts/helma/parameters_ref_L24_U5.3_h1_4w`, read from nividia `PROD_L24_U5.3_h1_4w/parameters`) with
L1 = L2 = 32, Beta = 32.d0, N_Phi = 512, **NSweep = 10**, N_delay = 64, ham_U per run, CPU_MAX set per link. `Mz = .F.`, `N_FL = 2`, `ham_h0 = 1.d0`, `Dtau = 0.1d0`, `Nwrap = 10`,
`Projector = .F.`, Checkerboard/Symm .T. as in the reference. `diff` of every `Temp_0/parameters` against the reference shows only
L1, L2, N_Phi, Beta, NSweep, CPU_MAX, N_delay, ham_U (checked on the generated directories). `setup_pp.py` writes `Temp_i/parameters` directly (no f90nml).

- U = 5.0, 5.1, 5.2, 5.25, 5.3, 5.35, 5.4, 5.45, 5.5, 5.6; 2 parts per U, 8 walkers (MPI ranks) per part, **16 walkers per U**, 160 ranks in total.
- Each part = 1 job chain, 1 GPU, 8 ranks, MPS thread share 25 %, `OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1`, `--bind-to core`; 4 links of 24 h in `preempt`,
  each link depends on the previous (`--dependency=afterany`); `CPU_MAX` of a link = remaining walltime - 0.3 h (a bin of 8 walkers takes 10 x 34 s = 5.7 min).
- Seeds: one integer per rank, `sha256(master|L|U|walker)`, all 160 distinct (checked, `seed_check.txt`).
- Runs: `/hnvme/workspace/b133ae29-ALF_GPU/runs/PROD_L32_h1_Uscan/L32_U<U>_p<0|1>/Temp_0/{data.h5,confout_*,info}`; job ids in `README_jobs.txt` there;
  Slurm logs `/hnvme/workspace/b133ae29-ALF_GPU/slurm/`.
- Budget: worst case 20 chains x 4 links x 24 h x 1 GPU = **1920 GPU-h** (limit 2500 incl. all chain links; 500 GPU-h reserve of the 3000 kept).

Projection (arithmetic on section 4.2, not measured in production): 0.232 sweeps/s per GPU at 8 ranks -> 20 000 sweeps per GPU per day; per U (2 GPUs) 40 000 sweeps/day
= 4 000 bins/day. Warm-up of 1000 sweeps per walker (100 bins) takes 1000 x 34.5 s = 9.6 h per chain and costs 16 000 sweeps per U (0.4 days of the two GPUs).
Four full links per U: ~160 000 sweeps (16 000 bins) of which 16 000 sweeps are warm-up.

### 7a. Change at about 17:00 (priority)

The 80 `preempt` jobs (none had started in ~2 h) were cancelled and the 20 chains resubmitted to **`h100`**: `sprio` priority 165 -> 1157
(partition term 1 -> 1000; fair-share 157 and job size 1 unchanged, QOS 0). `h100` jobs also cannot be cancelled by other jobs, so a link runs its full 24 h.
`sbatch --test-only` start estimate for 1 GPU in `h100`: 2026-11-13 17:19. Each running link submits the next one (`afterany` on itself, `prod_job.sh`), so only
20 jobs are ever pending (Slurm `bf_max_job_user=20`). Priority formula (`scontrol show config`): multifactor, weights QOS 100000, fair-share 10000, partition 10000,
age 1000, job size 1000, `PriorityFavorSmall=no`; `h100` PriorityTier 100 / `preempt` 1.

## 8. Open issues

- **Production has not started** at the time of writing (queue full). Start times, bins done per U and GPU-h are to be added from `squeue`/`sacct`/`info`.
- `CPU_MAX` limits CPU time, the link ends at the first bin after it; if CPU time per rank is below wall time the link may end later than planned and hit the walltime (cost: one bin).
- `preempt` can cancel a link after 2 h; the chain continues with the next link when it is scheduled (queue wait without GPU), the restart check above covers the file handling.
  A kill during an HDF5 write could corrupt `data.h5`: `run_pp.sh` validates it with `h5ls` at the start of each link and restores `data.h5.bak` (copy of the last readable file).
- 4-GPU node layout and finer GPU-utilisation measurements are open; 16 ranks per GPU would give +6 % sweeps/s at 2x the warm-up latency (not used).
- A30 comparison numbers are from old files with a different timer (see 4.1); the old PP_GPU.md numbers are `N_FL = 1` pi-flux and are not used here.

## 9. Files

`docs/gpu/scripts/helma/`: `setup_pp.py` (PP directory, seeds), `run_pp.sh` (MPS per GPU, restart rename, backup, `mpirun`), `rank_wrap.sh` (GPU binding),
`prod_job.sh` / `chain_submit.sh` / `submit_prod.sh` (chain), `bench_job.sh` / `bench_summarize.py` (benchmark), `smoke.sh`, `cmp_h5.py`, `restart_test.sh`,
`parameters_ref_L24_U5.3_h1_4w`. `configure.sh` (MACHINE `HELMAGPU`), `HDF5/install_hdf5.sh`.
