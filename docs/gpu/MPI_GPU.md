# MPI together with the GPU path (NVIDIA)

Status 2026-11-08, branch `gpu-delayed-submatrix-defect-production-old-Ntau-mpi`.

## What works

`. ./configure.sh PGI MPI GPU=nvidia HDF5 NO-INTERACTIVE` (HPC-X OpenMPI `mpifort`, nvfortran) builds and runs. Every MPI rank is an
independent Markov chain with its own device mirror and the host random numbers; the bins are averaged over the ranks of a group by
the unchanged host-side `MPI_REDUCE` (one `data.h5`, like the CPU runs). The device of a rank is `rank mod omp_get_num_devices()`,
so a node with several GPUs spreads the ranks over the devices. Tempering stays refused on the GPU path.

Run: `OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1 mpirun -n N --bind-to core ALF.out`.

Check: one rank of the MPI build reproduces the non-MPI run (same seeds, L = 16, first 3 bins) to 1e-12 (Ener_scal 7e-12, SpinZ_localtau
1.4e-12, SpinZ_local 7e-15).

## Rank-to-device binding (known limitation, 2026-11-13)

The device of a rank is `Mod(env_rank, omp_get_num_devices())` (`Prog/gpu_state_mod.F90:552`, called in `gpu_supported` before `la_init`), where
`env_rank` is the **global** `MPI_COMM_WORLD` rank (`env_mpi_ids`, `gpu_state_mod.F90:327`), not the node-local rank. The text above ("`rank mod
omp_get_num_devices()`, so a node with several GPUs spreads the ranks over the devices") is true only under these conditions:

- **Block placement**: consecutive ranks fill a node (Slurm `--distribution=block`, the default of `srun` for `--ntasks-per-node`), **and** the number of
  ranks per node is a multiple of the number of visible devices (or the nodes start at a multiple of it). Then global rank mod devices = local rank mod
  devices.
- **Or an explicit mask per rank**: a wrapper that gives every rank exactly one visible device (`ZE_AFFINITY_MASK=<tile>` on Intel, `CUDA_VISIBLE_DEVICES=<gpu>`
  on NVIDIA, see `scripts/helma/rank_wrap.sh`). Then `omp_get_num_devices() = 1` and the formula gives device 0, which is the masked one.

What goes wrong otherwise: with round-robin placement over several nodes (rank `r` on node `r mod n_nodes`), the ranks of one node are `k, k + n_nodes, ...`;
with 4 nodes and 4 devices per node every rank of a node gets the same device (0, 4, 8, ... all map to device 0 of the node) while the other devices stay idle.
With fewer ranks than devices per node the ranks land on the first devices only (SuperMUC: 4 ranks on 8 tiles use 2 cards, `sessions/supermuc_large_prompt/FINDINGS.md`
section 6, item 3; one rank per card with `ZE_AFFINITY_MASK=2*local_rank` was 1.8 times faster). The code neither checks nor warns. Not changed in the code; also listed in
`NOT_IMPLEMENTED.md` section 5.5. Check the device of every rank in `info` (`GPU path (gpu_state): ACTIVE ..., OpenMP device N`) and with the vendor tool of the machine.
The CUDA note in `alf_la_mod.F90` (bind ranks with `CUDA_VISIBLE_DEVICES`, the OpenACC current device must equal the OpenMP default device) is the reason the mask is preferred on NVIDIA.

## Several ranks on one A30 (L = 16, submatrix `N_delay = 32`, 30 sweeps, nividia, CUDA MPS daemon)

CPU time of rank 0 for 30 sweeps; one rank alone: 38.8 s (1.29 s per sweep).

| ranks | MPS thread share per rank | time per chain | vs 1 rank | aggregate sweeps/s |
|---|---|---|---|---|
| 1 | 100 % | 38.8 s | - | 0.73 |
| 1 | 25 % | 47.0 s | +21 % | 0.61 |
| 2 | 50 % | 45.3 s | +17 % | 1.25 |
| 4 | 25 % | 54.6 s | +41 % | 2.09 |
| 8 | 13 % | 81.4 s | +110 % | 2.83 |
| 16 | default | 145.6 s | +275 % | 3.18 |
| 16 | 50 % | | | 3.92 |
| 16 | 25 % | | | 4.68 |

MPS itself costs nothing (1 rank under MPS = 1 rank without). The A30 is shared, not multiplied: one chain already keeps it busy 58 % of
the time (`nvidia-smi dmon`), two chains 82 %. Acceptance (0.3375) and `Precision Green` (1.0e-14) do not depend on the number of ranks.

**Recommendation (owner decision 2026-11-08): one rank per A30**, so that the wall-clock time of a chain does not suffer. More ranks
only help the aggregate throughput (up to 4.7 sweeps/s with 16 ranks at 25 % thread share, each chain 3x slower). Without MPS the ranks
time-slice and are 13-15x slower (`CLAUDE.md`). Several GPUs per node scale per rank (not measured here).
