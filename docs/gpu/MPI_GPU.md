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
