#!/bin/bash
#
# ALF job script for SuperMUC-NG Phase 2 (LRZ, Intel Data Center GPU Max 1550 "PVC") with the GPU build of ALF
# (configure.sh IntelLLVM MPI GPU=intel HDF5, see docs/gpu/BUILD.md and docs/supermuc_large_prompt/FINDINGS.md).
# For the CPU-only SuperMUC-NG Phase 1 nodes use JobfileSuperMUC-NG.sh.
#
# Variables marked ##...## have to be replaced, e.g. ##Nranks## by 4.
#
# Node layout (confirmed): 2 x Xeon Platinum 8480L (112 cores, 224 hardware threads), 512 GB, 4 GPU cards = 8 tiles.
# Slurm has no GPU gres on this system: a job always gets whole nodes and sees all 4 cards (8 devices with
# ZE_FLAT_DEVICE_HIERARCHY=FLAT, the default). Every MPI rank is one independent Markov chain ("walker") with its
# own device mirror; ranks are bound to GPUs below.
#
#SBATCH --job-name ##NAME##
#SBATCH --output=out.%j.log
#SBATCH --error=err.%j.log
#Notification and type
#SBATCH --mail-type=ALL
#SBATCH --mail-user=##EMAIL##
# Wall clock limit (HH:MM:SS); partition test: max 1 h, 4 nodes; general: max 2 d, 96 nodes
#SBATCH --time=##TIME##
#SBATCH --no-requeue
#Setup of execution environment
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=##projectID##

#available partitions: test, general, large, preempt (check with sinfo)
#SBATCH --partition=##PARTITION##
#SBATCH --nodes=##Nnodes##
# ranks per node: 4 = one walker per GPU card (recommended, measured); up to 8 = one per tile (not measured)
#SBATCH --ntasks-per-node=##NtaskPnode##
#SBATCH --cpus-per-task=1

#Important
module load slurm_setup
# ifx >= 2025.3 is needed (omp_lib interop, see FINDINGS.md); compiler, MKL and HDF5 versions must match the build
module load intel/2025.3.0 intel-mpi/2021.15.0
module unload intel-mkl
module load intel-mkl/2025.3.0 hdf5/1.14.5-intel25

# GPU path: a missing device must stop the run instead of silently running on the CPU
export OMP_TARGET_OFFLOAD=MANDATORY
# Level Zero staging buffer 1 MB (default is smaller): the 32 KB D->H block copy of the submatrix update drops from ~380 us
# to ~22 us; one chain L=24 10.5 -> 7.05 s/sweep, L=32 31.5 -> 23.4, bitwise-identical results (team session 2026-11-12).
export LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576
# one host thread per walker (all heavy work is on the GPU)
export OMP_NUM_THREADS=1
# Rank -> GPU. omp_set_default_device(rank mod #visible devices) is done by ALF itself; to put ranks on different CARDS
# (rank mod 8 would put ranks 0..3 on cards 0 and 1 only) every rank sees exactly one tile:
# rank r -> device (r mod 4)*2 + r/4  (4 ranks: tiles 0,2,4,6 = one per card; 8 ranks: all tiles)
BIND='r=${MPI_LOCALRANKID:-$PMI_RANK}; export ZE_AFFINITY_MASK=$(( (r % 4) * 2 + r / 4 )); exec "$@"'

bash ./out_to_in.sh
mpiexec -n $SLURM_NTASKS bash -c "$BIND" _ ##EXECUTABLE##
