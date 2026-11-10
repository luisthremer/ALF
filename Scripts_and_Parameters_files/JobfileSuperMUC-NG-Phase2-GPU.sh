#!/bin/bash
#
# Jobscript for the GPU build of ALF (configure.sh IntelLLVM MPI GPU=intel HDF5) on SuperMUC-NG Phase 2 (LRZ).
# Each node has 4 Intel Data Center GPU Max 1550 (2 tiles each, 8 tiles) and 2 x 56 cores. Nodes are not shared.
#
# The following jobscript contains a few 'variables' marked by the ##...## pattern.
# The user has to provide the appropriate values, e.g. replace ##Nnodes## by 1.
#
# ALF GPU path on Phase 2: every MPI rank is an independent Markov chain with its own device mirror, the device of a rank is
# rank mod omp_get_num_devices(). With ZE_FLAT_DEVICE_HIERARCHY=FLAT every tile is a device (8 per node), so use
# NtaskPnode <= 8 (one rank per tile, the default) and Nthreads = 1.
#
# Partitions on Phase 2: test (<= 1 h, <= 4 nodes), general (<= 2 d, <= 96 nodes), large, preempt (no micro).
# Phase 2 is charged in GPU-h (the project budget is in GPU-h).
#
#SBATCH --job-name ##NAME##
#SBATCH --output=out.%j.log
#SBATCH --error=err.%j.log
#Notification and type
#SBATCH --mail-type=ALL
#SBATCH --mail-user=##EMAIL##
# Wall clock limit (HH:MM:SS):
#SBATCH --time=##TIME##
#SBATCH --no-requeue
#Setup of execution environment
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=##projectID##

#available partitions: test, general, large, preempt
#SBATCH --partition=##PARTITION##
#SBATCH --nodes=##Nnodes##
#SBATCH --ntasks-per-node=##NtaskPnode##
#SBATCH --cpus-per-task=1

#Important
module load slurm_setup

# Environment the executable was built with (oneAPI 2025.1; the HDF5 module matches its ifx / Intel MPI)
module load intel-toolkit/2025.1.0 hdf5/1.14.5-intel25-impi

# Every Level Zero tile is a separate OpenMP/SYCL device
export ZE_FLAT_DEVICE_HIERARCHY=FLAT
# A missing device is an error instead of a silent run on the host
export OMP_TARGET_OFFLOAD=MANDATORY
export OMP_NUM_THREADS=1

# Pinning of the MPI ranks: ##NtaskPnode## ranks, each next to its GPU tile
export I_MPI_PIN_DOMAIN=auto
export I_MPI_PIN_CELL=core

# Needs ##EXECUTABLE## = ALF.out built with GPU=intel and a 'seeds' file with one integer per MPI rank
mpiexec -n $SLURM_NTASKS ##EXECUTABLE##
