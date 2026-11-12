#!/bin/bash
#
# Jobscript for the GPU build of ALF (configure.sh IntelLLVM MPI GPU=intel HDF5) on SuperMUC-NG Phase 2 (LRZ).
# Each node has 4 Intel Data Center GPU Max 1550 (2 tiles each, 8 tiles) and 2 x 56 cores (112). Nodes are not shared.
#
# The following jobscript contains a few 'variables' marked by the ##...## pattern.
# The user has to provide the appropriate values, e.g. replace ##Nnodes## by 1.
#
# Measured layout (L = 24 pi-flux Hubbard with point defect, submatrix updates, see docs/gpu/SUPERMUC_NG2.md):
#  * ONE MPI rank (= one independent Markov chain) per GPU tile: NtaskPnode = 8, 14 cores (OpenMP threads for the host MKL) each.
#    More ranks than tiles do not help: PVC has no MPS, the ranks on one tile are time-sliced (aggregate throughput stays flat).
#  * N_delay = 128 in the parameters file (29 % faster per sweep than the usual 32 on this hardware).
#  * The device of a rank is rank mod 8 (ZE_FLAT_DEVICE_HIERARCHY=FLAT, every tile is a device), so the ranks of a node spread over
#    all 4 cards. With fewer than 8 ranks per node bind one rank per card: ZE_AFFINITY_MASK=<2*local rank>.
#  * The executable has to be built ahead of time for PVC (AOT, see docs/gpu/SUPERMUC_NG2.md); JIT costs about 1 minute per rank at start.
#
# Partitions on Phase 2: test (<= 1 h, <= 4 nodes), general (<= 2 d, <= 96 nodes), large, preempt (no micro).
# The project budget is in GPU-h: 1 node-hour = 4 GPU-h.
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
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=14

#Important
module load slurm_setup

# Environment the executable was built with (oneAPI 2025.1; the LRZ HDF5 module matches its ifx / Intel MPI).
# Newer oneAPI (2025.3) builds but fails at run time (OpenMP/SYCL interop), do not use it on this system.
module load intel-toolkit/2025.1.0 hdf5/1.14.5-intel25-impi

# Every Level Zero tile is a separate OpenMP/SYCL device
export ZE_FLAT_DEVICE_HIERARCHY=FLAT
# A missing device is an error instead of a silent run on the host
export OMP_TARGET_OFFLOAD=MANDATORY
# Host threads of every rank (host MKL of the submatrix flush)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_AFFINITY=granularity=fine,compact,1,0

# Pinning of the MPI ranks: every rank gets 14 cores next to its GPU tile
export I_MPI_PIN_DOMAIN=auto
export I_MPI_PIN_ORDER=compact

# Needs ##EXECUTABLE## = ALF.out built with GPU=intel (AOT) and a 'seeds' file with one integer per MPI rank
mpiexec -n $SLURM_NTASKS ##EXECUTABLE##
