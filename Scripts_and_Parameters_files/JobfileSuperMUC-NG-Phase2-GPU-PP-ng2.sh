#!/bin/bash
#
# Jobscript for a parameter scan (PARALLEL_PARAMS, one parameter set per U value, several independent walkers per set) with the GPU
# build of ALF on SuperMUC-NG Phase 2 (LRZ):  configure.sh IntelLLVM PP GPU=intel HDF5 NO-INTERACTIVE   (+ AOT, see docs/gpu/platforms/intel_pvc_supermuc_ng2.md)
#
# Set up the run directory (parameters, Temp_<i>/parameters, seeds) with make_PP_Uscan.py in this directory, e.g. for L = 28:
#   python3 make_PP_Uscan.py runL28 --L 28 --U 5.0 5.1 5.2 5.25 5.3 5.35 5.4 5.45 5.5 5.6 --ranks-per-set 16 --cpu-max ##CPUMAX##
# then submit THIS script from inside that directory (the job runs in the submission directory):  cd runL28; sbatch ../JobfileSuperMUC-NG-Phase2-GPU-PP.sh
#
# The following jobscript contains a few 'variables' marked by the ##...## pattern.
#  NtaskPnode = ranks per node: 8 (one per GPU tile, 14 threads), 16 or 32 (2 or 4 chains per tile, 7 or 3 threads; cheaper in GPU-h but each chain is
#  2-3x slower, see docs/gpu/platforms/intel_pvc_supermuc_ng2.md).  Nthreads = 112 / NtaskPnode.  Nnodes = (number of U values * ranks per set) / NtaskPnode, e.g. 10 * 16 / 8 = 20.
#  make_PP_Uscan.py --ranks-per-node N --jobfile ... fills NtaskPnode, Nthreads and Nnodes in.    CPUMAX (hours) must be smaller than TIME (the run stops
#  after writing its last bin; leave 10 to 15 minutes for the start-up and the final bin).
#
# Layout (measured, docs/gpu/platforms/intel_pvc_supermuc_ng2.md): ONE MPI rank per GPU tile (8 ranks, 14 cores each) gives the shortest wall time; 2 or 4 ranks per tile
# (no MPS on PVC, time slicing) raise the node throughput by 0-27 % but slow every chain down 2-3x. Use N_delay = 128. One node costs 4 GPU-h per hour.
#
# Time per sweep and chain on a full node (8 ranks), N_delay = 128:   L = 24: 8.2 s,  L = 28: 16.0 s,  L = 32: 26.3 s.
#
#SBATCH --job-name ##NAME##
#SBATCH --output=out.%j.log
#SBATCH --error=err.%j.log
#Notification and type
#SBATCH --mail-type=ALL
#SBATCH --mail-user=##EMAIL##
# Wall clock limit (HH:MM:SS), at most 48:00:00 on 'general':
#SBATCH --time=##TIME##
#SBATCH --no-requeue
#Setup of execution environment
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=##projectID##

#available partitions: test (<= 4 nodes, 1 h), general (<= 96 nodes, 2 d)
#SBATCH --partition=##PARTITION##
#SBATCH --nodes=##Nnodes##
#SBATCH --ntasks-per-node=##NtaskPnode##
#SBATCH --cpus-per-task=##Nthreads##

#Important
module load slurm_setup
module load intel-toolkit/2025.1.0 hdf5/1.14.5-intel25-impi

export ZE_FLAT_DEVICE_HIERARCHY=FLAT
export OMP_TARGET_OFFLOAD=MANDATORY
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_AFFINITY=granularity=fine,compact,1,0
export I_MPI_PIN_DOMAIN=auto
export I_MPI_PIN_ORDER=compact

# Consistency of the setup: number of ranks = number of Temp_<i> directories * mpi_per_parameter_set
NSETS=$(ls -d Temp_* 2>/dev/null | wc -l)
RPS=$(grep -i mpi_per_parameter_set parameters | head -1 | sed 's/.*=[ ]*//')
if [ -z "$RPS" ] || [ $((NSETS * RPS)) -ne "$SLURM_NTASKS" ]; then
  echo "ERROR: $NSETS Temp_* directories x mpi_per_parameter_set=$RPS does not match the $SLURM_NTASKS MPI ranks" 1>&2
  exit 1
fi
NSEEDS=$(wc -l < seeds)
if [ "$NSEEDS" -lt "$SLURM_NTASKS" ]; then echo "ERROR: seeds has $NSEEDS lines, $SLURM_NTASKS needed" 1>&2; exit 1; fi

# Needs ##EXECUTABLE## = ALF.out built with PP and GPU=intel (AOT)
mpiexec -n $SLURM_NTASKS ##EXECUTABLE##
