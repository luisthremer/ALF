#!/bin/bash
#
# ALF production job on SuperMUC-NG Phase 2 (Intel Max 1550 GPUs): PARALLEL_PARAMS build, one parameter set per U value with
# ##RanksPerSet## MPI walkers each, 4 walkers per GPU tile. For a single parameter set use JobfileSuperMUC-NG-Phase2-GPU.sh.
#
# Setup (login node, python3 only):   python3 Start_PP_GPU/make_pp_dirs.py RUNDIR ##L## ##ULIST## ##RanksPerSet## ##NSweep## 100000 ##NDelay## ##CPU_MAX_H##
#   e.g. L=28, U list 5.0,5.1,5.2,5.25,5.3,5.35,5.4,5.45,5.5,5.6, 16 ranks per set, NSweep 50, CPU_MAX 46.0 (hours, stops cleanly before the limit)
# Submit from RUNDIR:                 sbatch JobfileSuperMUC-NG-Phase2-GPU-PP.sh   (copy this file and fill the ##...## first)
# Restart after the time limit:       just sbatch again; the script moves confout_* to confin_* (out_to_in_pp.sh) and data.h5 gets appended.
#
# Ranks = (number of U values) x ##RanksPerSet##; nodes = ranks / 32 (32 ranks = 4 per tile x 8 tiles per node).
#
#SBATCH --job-name ##NAME##
#SBATCH --output=out.%j.log
#SBATCH --error=err.%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=##EMAIL##
# partition general: max 96 nodes, 2 days
#SBATCH --time=##TIME##
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=##projectID##
#SBATCH --partition=##PARTITION##
#SBATCH --nodes=##Nnodes##
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1

module load slurm_setup
module load intel/2025.3.0 intel-mpi/2021.15.0
module unload intel-mkl
module load intel-mkl/2025.3.0 hdf5/1.14.5-intel25

export OMP_TARGET_OFFLOAD=MANDATORY
export OMP_NUM_THREADS=##OMPT##
# expose 4 compute engines per tile: 4 walkers per tile then run concurrently (measured 3.5x node throughput; without it they time-slice)
export ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4
# no I_MPI_PIN_DOMAIN: pinning every rank to a single cpu slows the GPU runtime's helper threads by ~35 % (measured)
# device of rank r = r mod 8 (ALF), ranks are placed node by node (32 per node), so every tile gets 4 walkers

bash ./out_to_in_pp.sh 2>/dev/null || bash ./Start_PP_GPU/out_to_in_pp.sh
mpiexec -n $SLURM_NTASKS ##EXECUTABLE##
