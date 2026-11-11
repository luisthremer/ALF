#!/bin/bash
#SBATCH --job-name=lp_lat
#SBATCH --output=lp_lat.%j.log
#SBATCH --error=lp_lat.%j.err
#SBATCH --time=00:15:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
D=$SCRATCH/large_prompt/phase2; mkdir -p $D; cd $D
ifx -O2 -fiopenmp -fopenmp-targets=spir64 -qmkl -I$MKLROOT/include latbench.F90 -o latbench > build_lat.log 2>&1; tail -5 build_lat.log
export OMP_NUM_THREADS=1 ZE_AFFINITY_MASK=0
echo "== default (1 tile, 1 process)"; ./latbench
echo "== CCS 4"; ZEX_NUMBER_OF_CCS=0:4 ./latbench
echo "== command mode async_ordered (if known)"; LIBOMPTARGET_LEVEL_ZERO_COMMAND_MODE=async_ordered ./latbench
echo "== plugin env"; env | grep -i LIBOMPTARGET
