#!/bin/bash
#SBATCH --job-name=lp_dqr
#SBATCH --output=lp_dqr.%j.log
#SBATCH --error=lp_dqr.%j.err
#SBATCH --time=00:20:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
D=$SCRATCH/large_prompt/phase2; cd $D
ifx -O2 -fiopenmp -fopenmp-targets=spir64 -qmkl -I$MKLROOT/include devqr_bench.F90 -o devqr_bench > build_dqr.log 2>&1; grep -i -B2 -A3 "error" build_dqr.log | head -40
export ZE_AFFINITY_MASK=0
for t in 1 8; do echo "== host threads $t"; OMP_NUM_THREADS=$t MKL_NUM_THREADS=$t ./devqr_bench 2>&1 | head -60; done
