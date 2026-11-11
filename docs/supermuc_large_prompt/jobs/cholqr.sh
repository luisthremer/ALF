#!/bin/bash
#SBATCH --job-name=lp_cqr
#SBATCH --output=lp_cqr.%j.log
#SBATCH --error=lp_cqr.%j.err
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
ifx -O2 -fiopenmp -fopenmp-targets=spir64 -qmkl -I$MKLROOT/include cholqr_bench.F90 -o cholqr_bench > build_cqr.log 2>&1; grep -i "error" build_cqr.log | head; 
export ZE_AFFINITY_MASK=0
for t in 1 8; do echo "== host threads $t"; OMP_NUM_THREADS=$t MKL_NUM_THREADS=$t ./cholqr_bench; done
