#!/bin/bash
#SBATCH --job-name=tm_lbench
#SBATCH --output=tm_lbench.%j.log
#SBATCH --error=tm_lbench.%j.err
#SBATCH --time=00:20:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# oneMKL SYCL LAPACK per-call wall time on one tile (lapack_bench.cpp), N = 576 (L=24) and 1024 (L=32), 1 and 8 host threads
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
J=$SCRATCH/large_prompt/team/jobs; cd $J
icpx -fsycl -qmkl=parallel -O2 lapack_bench.cpp -o lapack_bench || exit 1
export ZE_AFFINITY_MASK=0
ONEAPI_DEVICE_SELECTOR=level_zero:gpu ./lapack_bench 576 3 2>&1 | tail -2   # V1 probe
export SYCL_UR_USE_LEVEL_ZERO_V2=1
for t in 1 8; do for n in 576 1024; do OMP_NUM_THREADS=$t MKL_NUM_THREADS=$t ./lapack_bench $n 20; done; done
