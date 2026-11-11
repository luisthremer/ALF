#!/bin/bash
#SBATCH --job-name=tm_c1
#SBATCH --output=tm_c1.%j.log
#SBATCH --error=tm_c1.%j.err
#SBATCH --time=00:20:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# Phase C quick probes: (1) per-launch cost vs array passing (descbench.F90), (2) oneMKL batched geqrf/getrf vs plain (lapack_bench).
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
J=$SCRATCH/large_prompt/team/jobs; cd $J
AOT="-fopenmp-targets=spir64_gen -Xopenmp-target-backend -device -Xopenmp-target-backend pvc"
export ZE_AFFINITY_MASK=0 OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=8
ifx -O2 -fiopenmp $AOT descbench.F90 -o descbench > b_desc.log 2>&1 || tail -20 b_desc.log
./descbench
LIBOMPTARGET_DEBUG=1 ./descbench 2>&1 | grep -c "first-private array"
icpx -fsycl -qmkl=parallel -O2 lapack_bench.cpp -o lapack_bench 2>&1 | tail -5
export SYCL_UR_USE_LEVEL_ZERO_V2=1
for n in 576 1024; do MKL_NUM_THREADS=8 ./lapack_bench $n 10 | grep -E "geqrf|getrf"; done
echo "done $(date -Is)"
