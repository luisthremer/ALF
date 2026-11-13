#!/bin/bash
#SBATCH --job-name=p8b_gp
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:10:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
cd $SCRATCH/large_prompt/team/p8/blk/jobs
echo "start $(date -Is) $(hostname)"
icpx -fsycl -O2 -std=c++17 -fsycl-targets=spir64_gen -Xs "-device pvc" -I$MKLROOT/include graph_probe.cpp -o graph_probe -qmkl -lmkl_sycl_blas 2>&1 | tail -20
ZE_AFFINITY_MASK=0 mpiexec -n 1 ./graph_probe 2>&1 | tail -20
echo "end $(date -Is)"
