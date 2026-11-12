#!/bin/bash
#SBATCH --job-name=lp_la
#SBATCH --output=lp_la.%j.log
#SBATCH --error=lp_la.%j.err
#SBATCH --time=00:20:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# per-routine times of the alf_la ONEMKL backend on one PVC tile (repo benchmark testsuite/alf_la.tests/t_timing, N = 400 and 1024)
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
T=$SCRATCH/large_prompt/build_gpu-ser/testsuite/tests/alf_la.tests/t_timing
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1
echo "host $(hostname) $(date -Is)"; ls -l $T
for th in 1 8; do echo "=== OMP_NUM_THREADS=$th (host threads only matter for the host routines)"; OMP_NUM_THREADS=$th $T 2>&1 | grep -v "^ *$" | cut -c1-120; done
echo "done $(date -Is)"
