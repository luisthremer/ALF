#!/bin/bash
#SBATCH --job-name=p6_ceil2
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH -x i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:30:00
# Gate 1b: is a batched/wide gemm result bitwise equal to the single gemm (and independent of B) for N = 36, 144, 576, 1024?  1 rank, short.
echo "start $(date -Is) node $(hostname) job $SLURM_JOB_ID"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P6=$SCRATCH/large_prompt/team/p6; cd $P6; mkdir -p out
export OMP_NUM_THREADS=3 OMP_TARGET_OFFLOAD=MANDATORY LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576
for pair in "36 144" "256 400" "576 1024"; do set -- $pair
  echo "== N $1 $2 $(date -Is)"; CB_N1=$1 CB_N2=$2 CB_FLOPS=3e11 timeout 600 mpiexec -n 1 $P6/ceil_bench > out/b$1_$2.$SLURM_JOB_ID.out 2>&1; echo rc=$?
done
echo "end $(date -Is)"
