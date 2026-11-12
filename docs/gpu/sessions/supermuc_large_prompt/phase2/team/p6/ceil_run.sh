#!/bin/bash
#SBATCH --job-name=p6_ceil
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
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:40:00
# Gate 1: zgemm ceiling on one tile. 1 rank (tile 0), 8 ranks (1/tile), 32 ranks (4/tile, CCS 4, as production).
echo "start $(date -Is) node $(hostname) job $SLURM_JOB_ID"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P6=$SCRATCH/large_prompt/team/p6; cd $P6; mkdir -p out
export OMP_NUM_THREADS=3 OMP_TARGET_OFFLOAD=MANDATORY
export LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576
CCS="ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4"
cr() { local tag=$1 n=$2; shift 2
  echo "== CASE $tag ranks $n [$*] $(date -Is)"
  ( for kv in "$@"; do export "$kv"; done; timeout 900 mpiexec -n $n $P6/ceil_bench > out/$tag.$SLURM_JOB_ID.out 2>&1; echo "rc=$?" )
  grep -c "^CEIL rank" out/$tag.$SLURM_JOB_ID.out
}
cr r1 1
cr r8 8
cr r8ccs 8 $CCS
cr r32 32 $CCS CB_FLOPS=1.5e12
echo "end $(date -Is)"
