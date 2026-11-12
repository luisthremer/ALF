#!/bin/bash
#SBATCH --job-name=p6_ceil3
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
# Gate 1c: rank-k (submatrix flush) gemms, single vs strided batch; 1 rank, 8 ranks (1/tile), 32 ranks (4/tile, CCS 4). N = 576 only.
echo "start $(date -Is) node $(hostname) job $SLURM_JOB_ID"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P6=$SCRATCH/large_prompt/team/p6; cd $P6; mkdir -p out
export OMP_NUM_THREADS=3 OMP_TARGET_OFFLOAD=MANDATORY LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 CB_N1=576 CB_N2=0 CB_FLOPS=3e11
CCS="ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4"
echo "== r1 $(date -Is)";  timeout 600 mpiexec -n 1 $P6/ceil_bench > out/k1.$SLURM_JOB_ID.out 2>&1; echo rc=$?
echo "== r8 $(date -Is)";  timeout 600 mpiexec -n 8 $P6/ceil_bench > out/k8.$SLURM_JOB_ID.out 2>&1; echo rc=$?
echo "== r32 $(date -Is)"; env $CCS timeout 900 mpiexec -n 32 $P6/ceil_bench > out/k32.$SLURM_JOB_ID.out 2>&1; echo rc=$?
echo "end $(date -Is)"
