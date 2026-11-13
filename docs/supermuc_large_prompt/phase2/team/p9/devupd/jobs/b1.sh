#!/bin/bash
#SBATCH --job-name=p9du_b1
#SBATCH --time=00:15:00
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
# p9/devupd feasibility bench (bitwise host vs device of the per-site update arithmetic; one-team block kernel timing)
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P=$SCRATCH/large_prompt/team/p9/devupd; cd $P/bench
FL="-std08 -cpp -O3 -fp-model=fast=2 -no-prec-div -xHost -unroll -finline-functions -fiopenmp -fopenmp-targets=spir64_gen -Xopenmp-target-backend -device -Xopenmp-target-backend pvc"
time ifx $FL -o dubench dubench.F90 2>&1 | tail -20
export OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 OMP_NUM_THREADS=1
for kk in 8 24; do echo "=== K $kk"; DU_K=$kk ./dubench 2>&1 | grep -v "^$"; done
for bb in 32 64; do echo "=== B $bb"; DU_NS=1000 DU_B=$bb ./dubench 2>&1 | grep BLK; done
echo "end $(date -Is)"
