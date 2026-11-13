#!/bin/bash
#SBATCH --job-name=p8_b7
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
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P8=$SCRATCH/large_prompt/team/p8/obs; cd $P8/bench
FL="-std08 -cpp -O3 -fp-model=fast=2 -no-prec-div -xHost -unroll -finline-functions -fiopenmp -fopenmp-targets=spir64_gen -Xopenmp-target-backend -device -Xopenmp-target-backend pvc"
for q in 5 6; do ifx $FL -DQBV=$q -o obsbench_q$q obsbench.F90 2>&1 | grep -E "error" ; done
export OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 OMP_NUM_THREADS=1
for q in 5 6; do echo "################ QB=$q"; for cfg in "24 1" "24 2" "32 1" "6 1" "12 1"; do set -- $cfg; echo "=== L $1 Norb $2"; OB_L=$1 OB_NORB=$2 OB_NREP=10 ./obsbench_q$q 2>&1 | grep -E "^OB (N=|TIME|chunk)|rror"; done; done
echo "end $(date -Is)"
