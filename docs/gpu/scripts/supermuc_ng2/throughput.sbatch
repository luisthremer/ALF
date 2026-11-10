#!/bin/bash
#SBATCH -J alf_perf7
#SBATCH --account=pn93la
#SBATCH --time=00:50:00
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=14
#SBATCH -o /hppfs/scratch/0C/di38zuw/alf/perf7_%j.out
module load slurm_setup
source ~/alf-gpu/env251.sh
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=14
export RUNROOT=/hppfs/scratch/0C/di38zuw/alf/perf7; mkdir -p $RUNROOT
BIN=$HOME/alf-gpu/bin/ALF-prodaot.out
RC=$HOME/alf-gpu/run_case.sh
export I_MPI_PIN_DOMAIN=auto I_MPI_PIN_ORDER=compact KMP_AFFINITY=granularity=fine,compact,1,0
for LL in 24 28 32; do echo "#### L=$LL 8 ranks x14 thr ND=128"; LL=$LL ND=128 bash $RC $BIN 8 1 2 L${LL} 2>&1 | grep RESULT; grep -m3 -E "GPU mirror total|Precision Green|Acceptance  " $RUNROOT/L${LL}/info | cut -c1-120; done
echo "#### done"
