#!/bin/bash
#SBATCH --job-name=p8_build9
#SBATCH --time=00:20:00
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
P=$SCRATCH/large_prompt/team/p8/obs; cd $P/jobs; export SRC=$P/src_head
BLD=$P/build_txt9 bash build_comb.sh txtnew > b_txt9.out 2>&1
ls -l $P/build_txt9/Prog/ALF.out
echo "end $(date -Is)"
