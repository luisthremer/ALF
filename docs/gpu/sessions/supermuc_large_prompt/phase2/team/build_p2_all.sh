#!/bin/bash
#SBATCH --job-name=p2_build
#SBATCH --output=p2_build.%j.log
#SBATCH --error=p2_build.%j.err
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
cd $SCRATCH/large_prompt/team/jobs
for v in gpu-txt aottim gpu-txt_def; do bash build_p2.sh $v > p2b_$v.log 2>&1 & done
wait
ls -l $SCRATCH/large_prompt/team/build_p2_*/Prog/ALF.out
