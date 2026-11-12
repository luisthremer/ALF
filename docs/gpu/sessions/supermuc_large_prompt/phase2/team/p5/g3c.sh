#!/bin/bash
#SBATCH --job-name=p5_g3c
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:59:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p5; cd $P/jobs
source p5lib.sh; export W=$C/work_g3c
for r in 2; do for v in V4 V7; do one T_L32_${v}_$r $v 2 LSIZE=32; done; done
echo done
