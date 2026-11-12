#!/bin/bash
#SBATCH --job-name=p5_g3a
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
source p5lib.sh; export W=$C/work_g3a
echo "host $(hostname) $(date -Is)"
echo "######## alternating V4 (async + devqr 1) vs V7 (async + devqr 2 + devunmqr), single tile (one chain)"
for r in 1 2 3; do for v in V4 V7; do one T_L24_${v}_$r $v 4; done; done
for r in 1 2; do for v in V4 V7; do one T_L32_${v}_$r $v 2 LSIZE=32; done; done
echo "done $(date -Is)"
