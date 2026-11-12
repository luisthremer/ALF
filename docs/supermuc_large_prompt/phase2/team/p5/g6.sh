#!/bin/bash
#SBATCH --job-name=p5_g6
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
bash $SCRATCH/large_prompt/team/p5/jobs/t_plant2.sh
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p5; cd $P/jobs
source p5lib.sh; export W=$C/work_g6
for pl in nobar onebuf; do echo "[plant $pl] node config"; node NQ_$pl V7 3 ALFSH_DEVQR_LIB=$C/plant/$pl/libalfdevqr.so; pg NQ_$pl; done
echo done
