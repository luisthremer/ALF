#!/bin/bash
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:59:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; C=$L0/comb; cd $C/jobs
# Phase 5 profile: V4 (staging+async+devqr) one tile, shim per-call timers (ALFSH_TIMING=1) + phase timers, L24 and L32
source c_lib.sh; export W=$C/work_prof
one P_L24_V4 V4 4 ALFSH_TIMING=1
grep -h ALFSH_TIMING $W/runs/P_L24_V4/run.log
one P_L32_V4 V4 2 LSIZE=32 ALFSH_TIMING=1
grep -h ALFSH_TIMING $W/runs/P_L32_V4/run.log
echo "done $(date -Is)"
