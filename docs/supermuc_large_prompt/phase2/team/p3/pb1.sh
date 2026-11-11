#!/bin/bash
#SBATCH --job-name=p3_pb1
#SBATCH --output=p3_pb1.%j.log
#SBATCH --error=p3_pb1.%j.err
#SBATCH --time=00:50:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# Phase 3b step 2: per-block time breakdown of the submatrix update phase (P3_TIMERS), L=24 one tile, 8 host threads
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p3; J=$L0/jobs; cd $P
bash build_intel_p3.sh aotttnew > b_aotttnew.out 2>&1; ls -l $P/build_aotttnew/Prog/ALF.out || { tail -30 b_aotttnew.out; exit 1; }
export W=$P/work_pb1 SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0
NEW=$P/build_aotttnew/Prog/ALF.out
cp $J/run_ref_sng.py .
for r in 1 2; do
  OMPT=8 python3 run_ref_sng.py pb$r $NEW 1 3 1
  grep -h "TEAM_TIMERS\|P3_TIMERS" $W/runs/pb$r/run.log
done
echo "done $(date -Is)"
