#!/bin/bash
#SBATCH --job-name=tm_phase
#SBATCH --output=tm_phase.%j.log
#SBATCH --error=tm_phase.%j.err
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# Wall time per sweep phase (ALF_TEAM_TIMERS build aottt) + shim LAPACK timers, one tile, 8 host threads. Per sweep by difference.
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; cd $J
bash build_intel.sh aottt > b_aottt.out 2>&1; ls -l $L0/build_aottt/Prog/ALF.out || { tail -30 b_aottt.out; exit 1; }
export W=$L0/work_phase SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 ALFSH_TIMING=1
B=$L0/build_aottt/Prog/ALF.out
for c in "L24 2" "L24 6" "L32 1" "L32 3"; do set -- $c
  if [ $1 = L32 ]; then export LSIZE=32; else unset LSIZE; fi
  OMPT=8 python3 run_ref_sng.py ${1}_s$2 $B 1 $2 1 | tail -1
  echo "--- ${1}_s$2 $(grep -h 'CPU Time' $W/runs/${1}_s$2/info)"; grep -hE "TEAM_TIMERS|ALFSH_TIMING (geqrf|unmqr|getrf) " $W/runs/${1}_s$2/run.log
done
echo "done $(date -Is)"
