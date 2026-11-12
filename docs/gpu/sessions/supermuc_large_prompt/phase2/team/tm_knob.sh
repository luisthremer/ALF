#!/bin/bash
#SBATCH --job-name=tm_knob
#SBATCH --output=tm_knob.%j.log
#SBATCH --error=tm_knob.%j.err
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# Runtime-knob scan (no numerics change), L=24 one tile, ref params. Per sweep = (CPU Time(6) - CPU Time(2))/4 from ALF's info.
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; cd $J
export W=$L0/work_knob SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 ALFSH_TIMING=1
B=${B:-$L0/build_aottim/Prog/ALF.out}
cpu() { grep -h "CPU Time" $W/runs/$1/info | awk '{print $NF}'; }
geq() { grep -h "ALFSH_TIMING geqrf " $W/runs/$1/run.log | awk '{print $6}'; }
run() { local tag=$1; shift
  for ns in 2 6; do env "$@" python3 run_ref_sng.py ${tag}_s$ns $B 1 $ns 1 > /dev/null; done
  python3 -c "a,b=$(cpu ${tag}_s2),$(cpu ${tag}_s6); g1,g2=$(geq ${tag}_s2),$(geq ${tag}_s6); print('RESULT %-28s s/sweep %.2f  geqrf s/sweep %.2f  rest %.2f'%('$tag',(b-a)/4,(g2-g1)/4,(b-a)/4-(g2-g1)/4))"
}
for v in ${VARIANTS:-def imm0 immall t4 t2}; do case $v in
  def)    run def OMPT=8 ;;
  imm0)   run imm0 OMPT=8 LIBOMPTARGET_LEVEL_ZERO_USE_IMMEDIATE_COMMAND_LIST=0 ;;
  immall) run immall OMPT=8 LIBOMPTARGET_LEVEL_ZERO_USE_IMMEDIATE_COMMAND_LIST=all ;;
  t4)     run t4 OMPT=4 ;;
  t2)     run t2 OMPT=2 ;;
  t16)    run t16 OMPT=16 ;;
  dbg)    LIBOMPTARGET_DEBUG=1 OMPT=8 python3 run_ref_sng.py dbg_s1 $B 1 1 1 >/dev/null; grep -iE "immediate|command list|plugin|option" $W/runs/dbg_s1/run.log | sort | uniq -c | sort -rn | head -30 ;;
esac; done
echo "done $(date -Is)"
