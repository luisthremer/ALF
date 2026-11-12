#!/bin/bash
#SBATCH --job-name=p2_tim
#SBATCH --output=p2_tim.%j.log
#SBATCH --error=p2_tim.%j.err
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# prototyper p2 timing: old (team/build_aottim) vs new (build_p2_aottim) alternately, 3 reps, L=24 one tile, 8 host threads
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; cd $J
export W=$L0/work_p2tim SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 ALFSH_TIMING=1
echo "host $(hostname) $(date -Is)"
cpu() { grep -h "CPU Time" $W/runs/$1/info | awk '{print $NF}'; }
run() { local tag=$1 B=$2
  for ns in 2 6; do OMPT=8 python3 run_ref_sng.py ${tag}_s$ns $B 1 $ns 1 > /dev/null; done
  python3 -c "a,b=$(cpu ${tag}_s2),$(cpu ${tag}_s6); print('RESULT %-10s s/sweep %.2f'%('$tag',(b-a)/4))"
}
for r in 1 2 3; do
  run old_r$r $L0/build_aottim/Prog/ALF.out
  run new_r$r $L0/build_p2_aottim/Prog/ALF.out
done
for r in $W/runs/*_s6/info; do echo "$r: $(grep -h 'Precision Green' $r) $(grep -h 'Acceptance  ' $r)"; done
echo "done $(date -Is)"
