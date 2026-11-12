#!/bin/bash
#SBATCH --job-name=tm_p1
#SBATCH --output=tm_p1.%j.log
#SBATCH --error=tm_p1.%j.err
#SBATCH --time=00:50:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=224
# Phase 1: oneMKL batched geqrf/getrf (ALFSH_BATCH=1) vs plain. Rebuild aottim + gpu-txt (shim changed); timing L24/L32 one tile
# (ALF CPU Time, 2/6 and 1/3 sweeps), shim per-call times, Precision Green; same-seed comparison vs cpu-txt with ALFSH_BATCH=1.
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; cd $J
bash build_intel.sh aottim > b_aottim.out 2>&1 & bash build_intel.sh gpu-txt > b_gpu-txt.out 2>&1 & wait
ls -l $L0/build_aottim/Prog/ALF.out $L0/build_gpu-txt/Prog/ALF.out || exit 1
export W=$L0/work_p1 SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 ALFSH_TIMING=1
B=$L0/build_aottim/Prog/ALF.out
cpu() { grep -h "CPU Time" $W/runs/$1/info | awk '{print $NF}'; }
for c in "L24 2 6" "L32 1 3"; do set -- $c
  if [ $1 = L32 ]; then export LSIZE=32; else unset LSIZE; fi
  for bt in 0 1; do
    for ns in $2 $3; do ALFSH_BATCH=$bt OMPT=8 python3 run_ref_sng.py ${1}_b${bt}_s$ns $B 1 $ns 1 > /dev/null; done
    a=$(cpu ${1}_b${bt}_s$2); b=$(cpu ${1}_b${bt}_s$3)
    echo "RESULT $1 batch=$bt s/sweep $(python3 -c "print('%.2f'%(($b-$a)/($3-$2)))")  $(grep -h 'Precision Green' $W/runs/${1}_b${bt}_s$3/info)"
    grep -hE "ALFSH_TIMING (geqrf|getrf)" $W/runs/${1}_b${bt}_s$3/run.log; grep -h "alfsh:" $W/runs/${1}_b${bt}_s$3/run.log | head -3
  done
done
unset LSIZE
echo "=== same-seed comparison, ALFSH_BATCH=1"
ALFSH_BATCH=1 bash tm_cmp.sh gpu-txt p1batch nodefect 2>&1 | grep -vE "^host"
echo "done $(date -Is)"
