#!/bin/bash
#SBATCH --job-name=tm_pti
#SBATCH --output=tm_pti.%j.log
#SBATCH --error=tm_pti.%j.err
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# PTI kernel trace of one L=24 chain (one tile, 8 host threads), 2 and 6 sweeps: per-sweep kernel count and device time by difference.
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; cd $J
PTI=/lrz/sys/intel/oneapi_2025.2.1/pti/latest
icpx -O2 -shared -fPIC ptitrace.cpp -I$PTI/include -L$PTI/lib -lpti_view -Wl,-rpath,$PTI/lib -o libptitrace.so || exit 1
export W=$L0/work_pti SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 ALFSH_TIMING=1
B=$L0/build_aottim/Prog/ALF.out
for c in "1 1" "0 1" "${TL:-1} 2" "${TL:-1} 6"; do set -- $c; ns=$2
  ZE_ENABLE_TRACING_LAYER=$1 OMPT=8 python3 run_ref_sng.py pti_s$ns $B 1 $ns 1 $J/pti_wrap.sh | tail -1
  echo "=== s$ns CPU Time $(grep -h 'CPU Time' $W/runs/pti_s$ns/info)"; grep -hE "^PTI|ALFSH_TIMING (geqrf|unmqr|getrf) " $W/runs/pti_s$ns/run.log
done
echo "done $(date -Is)"
