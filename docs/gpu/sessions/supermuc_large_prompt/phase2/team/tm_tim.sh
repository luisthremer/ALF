#!/bin/bash
#SBATCH --job-name=tm_tim
#SBATCH --output=tm_tim.%j.log
#SBATCH --error=tm_tim.%j.err
#SBATCH --time=00:50:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# (1) build aottim (AOT build of the current src, shim with opt-in ALFSH_TIMING), (2) per-routine shim call counts/wall time
# inside ALF, L=24 and L=32 one tile, (3) lapack_bench size scan with queue-wide waits.
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; cd $J
bash build_intel.sh aottim > b_aottim.out 2>&1; ls -l $L0/build_aottim/Prog/ALF.out || exit 1
export W=$L0/work_tim SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 ALFSH_TIMING=1
B=$L0/build_aottim/Prog/ALF.out
for c in "L24 t1 2 6" "L24 t8 2 6" "L32 t8 1 3"; do set -- $c
  for ns in $3 $4; do
    if [ $1 = L32 ]; then export LSIZE=32; else unset LSIZE; fi
    OMPT=${2#t} python3 run_ref_sng.py ${1}_${2}_s$ns $B 1 $ns 1 | tail -1
    echo "--- ${1}_${2}_s$ns"; grep -h ALFSH_TIMING $W/runs/${1}_${2}_s$ns/run.log
  done
done
unset LSIZE
icpx -fsycl -qmkl=parallel -O2 lapack_bench.cpp -o lapack_bench || exit 1
export SYCL_UR_USE_LEVEL_ZERO_V2=1
for t in 1 8; do for n in 512 576 640 768 1024; do OMP_NUM_THREADS=$t MKL_NUM_THREADS=$t ./lapack_bench $n 10 | grep -E "device|geqrf|gemm|getr"; done; done
echo "done $(date -Is)"
