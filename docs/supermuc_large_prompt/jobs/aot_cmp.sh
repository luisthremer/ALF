#!/bin/bash
#SBATCH --job-name=lp_aotc
#SBATCH --output=lp_aotc.%j.log
#SBATCH --error=lp_aotc.%j.err
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# AOT (spir64_gen -device pvc) vs JIT (spir64): same source, L=24, 32 ranks (4 per tile, CCS), 3 sweeps in 1 bin; plus a same-seed 4-rank check against JIT
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
L=$SCRATCH/large_prompt; J=$L/jobs; cd $J
export W=$L/work_aotc SEEDS=$L/ref_L24/seeds PARAMS=$L/ref_L24/parameters OMP_TARGET_OFFLOAD=MANDATORY ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4; mkdir -p $W
echo "host $(hostname) $(date -Is)"
for v in jit aot jit aot; do
  B=$L/build_mpi/Prog/ALF.out; [ $v = aot ] && B=$L/build_aot/Prog/ALF.out
  python3 run_ref_sng.py aotc_${v}_$RANDOM $B 32 3 1 | tail -1 | sed "s/^/[$v] /"
done
echo "=== same seeds, 4 ranks, 2 sweeps: JIT vs AOT outputs"
python3 run_ref_sng.py aotc_j4 $L/build_mpi/Prog/ALF.out 4 2 1 | tail -1; python3 run_ref_sng.py aotc_a4 $L/build_aot/Prog/ALF.out 4 2 1 | tail -1
for f in Ener_scal; do for r in aotc_j4 aotc_a4; do echo "$r $(grep -h 'Acceptance ' $W/runs/$r/info | head -1)"; done; done
echo "done $(date -Is)"
