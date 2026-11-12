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
# L=32 single chain (5 variants x 2) and node throughput (32 ranks, 4/tile, CCS 4, OMPT=3; rank-0 in-loop timer AND CPU-Time difference)
source c_lib.sh; export W=$C/work_two
for r in 1 2; do for v in V0 V1 V2 V3 V4; do one L32_${v}_$r $v 2 LSIZE=32; done; done
cpu() { grep -h "CPU Time" $W/runs/$1/info | awk "{print \$NF}"; }
CCS="ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4"
for v in V0 V1 V2 V3 V4; do
  for ns in 1 3; do env $(var $v) $CCS OMPT=3 python3 $J/run_ref_sng.py N_${v}_s$ns $C/build_aottt/Prog/ALF.out 32 $ns 1 > /dev/null; done
  tim N_${v}_s3 | sed "s/RESULT/RESULT(rank0)/"
  python3 -c "x,y=$(cpu N_${v}_s1),$(cpu N_${v}_s3); m=(y-x)/2; print(\"RESULT N_$v CPU-Time-diff s/sweep/chain %.2f sweeps/s/node %.3f per card %.3f\"%(m,32/m,8/m))"
done
echo "done $(date -Is)"
