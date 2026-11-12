#!/bin/bash
#SBATCH --job-name=p5_g3b
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:59:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p5; cd $P/jobs
source p5lib.sh; export W=$C/work_g3b
echo "host $(hostname) $(date -Is)"
node8() { local t=$1 v=$2 ns=$3; shift 3; env $(var $v) "$@" OMPT=14 python3 $J/run_ref_sng.py $t $C/build_aottt/Prog/ALF.out 8 $ns 1 > /dev/null; }
echo "######## node throughput by CPU-Time difference (1 vs 3 sweeps), alternating V4 / V7, 2 pairs; layout A: 32 ranks, 4/tile, CCS 4, OMP 3 (8 chains/card); layout B: 8 ranks, 1/tile, OMP 14, no CCS (2 chains/card)"
for r in 1 2; do for v in V4 V7; do
  node TA_${v}_${r}_s1 $v 1; node TA_${v}_${r}_s3 $v 3; pg TA_${v}_${r}_s3
  python3 -c "x,y=$(cpu TA_${v}_${r}_s1),$(cpu TA_${v}_${r}_s3); m=(y-x)/2; print(\"RESULT layoutA_32x4 $v pair $r: CPU-Time-diff s/sweep/chain %.2f  sweeps/s per card %.3f\"%(m,8/m))"
  node8 TB_${v}_${r}_s1 $v 1; node8 TB_${v}_${r}_s3 $v 3; pg TB_${v}_${r}_s3
  python3 -c "x,y=$(cpu TB_${v}_${r}_s1),$(cpu TB_${v}_${r}_s3); m=(y-x)/2; print(\"RESULT layoutB_8x1 $v pair $r: CPU-Time-diff s/sweep/chain %.2f  sweeps/s per card %.3f\"%(m,2/m))"
done; done
echo "done $(date -Is)"
