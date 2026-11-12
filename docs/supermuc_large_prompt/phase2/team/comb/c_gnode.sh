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
# Gate after the devqr race fix (29005a0f): node configuration (32 ranks, 4/tile, CCS 4, OMPT=3), Precision Green of all ranks,
# per-card sweeps/s by CPU-Time difference, V2 (async) vs V4 (async+devqr) alternating, 2 pairs.
# Planted defect: the pre-fix libalfdevqr.so (race) must FAIL the Precision Green check (expected max >> 1e-8).
source c_lib.sh; export W=$C/work_gnode
cpu() { grep -h "CPU Time" $W/runs/$1/info | awk "{print \$NF}"; }
pg() { grep -h "Precision Green" $W/runs/$1/info | awk -v t=$1 '{m=$NF; printf "PG %-12s mean %s max %s -> %s\n", t, $(NF-1), m, (m+0<1e-8)?"PASS":"FAIL"}'; }
CCS="ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4"
echo "######## PLANTED: pre-fix devqr lib, node config, 3 sweeps"
env $(var V4) ALFSH_DEVQR_LIB=$C/plant/libalfdevqr_prefix.so $CCS OMPT=3 python3 $J/run_ref_sng.py P_pre_s3 $C/build_aottt/Prog/ALF.out 32 3 1 > /dev/null; pg P_pre_s3
for r in 1 2; do for v in V2 V4; do
  for ns in 1 3; do env $(var $v) $CCS OMPT=3 python3 $J/run_ref_sng.py N_${v}_${r}_s$ns $C/build_aottt/Prog/ALF.out 32 $ns 1 > /dev/null; done
  pg N_${v}_${r}_s3
  python3 -c "x,y=$(cpu N_${v}_${r}_s1),$(cpu N_${v}_${r}_s3); m=(y-x)/2; print(\"RESULT N_${v}_$r CPU-Time-diff s/sweep/chain %.2f sweeps/s/node %.3f per card %.3f\"%(m,32/m,8/m))"
done; done
echo "done $(date -Is)"
