# Node configuration (32 ranks, 4 per tile, ZEX_NUMBER_OF_CCS 4, OMPT=3, staging 1 MB): Precision Green (rank 0) and per-card sweeps/s by CPU-Time difference,
# devqr off/on alternating, 2 pairs (as team/comb/jobs/c_two.sh). usage (queue script): bash node.sh
source $SCRATCH/large_prompt/team/p4/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p4; cd $P/jobs
export W=$P/work_node SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters OMP_TARGET_OFFLOAD=MANDATORY ALFSH_QRNB=16 LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576
B=$P/build_aotttnew/Prog/ALF.out
CCS="ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4"
cpu() { grep -h "CPU Time" $W/runs/$1/info | awk '{print $NF}'; }
for r in 1 2; do for v in off on; do
  if [ $v = on ]; then D=1; else D=0; fi
  for ns in 1 3; do env ALFSH_DEVQR=$D $CCS OMPT=3 python3 run_ref_sng.py N_${v}${r}_s$ns $B 32 $ns 1 | tail -1; done
  python3 -c "x,y=$(cpu N_${v}${r}_s1),$(cpu N_${v}${r}_s3); m=(y-x)/2; print('RESULT devqr=$v pair=$r CPU-Time-diff s/sweep/chain %.2f  sweeps/s/node %.3f  per card %.3f'%(m,32/m,8/m))"
  echo "   $(grep -h 'Precision Green' $W/runs/N_${v}${r}_s3/info | head -1)"
done; done
echo "done $(date -Is)"
