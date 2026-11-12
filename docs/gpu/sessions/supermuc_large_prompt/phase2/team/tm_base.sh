#!/bin/bash
#SBATCH --job-name=tm_base
#SBATCH --output=tm_base.%j.log
#SBATCH --error=tm_base.%j.err
#SBATCH --time=00:50:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# Baseline timing, reference parameters (ref_L24, N_delay=32, Nwrap=10). Marginal s/sweep = (wall(4 sweeps) - wall(1 sweep))/3.
# usage: sbatch tm_base.sh [BUILD (default aot)] [sections, default "tile card L32 node8 node32"]
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
B=${1:-aot}; SECS=${2:-tile card L32 node8 node32}
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; cd $J; BIN=$L0/build_$B/Prog/ALF.out
export W=$L0/work_base_$B SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters OMP_TARGET_OFFLOAD=MANDATORY; mkdir -p $W
echo "host $(hostname) $(date -Is) bin $BIN"; xpu-smi discovery 2>&1 | grep -E "Device Name|Device ID" | head -4
marg() { # tag nranks n1 n2 [env...]
  local t=$1 nr=$2 a=$3 b=$4; shift 4
  local w1=$(env "$@" python3 run_ref_sng.py ${t}_s$a $BIN $nr $a 1 | sed -n 's/.*wall=\([0-9.]*\)s.*/\1/p')
  local w2=$(env "$@" python3 run_ref_sng.py ${t}_s$b $BIN $nr $b 1 | sed -n 's/.*wall=\([0-9.]*\)s.*/\1/p')
  python3 -c "w1,w2=$w1,$w2; m=(w2-w1)/($b-$a); print('RESULT %-14s ranks=%-3d wall(%d)=%.1f wall(%d)=%.1f marginal s/sweep/chain=%.2f  sweeps/s total=%.3f'%('$t',$nr,$a,w1,$b,w2,m,$nr/m))"
}
for s in $SECS; do case $s in
  tile)   marg L24_tile_t1 1 1 4 OMPT=1 ZE_AFFINITY_MASK=0
          marg L24_tile_t8 1 1 4 OMPT=8 ZE_AFFINITY_MASK=0 ;;
  card)   marg L24_card_t8 1 1 4 OMPT=8 ZE_FLAT_DEVICE_HIERARCHY=COMPOSITE ZE_AFFINITY_MASK=0 ;;
  L32)    marg L32_tile_t8 1 1 3 OMPT=8 ZE_AFFINITY_MASK=0 LSIZE=32 NDELAY=32 ;;
  node8)  marg L24_node8_t8 8 1 4 OMPT=8 ;;
  warm)   OMPT=1 ZE_AFFINITY_MASK=0 python3 run_ref_sng.py warm $BIN 1 1 1 | tail -1 ;;
  t1)     marg L24_tile_t1 1 2 6 OMPT=1 ZE_AFFINITY_MASK=0 ;;
  t2)     marg L24_tile_t2 1 2 6 OMPT=2 ZE_AFFINITY_MASK=0 ;;
  t8)     marg L24_tile_t8 1 2 6 OMPT=8 ZE_AFFINITY_MASK=0 ;;
  c1)     marg L24_card_t1 1 2 6 OMPT=1 ZE_FLAT_DEVICE_HIERARCHY=COMPOSITE ZE_AFFINITY_MASK=0 ;;
  L32t1)  marg L32_tile_t1 1 1 4 OMPT=1 ZE_AFFINITY_MASK=0 LSIZE=32 NDELAY=32 ;;
  n8t1)   marg L24_node8_t1 8 2 6 OMPT=1 ;;
  n16t1)  marg L24_node16_t1 16 2 5 OMPT=1 ;;
  node32) marg L24_node32_t3 32 1 3 OMPT=3 ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4 ;;
esac; done
for r in $W/runs/*/info; do echo "$r: $(grep -h "Precision Green" $r)"; done 2>/dev/null
echo "done $(date -Is)"
