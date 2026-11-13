#!/bin/bash
#SBATCH --job-name=p8b_g1
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
# usage: sbatch g1.sh PART      PART = v7 | v2
#  v7: G1 bitwise (V7 switches; same binary, same seed, fresh start, ALFSH_BLKPACK unset vs 1, CMP_TOL=0) at L6F (N_delay 8, 4), L6T, L12F, L24F (32, 64);
#      then BLKPACK=2 (NOT bitwise by design): same comparison reported at tol 1e-10, worst relative deviation printed
#  v2: G1 bitwise with V2 (no DEVQR switches), default build vs previous binary (bitwise), ASYNCBLK=1 default vs previous (bitwise), G2 planted defects
PART=${1:-v7}
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; cd $L0/p8/blk/jobs; source p8lib.sh
export W=$P8/work_g1_$PART; mkdir -p $W
echo "host $(hostname) $(date -Is) part $PART"
cmp1() { # tag, V, plant, keys...   (A = BLKPACK unset, B = BLKPACK=$PK via wrap_p.sh)
  local tag=$1 v=$2 pl=$3; shift 3
  env $(var $v) PLANT_N=$pl PACK_N=${PK:-1} W=$W/$v WRAP=$PB/wrap_p.sh CPUBIN=$B_AOT GPUBIN=$B_AOT CMP_TOL=${TOL:-0} SEEDS=$SEEDS_CMP ZE_AFFINITY_MASK=0 python3 $J/run_cmp_gg.py ${v}_$tag "$@" 2>&1 < /dev/null | grep -E "rc=|PASS|FAIL|acceptance|Traceback|rror" | cut -c1-200
  echo "   BLKPACK active in B: $(grep -c 'BLKPACK=1, packed' $W/$v/runs/${v}_$tag/gpu/run.log)  in A: $(grep -c 'BLKPACK=1, packed' $W/$v/runs/${v}_$tag/cpu/run.log)"
}
K6F="L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18"
K6F4="L=6 mz=F scheme=submatrix ndelay=4 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18"
K6T="L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18"
K12F="L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72"
K24F32="L=24 mz=F scheme=submatrix ndelay=32 nsweep=3 nbin=1 U=5.3d0 symm=T beta=24.d0 nphi=288"
K24F64="L=24 mz=F scheme=submatrix ndelay=64 nsweep=3 nbin=1 U=5.3d0 symm=T beta=24.d0 nphi=288"
if [ $PART = v7 ]; then
  echo "######## G1 bitwise V7: $(var V7)"
  for c in "L6F_nd8:$K6F" "L6F_nd4:$K6F4" "L6T_nd8:$K6T" "L12F_nd32:$K12F" "L24F_nd32:$K24F32" "L24F_nd64:$K24F64"; do
    tag=${c%%:*}; echo "## V7 $tag"; cmp1 $tag V7 "" ${c#*:}; done
  echo "######## BLKPACK=2 (not bitwise by design) vs ALFSH_BLKPACK unset, V7, tol 1e-10"
  for c in "L6F_nd8:$K6F" "L12F_nd32:$K12F" "L24F_nd32:$K24F32"; do
    tag=${c%%:*}; echo "## V7 pack2 $tag"; PK=2 TOL=1e-10 cmp1 p2_$tag V7 "" ${c#*:}; done
else
  echo "######## G1 bitwise V2: $(var V2)"
  for c in "L6F_nd8:$K6F" "L6F_nd4:$K6F4" "L6T_nd8:$K6T" "L12F_nd32:$K12F" "L24F_nd32:$K24F32"; do
    tag=${c%%:*}; echo "## V2 $tag"; cmp1 $tag V2 "" ${c#*:}; done
  echo "######## default build vs previous binary (p5 txtnew), no ALFSH switches (bitwise)"
  for a in "g_L6F $K6F" "g_L12F $K12F"; do
    env -u ALFSH_DEVQR -u ALFSH_DEVUNMQR W=$W/gg CPUBIN=$B_PREV GPUBIN=$B_TXT CMP_TOL=0 SEEDS=$SEEDS_CMP ZE_AFFINITY_MASK=0 python3 $J/run_cmp_gg.py $a 2>&1 < /dev/null | grep -E "rc=|PASS|FAIL"
  done
  echo "######## same with ALFSH_ASYNCBLK=1 (V2 switches, BLKPACK unset), bitwise"
  for a in "ga_L6F $K6F" "ga_L12F $K12F"; do
    env $(var V2) W=$W/gga CPUBIN=$B_PREV GPUBIN=$B_TXT CMP_TOL=0 SEEDS=$SEEDS_CMP ZE_AFFINITY_MASK=0 python3 $J/run_cmp_gg.py $a 2>&1 < /dev/null | grep -E "rc=|PASS|FAIL"
  done
  echo "######## G2 planted defects (V7, B = BLKPACK=1 + plant; each must FAIL): 1 wrong M offset (flavour 2), 2 wrong field-value offset, 3 no wait after the gather, 4 M of the wrong flavour, 5 last field value not uploaded"
  for pl in 1 2 3 4 5; do
    for c in "L6F_nd8:$K6F" "L12F_nd32:$K12F"; do
      tag=${c%%:*}; echo "## plant $pl $tag"; cmp1 pl${pl}_$tag V7 $pl ${c#*:} | cut -c1-170
    done
  done
fi
echo "done $(date -Is)"
