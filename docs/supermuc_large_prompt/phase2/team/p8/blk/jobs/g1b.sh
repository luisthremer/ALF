#!/bin/bash
#SBATCH --job-name=p8b_g1b
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
# G1 (bitwise, same binary, same seed, fresh start, ALFSH_BLKPACK unset vs =1, CMP_TOL=0) and G2 (planted defects must FAIL)
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; cd $L0/p8/blk/jobs; source p8lib.sh
export W=$P8/work_g1b; mkdir -p $W
echo "host $(hostname) $(date -Is)"
cmp1() { # tag, V, extra env words, keys...   (A = BLKPACK unset, B = BLKPACK=1 via wrap_p.sh)
  local tag=$1 v=$2 pl=$3; shift 3
  env $(var $v) PLANT_N=$pl W=$W/$v WRAP=$PB/wrap_p.sh CPUBIN=$B_AOT GPUBIN=$B_AOT CMP_TOL=0 SEEDS=$SEEDS_CMP ZE_AFFINITY_MASK=0 python3 $J/run_cmp_gg.py ${v}_$tag "$@" 2>&1 < /dev/null | grep -E "rc=|PASS|FAIL|acceptance|Traceback|rror"
  echo "   BLKPACK active in B: $(grep -c 'BLKPACK=1, packed' $W/$v/runs/${v}_$tag/gpu/run.log)  in A: $(grep -c 'BLKPACK=1, packed' $W/$v/runs/${v}_$tag/cpu/run.log)"
}
CASES_S="L6F_nd8:L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18
L6F_nd4:L=6 mz=F scheme=submatrix ndelay=4 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18
L6T_nd8:L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18
L12F_nd32:L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72
L24F_nd32:L=24 mz=F scheme=submatrix ndelay=32 nsweep=3 nbin=1 U=5.3d0 symm=T beta=24.d0 nphi=288
L24F_nd64:L=24 mz=F scheme=submatrix ndelay=64 nsweep=3 nbin=1 U=5.3d0 symm=T beta=24.d0 nphi=288"
echo "######## G1 bitwise: ALFSH_BLKPACK unset (A) vs 1 (B), same binary $B_AOT"
for v in V7 V2; do  # (job g12 already did V7/V2 L6F_nd8: the here-string loop was cut by mpiexec reading stdin)
  echo "#### $v: $(var $v)"
  echo "$CASES_S" | while IFS=: read tag keys; do echo "## $v $tag"; cmp1 $tag $v "" $keys; done
done
echo "done $(date -Is)"
