#!/bin/bash -l
#SBATCH -J ib_s2
#SBATCH -A pn93la
#SBATCH -p preempt
#SBATCH --qos=preempt
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH -t 00:20:00
#SBATCH -x i20r02c06s02
#SBATCH -o %x.%j.log
#SBATCH -e %x.%j.err
# Update-block emulation with the real flush structure (IB_BLK_REAL: per flavour gemm1 N x k x k + gemm2 N x N x k, one G pass),
# k = 16 and 32; 1 process alone (W = 1..16 batched) and 4 processes per tile with CCS 4 (production layout).
D=$SCRATCH/large_prompt/batched
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
B=$D/bin/ib_bench; O=$D/out/$SLURM_JOB_ID; mkdir -p $O
export ZE_FLAT_DEVICE_HIERARCHY=FLAT IB_BLK_REAL=1
echo "start $(date -Is) $(hostname)"
for k in 16 32; do
  ZE_AFFINITY_MASK=0 IB_BLK_K=$k $B block 0.4 > $O/r1_k$k.txt 2>&1
  pids=""
  for p in 1 2 3 4; do
    ZE_AFFINITY_MASK=0 ZEX_NUMBER_OF_CCS=0:4 IB_BLK_K=$k taskset -c $((p*8))-$((p*8+7)) $B block 0.4 > $O/blk4_ccs4_k${k}_$p.txt 2>&1 &
    pids="$pids $!"
  done
  wait $pids
done
echo "end $(date -Is)"
