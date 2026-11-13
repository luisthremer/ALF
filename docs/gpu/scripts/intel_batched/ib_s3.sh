#!/bin/bash -l
#SBATCH -J ib_s3
#SBATCH -A pn93la
#SBATCH -p preempt
#SBATCH --qos=preempt
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH -t 00:15:00
#SBATCH -x i20r02c06s02
#SBATCH -o %x.%j.log
#SBATCH -e %x.%j.err
# Observables corr kernel (SpinMz form): today's mapping (A) vs column-coalesced two-pass (B), L = 24 and 32;
# 1 process alone and 4 processes per tile with CCS 4.
D=$SCRATCH/large_prompt/batched
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
B=$D/bin/ib_bench; O=$D/out/$SLURM_JOB_ID; mkdir -p $O
export ZE_FLAT_DEVICE_HIERARCHY=FLAT
echo "start $(date -Is) $(hostname)"
ZE_AFFINITY_MASK=0 $B corr 0.5 > $O/corr_r1.txt 2>&1
pids=""
for p in 1 2 3 4; do
  ZE_AFFINITY_MASK=0 ZEX_NUMBER_OF_CCS=0:4 taskset -c $((p*8))-$((p*8+7)) $B corr 0.5 > $O/corr_ccs4_$p.txt 2>&1 &
  pids="$pids $!"
done
wait $pids
cat $O/corr_r1.txt $O/corr_ccs4_*.txt
echo "end $(date -Is)"
