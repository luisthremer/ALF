#!/bin/bash
#SBATCH --job-name=p8b_t4n
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:45:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# G4 (node layout) + G3: 32 ranks, 4/tile, CCS 4, OMP 3, V7, L24, 3-sweep runs (the CPU-time difference of g4n was too noisy: +-10 % drift between runs);
# BLKPACK 0 (old) / 1 / 2 alternating (ABBA) in one job; metric = in-loop TEAM_TIMERS (mean over the 32 ranks, per sweep: sweep, update, stab, cgr, obs) = inloop.py;
# Precision Green of rank 0 (info); old vs new run directories compared bitwise.
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; cd $L0/p8/blk/jobs; source p8lib.sh
export W=$P8/work_t4n; mkdir -p $W
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
runn() { local m=$1 nd=$2 tag=$3
  env $(var V7) $CCS OMPT=3 LSIZE=24 NDELAY=$nd ALFSH_BLKPACK=$m python3 $J/run_ref_sng.py $tag $BIN 32 3 1 > /dev/null
  echo "m=$m nd=$nd"; python3 $PB/inloop.py $W/runs/$tag; pg $tag; }
for nd in 64 32; do for m in 0 1 2 2 1 0; do runn $m $nd n${m}_nd${nd}_$RANDOM; done; done
echo "=== bitwise old vs new (same nd): compare.py tol 0"
for nd in 64 32; do a=$(ls -d $W/runs/n0_nd${nd}_* | head -1); b=$(ls -d $W/runs/n1_nd${nd}_* | head -1); echo "nd=$nd $a $b"; python3 $J/compare.py $a $b --tol 0 | tail -2; done
echo "end $(date -Is)"
