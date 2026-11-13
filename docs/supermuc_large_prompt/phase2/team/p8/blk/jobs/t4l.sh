#!/bin/bash
#SBATCH --job-name=p8b_t4l
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
# G4 one chain on tile 0 (in-loop phase + block timers, 3 sweeps), V7, L24, nd 32: BLKPACK 0/1/2 alternating (ABBA), then ablations of BLKPACK=1:
# ALFSB_GDEV=1 (gather via a device buffer + copy), ALFSB_NOENQ=1 (queue::parallel_for/memcpy instead of the enqueue functions)
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; cd $L0/p8/blk/jobs; source p8lib.sh
export W=$P8/work_t4l; mkdir -p $W
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
r1() { local t=$1 nd=$2; shift 2; one $t V7 3 LSIZE=24 NDELAY=$nd "$@" > $W/$t.out 2>&1; echo "== $t nd=$nd $*"; grep -h -E "RESULT.*upd|Precision" $W/$t.out; p3 $t 3; }
for m in 0 1 2 2 1 0; do r1 c${m}_nd32_$RANDOM 32 ALFSH_BLKPACK=$m; done
r1 abl_gdev_$RANDOM 32 ALFSH_BLKPACK=1 ALFSB_GDEV=1
r1 abl_noenq_$RANDOM 32 ALFSH_BLKPACK=1 ALFSB_NOENQ=1
r1 abl_both_$RANDOM 32 ALFSH_BLKPACK=1 ALFSB_NOENQ=1 ALFSB_GDEV=1
for m in 0 1 0 1; do r1 e${m}_nd32_$RANDOM 32 ALFSH_BLKPACK=$m; done
echo "end $(date -Is)"
