#!/bin/bash
#SBATCH --job-name=p3_pb2
#SBATCH --output=p3_pb2.%j.log
#SBATCH --error=p3_pb2.%j.err
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE in ALF itself: L=24 one tile 8 threads 3 sweeps, alternating; confout_0 must be identical
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p3; cd $P
export W=$P/work_pb2 SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0
NEW=$P/build_aotttnew/Prog/ALF.out
one() { # tag staging
  if [ "$2" = 0 ]; then OMPT=8 python3 run_ref_sng.py $1 $NEW 1 3 1 > /dev/null; else LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=$2 OMPT=8 python3 run_ref_sng.py $1 $NEW 1 3 1 > /dev/null; fi
  grep -h "TEAM_TIMERS" $W/runs/$1/run.log | awk -v t=$1 '{n=$(NF-7); printf "RESULT %-10s sweep %6.2f  upd %5.2f  stab %5.2f  cgr %5.2f  obs %5.2f s/sweep\n", t, $(NF-6)/n, ($(NF-5)+$(NF-4))/n, ($(NF-3)+$(NF-2))/n, $(NF-1)/n, $NF/n}'
  grep -h "P3_TIMERS" $W/runs/$1/run.log | awk '{printf "   P3 per sweep: blk_begin %5.2f (dev %5.2f)  blk_end %5.2f  slice %5.2f\n", $(NF-5)/3,$(NF-4)/3,$(NF-3)/3,$(NF-1)/3}'
}
for r in 1 2; do
  one base$r 0; one s1M$r 1048576; one s4M$r 4194304; one s16M$r 16777216
done
for t in s1M1 s4M1 s16M1 base2; do for f in confout_0; do cmp $W/runs/base1/$f $W/runs/$t/$f && echo "IDENTICAL $t $f"; done; done
ls $W/runs/base1 | head -20
echo "done $(date -Is)"
