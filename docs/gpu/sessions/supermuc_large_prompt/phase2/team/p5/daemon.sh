#!/bin/bash
#SBATCH --job-name=p5_daemon
#SBATCH --output=daemon.%j.log
#SBATCH --error=daemon.%j.err
#SBATCH --time=00:58:00
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# Dev-loop helper: runs queued scripts team/p5/queue/NAME.sh one after the other (output queue/NAME.out, marker NAME.done),
# so that edit-compile-test cycles do not pay a Slurm queue wait each. Leaves after 12 min without work or at 56 min.
# submit: sbatch -p preempt --qos=preempt daemon.sh
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P=$SCRATCH/large_prompt/team/p5; Qd=$P/queue; mkdir -p $Qd; cd $P
echo "daemon start $(hostname) $(date -Is)"; t0=$(date +%s); idle=$t0
while :; do
  now=$(date +%s)
  [ $((now - t0)) -gt 3360 ] && break
  [ $((now - idle)) -gt 720 ] && break
  f=$(ls $Qd/*.sh 2>/dev/null | head -1)
  if [ -n "$f" ]; then
    b=${f%.sh}; mv $f $b.run
    echo "run $(basename $b) $(date -Is)"
    ( cd $P; bash $b.run > $b.out 2>&1 ); mv $b.run $b.done
    idle=$(date +%s)
  else sleep 3; fi
done
echo "daemon end $(date -Is)"
