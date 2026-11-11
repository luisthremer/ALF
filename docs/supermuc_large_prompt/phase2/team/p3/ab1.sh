#!/bin/bash
#SBATCH --job-name=p3_ab1
#SBATCH --output=p3_ab1.%j.log
#SBATCH --error=p3_ab1.%j.err
#SBATCH --time=00:20:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# Gate A microbenchmark: sync vs nowait+depend chain
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P=$SCRATCH/large_prompt/team/p3; cd $P
AOT="-fopenmp-targets=spir64_gen -Xopenmp-target-backend -device -Xopenmp-target-backend pvc"
ifx -O2 -cpp -fiopenmp $AOT -qmkl -I$MKLROOT/include asyncbench.F90 -o ab_sync > b_sync.log 2>&1 || tail b_sync.log
ifx -O2 -cpp -fiopenmp $AOT -qmkl -I$MKLROOT/include -DASYNC asyncbench.F90 -o ab_async > b_async.log 2>&1 || tail b_async.log
ls -l ab_sync ab_async
export OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 OMP_NUM_THREADS=8
RT=$(dirname $(dirname $(which ifx)))/lib
echo "== knobs in the level zero RTL"
for f in $RT/libomptarget.rtl.level_zero*.so $RT/libomptarget*.so; do [ -f $f ] && strings $f | grep -E "^(LIBOMPTARGET_LEVEL_ZERO|LIBOMPTARGET_[A-Z_]*(BATCH|QUEUE|COMMAND|IMMEDIATE|ASYNC|ORDER)[A-Z_]*)" | sort -u; done 2>/dev/null | head -60
echo "== sync"
for r in 1 2 3; do ./ab_sync 300 0 | grep -E "checksum|RESULT"; done
echo "== async (nowait+depend), default"
for r in 1 2 3; do ./ab_async 300 0 | grep -E "checksum|RESULT"; done
echo "done $(date -Is)"
