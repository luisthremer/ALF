#!/bin/bash
#SBATCH --job-name=p3_ab3
#SBATCH --output=p3_ab3.%j.log
#SBATCH --error=p3_ab3.%j.err
#SBATCH --time=00:15:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# Gate A follow-up: real host work (sin loop) between launches; variant without depend clauses (order NOT guaranteed by OpenMP)
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P=$SCRATCH/large_prompt/team/p3; cd $P
AOT="-fopenmp-targets=spir64_gen -Xopenmp-target-backend -device -Xopenmp-target-backend pvc"
F="-O2 -cpp -fiopenmp $AOT -qmkl -I$MKLROOT/include"
ifx $F asyncbench.F90 -o ab_sync > b_sync.log 2>&1
ifx $F -DASYNC asyncbench.F90 -o ab_async > b_async.log 2>&1
ifx $F -DASYNC -DNODEP asyncbench.F90 -o ab_nodep > b_nodep.log 2>&1
ls -l ab_sync ab_async ab_nodep
export OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 OMP_NUM_THREADS=8
for hw in 0 100 400 1000; do for b in sync async nodep; do for r in 1 2; do
  echo "$b hw=$hw: $(timeout 120 ./ab_$b 300 $hw 2>&1 | grep -E "checksum|RESULT|rror" | tr '\n' ' ' | cut -c1-300)"
done; done; done
echo "== host-only cost of the sin loop (sync binary, no device work is not separable; use the hw=0 vs hw=N difference)"
echo "done $(date -Is)"
