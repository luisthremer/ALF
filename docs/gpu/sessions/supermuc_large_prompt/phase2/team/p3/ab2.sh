#!/bin/bash
#SBATCH --job-name=p3_ab2
#SBATCH --output=p3_ab2.%j.log
#SBATCH --error=p3_ab2.%j.err
#SBATCH --time=00:25:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# Gate A, mechanisms (c): runtime knobs x {sync, async} x {no host work, ~10 us host work per op}
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P=$SCRATCH/large_prompt/team/p3; cd $P
export OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 OMP_NUM_THREADS=8
run() { # label envs...
  local lab=$1; shift
  for b in sync async; do for hw in 0 7000; do
    r=$(env "$@" timeout 120 ./ab_$b 300 $hw 2>&1 | grep -E "checksum|RESULT|rror|ailed" | tr '\n' ' ' | sed 's/hostwork=\([0-9]*\) checksum= *\([0-9.E+]*\).*RESULT total *\([0-9.]*\) us\/op.*issue-only *\([0-9.]*\).*/hw=\1 chk=\2 total=\3 issue=\4/')
    echo "$lab $b: $r"
  done; done
}
run default A=1
run default2 A=1
for m in sync async async_ordered; do run cmdmode_$m LIBOMPTARGET_LEVEL_ZERO_COMMAND_MODE=$m; done
for bt in copy compute all; do run batch_$bt LIBOMPTARGET_LEVEL_ZERO_COMMAND_BATCH=$bt; done
run batch_all_64 LIBOMPTARGET_LEVEL_ZERO_COMMAND_BATCH=all,64
run nonimm LIBOMPTARGET_LEVEL_ZERO_USE_IMMEDIATE_COMMAND_LIST=0
run nonimm_async LIBOMPTARGET_LEVEL_ZERO_USE_IMMEDIATE_COMMAND_LIST=0 LIBOMPTARGET_LEVEL_ZERO_COMMAND_MODE=async
run nosyclflush LIBOMPTARGET_LEVEL_ZERO_NO_SYCL_FLUSH=1
run helper1 LIBOMP_NUM_HIDDEN_HELPER_THREADS=1
run helper0_threads OMP_NUM_THREADS=1
run multiq LIBOMPTARGET_LEVEL_ZERO_USE_MULTIPLE_COMPUTE_QUEUES=1
run batch_all_async LIBOMPTARGET_LEVEL_ZERO_COMMAND_BATCH=all LIBOMPTARGET_LEVEL_ZERO_COMMAND_MODE=async
echo "done $(date -Is)"
