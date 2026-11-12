#!/bin/bash
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --exclude=i20r02c06s02
#SBATCH --nodes=1
#SBATCH --time=00:40:00
#SBATCH --ntasks-per-node=64
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
source $SCRATCH/large_prompt/team/p7/jobs/p7lib.sh
export W=$P7/work_j2; mkdir -p $W $P7/data
# 32 ranks, 4 per tile, CCS 4, 3 OpenMP threads: unperturbed timing and PLUGIN_PROFILE (1 and 3 sweeps) for V7 and V2
for v in V7 V2; do
  echo "=== $v 32 ranks timing 3 sweeps (ALFSH_TIMING)"; r32 ${v}_n3 $v 3 ALFSH_TIMING=1; tm ${v}_n3; pg ${v}_n3
  echo "=== $v 32 ranks PLUGIN_PROFILE 1 and 3 sweeps"; r32 ${v}_pp1 $v 1 LIBOMPTARGET_PLUGIN_PROFILE=T; r32 ${v}_pp3 $v 3 LIBOMPTARGET_PLUGIN_PROFILE=T; tm ${v}_pp3
done
mkdir -p $P7/data/j2; for d in $W/runs/*; do t=$(basename $d); mkdir -p $P7/data/j2/$t; cp $d/rank*.log $d/info $P7/data/j2/$t/ 2>/dev/null; done
echo "end $(date -Is)"
