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
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-node=64
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
source $SCRATCH/large_prompt/team/p7/jobs/p7lib.sh
export W=$P7/work_j7; mkdir -p $W $P7/data
# N_delay scan (submatrix, V7): one chain on tile 0 (3 sweeps) and node layout 32 ranks (3 sweeps); order alternated, same seeds
echo "=== single tile"
for nd in 32 64 128 32 64 128; do NDELAY=$nd r1 s_nd${nd}_$RANDOM V7 3 >/dev/null; done
for d in $W/runs/s_nd*; do echo "$(basename $d) $(tm $(basename $d))"; done
echo "=== node 32 ranks"
for nd in 32 64 128 32 64; do NDELAY=$nd r32 n_nd${nd}_$RANDOM V7 3 >/dev/null; done
for d in $W/runs/n_nd*; do t=$(basename $d); echo "$t $(tm $t) $(pg $t)"; done
mkdir -p $P7/data/j7; for d in $W/runs/*; do t=$(basename $d); mkdir -p $P7/data/j7/$t; cp $d/rank0.log $d/info $P7/data/j7/$t/ 2>/dev/null; done
echo "end $(date -Is)"
