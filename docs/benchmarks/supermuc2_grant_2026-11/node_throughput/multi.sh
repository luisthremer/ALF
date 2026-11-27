#!/bin/bash
# Several layouts of one side and L in parallel, one full node each (one job slot for several nodes; test allows 5 jobs/user).
# usage: sbatch -N <nodes> -J nt_<side>_L<L>_m -p test|general -t HH:MM:SS multi.sh SIDE L LAYOUT [LAYOUT ...]
#   LAYOUT as in run_layouts.sh (cpu: "T:1", gpu: "W:T"). More layouts than nodes: the extra steps wait for a free node.
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --exclude=i20r02c06s02,i20r02c06s06,i20r03c06s06
SIDE=$1; L=$2; shift 2
J=/hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27/jobs/nt
module load slurm_setup >/dev/null 2>&1
echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NODELIST side $SIDE L $L layouts $* start $(date -Is)"
for lay in "$@"; do
  LAYOUTS_OVERRIDE=$lay NT_LOCAL=1 srun -N1 -n1 -c 224 --cpu-bind=none --mem=0 --exclusive \
    bash $J/run_layouts.sh $SIDE $L > $J/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.${lay/:/_}.log 2>&1 &
  sleep 2
done
wait; echo "end $(date -Is)"
