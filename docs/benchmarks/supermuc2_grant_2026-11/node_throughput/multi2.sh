#!/bin/bash
# Mixed runs in one job, one full node per run (srun steps; a step waits when all nodes are busy, so list long runs first).
# usage: [BETA=L RUNDIR=runs_betaL] sbatch -N <nodes> -J <name> -p test -t HH:MM:SS multi2.sh ITEM [ITEM ...]
#   ITEM = SIDE:L:LAYOUT:NBIN, LAYOUT as in run_layouts.sh with "_" instead of ":", e.g. cpu:36:2_1:3 gpub:24:4_1:6
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --exclude=i20r02c06s02,i20r02c06s06,i20r03c06s06
J=/hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27/jobs/nt
module load slurm_setup >/dev/null 2>&1
echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NODELIST beta ${BETA:-4} rundir ${RUNDIR:-runs} items $* start $(date -Is)"
for it in "$@"; do
  IFS=: read side L lay nbin var extra <<< "$it"   # optional 5th/6th field (tuning): VARIANT name, comma list of env settings (PINMODE=excl,NDELAY=128,ALFSH_QRNB=32)
  env ${extra//,/ } VARIANT=$var BETA=${BETA:-4} RUNDIR=${RUNDIR:-runs} NBIN=$nbin LAYOUTS_OVERRIDE=${lay/_/:} NT_LOCAL=1 RUN_TIMEOUT=${RUN_TIMEOUT:-9000} \
    srun -N1 -n1 -c 224 --cpu-bind=none --mem=0 --exclusive bash $J/run_layouts.sh $side $L > $J/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.${side}_L${L}_${lay}${var:+_$var}.log 2>&1 &
  sleep 2
done
wait; echo "end $(date -Is)"
