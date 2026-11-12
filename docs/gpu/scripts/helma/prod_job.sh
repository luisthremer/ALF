#!/bin/bash -l
# One production chain link. Submit through chain_submit.sh (sets partition, gres, time, dependency).
# usage: prod_job.sh RUNDIR NGPU RANKS_PER_GPU MPS_PCT MARGIN_HOURS REMAIN PARTITION WALLTIME
# REMAIN = number of further links; when this link starts it submits the next one (afterany on itself), so only one link per chain is ever
# pending (Slurm backfill here considers only bf_max_job_user=20 jobs per user).
# CPU_MAX of this link = remaining walltime - margin (margin covers the bin in progress, 10 sweeps; ALF stops at the next bin end after CPU_MAX).
#SBATCH -A b133ae
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
RUNDIR=$1; NGPU=$2; RPG=$3; PCT=$4; MARGIN=${5:-0.4}; REMAIN=${6:-0}; PART=${7:-h100}; WT=${8:-24:00:00}
HERE=/hnvme/workspace/b133ae29-ALF_GPU/ALF/docs/gpu/scripts/helma
[ -e "$RUNDIR/STOP" ] && { echo "STOP file present"; exit 0; }
if [ "$REMAIN" -gt 0 ]; then
  nid=$(sbatch --parsable -D /hnvme/workspace/b133ae29-ALF_GPU/slurm -p $PART --gres=gpu:h100:$NGPU --cpus-per-task=$((32*NGPU)) -t $WT -J $(basename $RUNDIR) \
        --dependency=afterany:$SLURM_JOB_ID -o /hnvme/workspace/b133ae29-ALF_GPU/slurm/$(basename $RUNDIR)_%j.out \
        $HERE/prod_job.sh $RUNDIR $NGPU $RPG $PCT $MARGIN $((REMAIN-1)) $PART $WT)
  echo "$(basename $RUNDIR) next link $nid (afterany:$SLURM_JOB_ID) remaining after it $((REMAIN-1))" >> $RUNDIR/chain.log
fi
rem=$(squeue -h -j $SLURM_JOB_ID -o %L)          # [D-]HH:MM:SS
sec=$(python3 - "$rem" <<'P'
import sys
s=sys.argv[1]; d=0
if "-" in s: d,s=s.split("-"); d=int(d)
p=[int(x) for x in s.split(":")]; p=[0]*(3-len(p))+p
print(d*86400+p[0]*3600+p[1]*60+p[2])
P
)
cpumax=$(python3 -c "print('%.3f' % max(0.05, $sec/3600.0 - $MARGIN))")
echo "link $SLURM_JOB_ID on $(hostname) start $(date) remaining ${sec}s CPU_MAX=$cpumax h" >> $RUNDIR/links.log
$HERE/run_pp.sh $RUNDIR $NGPU $RPG $PCT $cpumax
rc=$?
echo "link $SLURM_JOB_ID end $(date) rc=$rc" >> $RUNDIR/links.log
exit $rc
