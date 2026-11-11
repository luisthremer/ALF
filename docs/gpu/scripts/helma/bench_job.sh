#!/bin/bash -l
# Benchmark on one node: sbatch [-p preempt|h100] --gres=gpu:h100:NGPU --cpus-per-task=32*NGPU -t HH:MM:SS bench_job.sh NGPU "L:RPG:PCT:ND ..." [NSWEEP_PER_BIN NBIN]
# Each config runs NGPU parameter groups (U = 5.3 5.35 5.4 5.45), RPG walkers per group, group g on GPU g (all ranks of a group share
# one GPU), MPS thread share PCT (0 = off), N_delay ND, L, NSweep=2 x NBin=3 = 6 sweeps unless given.
#SBATCH -A b133ae
#SBATCH --export=NONE
#SBATCH -J alfbench
unset SLURM_EXPORT_ENV
export RUN_TIMEOUT=${RUN_TIMEOUT:-420}   # s per config (no-MPS runs with many ranks are >10x slower)
NGPU=$1; CFGS=$2; NS=${3:-2}; NB=${4:-3}
WS=/hnvme/workspace/b133ae29-ALF_GPU; S=$WS/ALF/docs/gpu/scripts/helma
OUT=$WS/bench/$SLURM_JOB_ID; mkdir -p $OUT
UL=5.3,5.35,5.4,5.45; UL=$(echo $UL | cut -d, -f1-$NGPU)
for cfg in $CFGS; do
  IFS=: read L RPG PCT ND <<< "$cfg"
  tag=L${L}_g${NGPU}_r${RPG}_m${PCT}_d${ND}
  python3 $S/setup_pp.py $OUT/$tag $L $NS $ND 0 $RPG $UL --nbin $NB --master bench >/dev/null || continue
  $S/run_pp.sh $OUT/$tag $NGPU $RPG $PCT > /dev/null 2>&1
  echo "done $tag rc=$(cat $OUT/$tag/DONE)"
  python3 $S/bench_summarize.py $OUT/$tag $NGPU $RPG $PCT $ND $L $((NS*NB)) | tee -a $OUT/results.csv
done
