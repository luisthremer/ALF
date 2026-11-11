#!/bin/bash
# usage: chain_submit.sh RUNDIR NGPU RPG MPS_PCT PARTITION WALLTIME NLINKS [MARGIN_H]
# Submits NLINKS links of prod_job.sh for RUNDIR, each depending (afterany) on the previous one, so a preempted/killed/timed-out link is
# followed by the next. All jobs restart from confout/data.h5 in RUNDIR. Prints the job ids. touch RUNDIR/STOP to end the chain softly.
RUNDIR=$1; NGPU=$2; RPG=$3; PCT=$4; PART=$5; WT=$6; N=$7; MARGIN=${8:-0.4}
HERE=/hnvme/workspace/b133ae29-ALF_GPU/ALF/docs/gpu/scripts/helma
SL=/hnvme/workspace/b133ae29-ALF_GPU/slurm; mkdir -p $SL
name=$(basename $RUNDIR); prev=""
for i in $(seq 1 $N); do
  dep=""; [ -n "$prev" ] && dep="--dependency=afterany:$prev"
  id=$(sbatch --parsable -p $PART --gres=gpu:h100:$NGPU --cpus-per-task=$((32*NGPU)) -t $WT -J $name $dep -o $SL/${name}_%j.out \
        $HERE/prod_job.sh $RUNDIR $NGPU $RPG $PCT $MARGIN) || exit 1
  echo "$name link$i $id $PART gpus=$NGPU $WT" | tee -a $RUNDIR/chain.log
  prev=$id
done
