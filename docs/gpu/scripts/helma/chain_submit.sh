#!/bin/bash
# usage: chain_submit.sh RUNDIR NGPU RPG MPS_PCT PARTITION WALLTIME NLINKS [MARGIN_H]
# Submits link 1 of NLINKS; every running link submits the next one with --dependency=afterany on itself, so a preempted/killed/timed-out
# link is still followed by the next. All jobs restart from confout/data.h5 in RUNDIR. Prints the job ids. touch RUNDIR/STOP to end the chain softly.
RUNDIR=$1; NGPU=$2; RPG=$3; PCT=$4; PART=$5; WT=$6; N=$7; MARGIN=${8:-0.4}
HERE=/hnvme/workspace/b133ae29-ALF_GPU/ALF/docs/gpu/scripts/helma
SL=/hnvme/workspace/b133ae29-ALF_GPU/slurm; mkdir -p $SL
name=$(basename $RUNDIR)
# only link 1 is submitted; each link submits the next one when it starts (see prod_job.sh)
id=$(sbatch --parsable -D /hnvme/workspace/b133ae29-ALF_GPU/slurm -p $PART --gres=gpu:h100:$NGPU --cpus-per-task=$((32*NGPU)) -t $WT -J $name -o $SL/${name}_%j.out \
      $HERE/prod_job.sh $RUNDIR $NGPU $RPG $PCT $MARGIN $((N-1)) $PART $WT) || exit 1
echo "$name link1 $id $PART gpus=$NGPU $WT (+$((N-1)) further links self-submitted)" | tee -a $RUNDIR/chain.log
