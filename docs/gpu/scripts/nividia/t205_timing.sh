#!/bin/bash
# Ticket 205, nividia: s/sweep GPU (device update, interim copies) vs the same GPU build on the host path
# (ALF_GPU_NOUPDATE=1) vs the GPU=none build, for N = 36, 144, 400, with 8 sweeps (NSweep 4 x NBin 2) per run.
#   Run only while `nvidia-smi` shows no other process (the A30 is shared with ticket 206).
#   usage: t205_timing.sh [CONFIG ...]     (default hub6x6_cb hub12x12_cb hub20x20_cb) -> logs/timing.txt
W=$HOME/alf-gpu/work/t205; export TMPDIR=$W/tmp; . ~/alf-gpu/env/gnu.sh
cfgs=${@:-hub6x6_cb hub12x12_cb hub20x20_cb}
cd $W/src_new/testsuite/gpu_reference || exit 1
: "${HOME:?}"
rm -rf -- "$W/runs/time"; mkdir -p "$W/runs"
python3 prep.py $W/runs/time $cfgs --nsweep 4 --nbin 2 --alf-dir $W/src_new > /dev/null
echo "A30 processes before the runs:" > $W/logs/timing.txt; nvidia-smi --query-compute-apps=pid,name --format=csv,noheader >> $W/logs/timing.txt
for c in $cfgs; do
  for v in none host gpu; do
    case $v in
      none) bin=nv_none; envs="";;
      host) bin=nv_def;  envs="ALF_GPU_NOUPDATE=1";;
      gpu)  bin=nv_def;  envs="";;
    esac
    BASE=time TAG=tm_$v bash $W/t205_run.sh $bin $c $envs >> $W/logs/timing.txt
    echo "   $v $c: $(grep 'CPU Time' $W/runs/tm_$v/$c/info | tr -s ' ')" >> $W/logs/timing.txt
  done
done
cat $W/logs/timing.txt
