#!/bin/bash
# Ticket 204, nividia: timing (ALF_GPU_TIMING=1, stops after the hook) of binary NAME on a list of configs; one summary file.
#   usage: t204_timing_all.sh NAME TAG CONFIG [CONFIG ...]  [extra env via TIMING_ENV="A=1 B=2"]
#   -> logs/timing_NAME_TAG.txt (the "GPU wrap timing" lines), runs/NAME_TAG/<config>/run.log
W=$HOME/alf-gpu/work/t204
N=$1; TG=$2; shift 2
out=$W/logs/timing_${N}_$TG.txt; : > "$out"
for c in "$@"; do
  r=$(TAG=${N}_$TG bash $W/t204_run.sh $N $c ALF_GPU_TIMING=1 ALF_GPU_EXIT=1 $TIMING_ENV)
  echo "== $r" >> "$out"
  grep -E "GPU wrap timing|FAILED|rror" $W/runs/${N}_$TG/$c/run.log | cut -c1-300 >> "$out"
done
echo done >> "$out"
