#!/bin/bash
# Ticket 204, nividia: self-test (ALF_GPU_SELFTEST=2, stops after the hook) of binary NAME on a list of configs; one summary file.
#   usage: t204_selftest_all.sh NAME CONFIG [CONFIG ...]     -> logs/selftest_NAME.txt (key lines), runs/NAME/<config>/run.log
W=$HOME/alf-gpu/work/t204
N=$1; shift
out=$W/logs/selftest_$N.txt; : > "$out"
for c in "$@"; do
  r=$(TAG=$N bash $W/t204_run.sh $N $c ALF_GPU_SELFTEST=2 ALF_GPU_EXIT=1)
  echo "== $r" >> "$out"
  grep -E "FAILED|PASSED|structure|GPU wrap self-test (hopping|Op_V|fused|sequences|up/down)" $W/runs/$N/$c/run.log | cut -c1-400 >> "$out"
done
echo done >> "$out"
