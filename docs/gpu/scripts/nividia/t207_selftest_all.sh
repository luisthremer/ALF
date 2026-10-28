#!/bin/bash
# Ticket 207, nividia: tau-pass self-test (ALF_GPU_SELFTEST=4, stops after the hook) of binary NAME on a list of configs.
#   usage: t207_selftest_all.sh NAME CONFIG [CONFIG ...]     -> logs/selftest_NAME.txt (key lines), runs/NAME/<config>/run.log
# Extra environment assignments can be passed through EXTRA (e.g. EXTRA='ALF_GPU_FAULT_TAU=perm').
W=$HOME/alf-gpu/work/t207
N=$1; shift
out=$W/logs/selftest_$N.txt; : > "$out"
for c in "$@"; do
  r=$(TAG=$N bash $W/t207_run.sh $N $c ALF_GPU_SELFTEST=4 ALF_GPU_EXIT=1 $EXTRA)
  echo "== $r" >> "$out"
  grep -E "FAILED|PASSED|GPU tau|worst slice" $W/runs/$N/$c/run.log | cut -c1-400 >> "$out"
done
echo done >> "$out"
