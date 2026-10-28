#!/bin/bash
# Ticket 205, nividia: one-slice update test (ALF_GPU_UPDTEST=1, stops after the test) of binary NAME on a list of configs.
#   usage: t205_updtest_all.sh NAME CONFIG [CONFIG ...]     -> logs/updtest_NAME.txt
W=$HOME/alf-gpu/work/t205
N=$1; shift
out=$W/logs/updtest_$N.txt; : > "$out"
for c in "$@"; do
  r=$(TAG=ut_$N bash $W/t205_run.sh $N $c ALF_GPU_UPDTEST=1 ALF_GPU_EXIT=1)
  echo "== $r" >> "$out"
  grep -E "FAILED|PASSED|GPU update (self-test|path)" $W/runs/ut_$N/$c/run.log | sort -u | cut -c1-400 >> "$out"
done
echo done >> "$out"
