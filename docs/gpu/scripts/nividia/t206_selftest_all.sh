#!/bin/bash
# Ticket 206, nividia: stabilization self-test (ALF_GPU_SELFTEST=3, stops after the hook) of binary NAME on a list of configs.
#   usage: t206_selftest_all.sh NAME CONFIG [CONFIG ...]     -> logs/selftest_NAME.txt (key lines), runs/NAME/<config>/run.log
# Optional environment: NSWEEP-independent; BETA=<value> would need a re-prepared run directory (see t206_prep_beta.sh).
W=$HOME/alf-gpu/work/t206
N=$1; shift
out=$W/logs/selftest_$N.txt; : > "$out"
for c in "$@"; do
  r=$(TAG=$N bash $W/t206_run.sh $N $c ALF_GPU_SELFTEST=3 ALF_GPU_EXIT=1)
  echo "== $r" >> "$out"
  grep -E "FAILED|PASSED|GPU stab" $W/runs/$N/$c/run.log | cut -c1-400 >> "$out"
done
echo done >> "$out"
