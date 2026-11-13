#!/bin/bash
# node-hours of all p8 jobs from the start/end lines of their logs (cluster side)
cd $SCRATCH/large_prompt/team/p8/obs/jobs
for f in p8_*.log; do
  s=$(grep -m1 "^start " $f | awk '{print $2}'); e=$(grep "^end " $f | tail -1 | awk '{print $2}')
  [ -z "$s" ] && continue; [ -z "$e" ] && e=$(date -Is -r $f)
  a=$(date -d "$s" +%s); b=$(date -d "$e" +%s); echo "$f $((b-a))"
done | awk '{t+=$2; printf "%s %d s\n",$1,$2} END{printf "TOTAL %.1f min = %.3f node-h\n", t/60, t/3600}'
