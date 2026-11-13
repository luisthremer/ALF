#!/bin/bash
# mean over ranks of the TEAM_TIMERS of the 3-sweep runs (s/sweep): sweep, obs; and the CPU-time of rank 0
cd $SCRATCH/large_prompt/team/p8/obs/work_g4/runs
for t in $(ls | grep -E "_s3$|^A_" | sort); do
  grep -h TEAM_TIMERS $t/run.log | awk -v t=$t '{n=$(NF-7); s+=$(NF-6)/n; o+=$NF/n; u+=($(NF-5)+$(NF-4))/n; st+=($(NF-3)+$(NF-2))/n; c+=$(NF-1)/n; k++} END{printf "%-14s ranks %2d  sweep %6.3f  upd %5.3f  stab %5.3f  cgr %5.3f  obs %5.3f\n", t, k, s/k, u/k, st/k, c/k, o/k}'
done
