#!/bin/bash
# 208-3, all remaining checks in one go (runs on the remote machine, detached).  usage: final.sh nv|in
# waits for an earlier t208-3 chain, then: builds (GPU + GPU=none) -> T0 (self-test, same-seed) -> production observables (def, LOG)
# -> [nividia only] GPU=none ctest -> logs/SUMMARY.txt (ends with FINAL_DONE).  Work dir ~/alf-gpu/work/t208-3f.
P=$1; W=$HOME/alf-gpu/work/t208-3f; cd $W || exit 1; S=$W/logs/SUMMARY.txt
while pgrep -u $USER -f '[s]cripts/t0.sh|[s]cripts/builds_|[s]cripts/obs.sh|[s]cripts/ctest.sh|[t]208-3/bin/ALF' > /dev/null; do sleep 30; done
if [ $P = nv ]; then
  bash scripts/builds_nvidia.sh > logs/builds.txt 2>&1; IR=$HOME/alf-gpu/work/t207-1/runs/base
else
  bash scripts/builds_intel.sh > logs/builds.txt 2>&1; IR=$HOME/alf-gpu/work/t302c-2b/runs/base
fi
bash scripts/t0.sh $P $IR > logs/t0.out 2>&1
bash scripts/obs.sh $P $W/ref > logs/obs.out 2>&1
[ $P = nv ] && bash scripts/ctest.sh > logs/ctest.out 2>&1
{ echo "== machine $P  $(hostname)  $(date -Is)"
  echo "== builds";           grep -E '^built|failed' logs/builds.txt
  echo "== T0 (self-test, same-seed)"; grep -E '^st_|^same|T0_DONE' logs/t0.txt | cut -c1-200
  echo "== observables";      grep -E '^run|sum t|RESULT|OBS_DONE' logs/obs.txt | cut -c1-120
  echo "== transfers (obs_def)"; grep -ahE 'INTERIM 205|interim bytes' runs/obs_def/run.log | cut -c1-300 | tail -3
  echo "== ctest";            grep -hE 'tests passed' logs/ctest_*.txt 2>/dev/null
  echo FINAL_DONE; } > $S 2>&1
