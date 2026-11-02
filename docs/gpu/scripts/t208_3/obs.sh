#!/bin/bash
# 207 close, check 4: one production-size GPU run (hub12x12_cb, beta 4, 250 bins x 70 sweeps, Ltau 1, seeds rotated by 37 lines so the
# stream differs from the CPU reference) per path, then compare2.py against the CPU reference run (16-rank GNU, same parameters).
#   usage: obs.sh P(nv|in) REFROOT    REFROOT/{ref_gnu-DEFAULT,ref_gnu-LOG}/hub12x12_cb = CPU reference dirs (data.h5, parameters, seeds)
# -> $W/runs/obs_<path>/ ; summary lines appended to $W/logs/obs.txt
W=$HOME/alf-gpu/work/t208-3; P=$1; RR=$2; export TMPDIR=$W/tmp OMP_NUM_THREADS=1
if [ $P = nv ]; then . ~/alf-gpu/env/nvhpc.sh; else . ~/alf-gpu/env/intel.sh 2>/dev/null; export ONEAPI_DEVICE_SELECTOR=level_zero:0; fi
PY=~/alf-gpu/venv/bin/python3
for k in def log; do
  K=$([ $k = def ] && echo DEFAULT || echo LOG); D=$W/runs/obs_$k; R=$RR/ref_gnu-$K/hub12x12_cb
  if [ ! -e $D/.rc ]; then
    mkdir -p $D && cp $R/parameters $D/ && ( tail -n +38 $R/seeds; head -n 37 $R/seeds ) > $D/seeds && cd $D
    s=$(date +%s); OMP_TARGET_OFFLOAD=MANDATORY flock $HOME/alf-gpu/gpu.lock timeout 14400 $W/bin/ALF_${P}_$k.out > run.log 2>&1; echo $? > .rc
    echo "run $P $k rc=$(cat .rc) wall=$(( $(date +%s) - s ))s" >> $W/logs/obs.txt; cd $W
  fi
  echo "== compare2 $P $k vs CPU ref $K" >> $W/logs/obs.txt
  ( cd $W/src/testsuite/gpu_reference && $PY -I compare2.py $R $D --json $D/compare2.json ) >> $W/logs/obs.txt 2>&1; echo "compare2 rc=$?" >> $W/logs/obs.txt
done
echo OBS_DONE >> $W/logs/obs.txt
