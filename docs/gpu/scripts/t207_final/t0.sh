#!/bin/bash
# 207 close (decision 41), T0: self-test (ALF_GPU_SELFTEST=4) and same-seed GPU vs GPU=none, default and LOG.
#   usage: t0.sh P(nv|in) INPUTROOT      P: binary prefix (nv_def, nvn_def | in_def, inn_def); INPUTROOT holds hub6x6_cb, hub6x6_cb_fd
# -> $W/runs/t0_<id>/, summary $W/logs/t0.txt.  GPU runs inside flock ~/alf-gpu/gpu.lock.
W=$HOME/alf-gpu/work/t207-final; P=$1; IR=$2; export TMPDIR=$W/tmp OMP_NUM_THREADS=1
if [ $P = nv ]; then . ~/alf-gpu/env/nvhpc.sh; else . ~/alf-gpu/env/intel.sh 2>/dev/null; export ONEAPI_DEVICE_SELECTOR=level_zero:0; fi
PY=~/alf-gpu/venv/bin/python3
run() { # id bin input gpu(1/0) env...
  id=$1; b=$2; c=$3; g=$4; shift 4; D=$W/runs/t0_$id; [ -e $D ] && return
  mkdir -p $D && cp -r $IR/$c/. $D/ && cd $D && mkdir -p qref
  [ "$c" = hub6x6_cb_fd ] && sed -i -E 's/^(NSweeps *= *)[0-9]+/\11/; s/^(NBins *= *)[0-9]+/\160/' parameters
  L=(); [ $g = 1 ] && L=(flock $HOME/alf-gpu/gpu.lock)
  s=$(date +%s); env "$@" "${L[@]}" timeout 3600 $W/bin/ALF_$b.out > run.log 2>&1; echo $? > .rc
  echo "$id rc=$(cat .rc) wall=$(( $(date +%s) - s ))s $(grep -hoE 'GPU tau self-test (PASSED|FAILED)' run.log | sort -u | tr '\n' ' ') $*" >> $W/logs/t0.txt; cd $W
}
for k in def log; do
  run st_$k ${P}_$k$([ $k = def ] && echo _pp) hub6x6_cb 1 ALF_GPU_DEBUG=1 ALF_GPU_EXIT=1 ALF_GPU_SELFTEST=4
  run sd_gpu_$k ${P}_$k hub6x6_cb_fd 1 ALF_NONE=1
  run sd_none_$k ${P}n_$k hub6x6_cb_fd 0 ALF_NONE=1
  echo "same-seed $k: $($PY -I $W/scripts/t207_firstdiv.py $W/runs/t0_sd_gpu_$k/data.h5 $W/runs/t0_sd_none_$k/data.h5 2>&1 | tail -4 | tr '\n' ' ')" >> $W/logs/t0.txt
done
echo T0_DONE >> $W/logs/t0.txt
