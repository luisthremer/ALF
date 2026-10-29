#!/bin/bash
# Ticket 207, nividia: statistical comparison (compare2.py) of the full device path (updates + stabilization + tau pass) against GPU=none.
#   usage: t207_statrun.sh CONFIG GPUBIN CPUBIN      (runs/stat/<CONFIG>_gpu and _cpu from t207_stat_prep.sh; different seeds)
# The CPU run goes in the background (1 core), the GPU run in the foreground (one GPU process at a time); then compare2.py.
W=$HOME/alf-gpu/work/t207; export TMPDIR=$W/tmp
. ~/alf-gpu/env/nvhpc.sh
C=$1; G=$2; P=$3
case "$C" in ""|*..*|*/*) echo "unsafe name"; exit 1;; esac
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
nvidia-smi --query-gpu=memory.used --format=csv,noheader
( cd "$W/runs/stat/${C}_cpu" && s=$(date +%s.%N) && $W/bin/ALF_$P.out > run.log 2>&1; e=$(date +%s.%N); printf 'cpu %s wall=%.1fs\n' $P $(echo "$e - $s" | bc) ) &
cd "$W/runs/stat/${C}_gpu" || exit 1
s=$(date +%s.%N); timeout ${TLIM:-900} $W/bin/ALF_$G.out > run.log 2>&1; rc=$?; e=$(date +%s.%N)
printf 'gpu %s rc=%s wall=%.1fs\n' $G $rc $(echo "$e - $s" | bc)
wait
cd $W/src_new/testsuite/gpu_reference || exit 1
python3 -I compare2.py "$W/runs/stat/${C}_gpu" "$W/runs/stat/${C}_cpu" --json "$W/logs/compare2_$C.json" | tee "$W/logs/compare2_$C.txt"
echo "compare2 exit code: ${PIPESTATUS[0]}"
