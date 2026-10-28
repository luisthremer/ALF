#!/bin/bash
# Ticket 207, nividia: GPU precision runs of the final code (tags fin_def = nv_def, fin_log = nv_log) on prepared configs, serially
# (one GPU process at a time); the CPU runs (nv_none, nv_none_pp, nv_none_log of the same config) must exist from t207_precrun.sh.
#   usage: t207_gpuprec.sh MODES CONFIG [CONFIG ...]     MODES = def | log | both     -> logs/gpuprec.txt
W=$HOME/alf-gpu/work/t207
M=$1; shift
for C in "$@"; do
  for m in def log; do
    case "$M" in both|$m) ;; *) continue;; esac
    nvidia-smi --query-gpu=memory.used --format=csv,noheader
    TAG=fin_$m timeout 900 bash $W/t207_run.sh nv_$m $C
  done
  echo "== $C" >> $W/logs/gpuprec.txt
  bash $W/t207_prectab.sh $C fin_def fin_log nv_none nv_none_pp nv_none_log 2>/dev/null | tee -a $W/logs/gpuprec.txt
done
