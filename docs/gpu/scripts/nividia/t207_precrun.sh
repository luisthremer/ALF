#!/bin/bash
# Ticket 207, nividia: precision runs of one config (identical seeds), GPU runs strictly one after the other, CPU runs in the
# background (<= 3 at a time).  usage: t207_precrun.sh CONFIG   -> runs/<bin>/CONFIG, table in logs/prectab_CONFIG.txt
#   GPU: nv_def (default), nv_log (LOG);  CPU: nv_none (GPU=none), nv_none_pp (QR_PREPIVOT), nv_none_log (CPU LOG)
W=$HOME/alf-gpu/work/t207; C=$1
for b in nv_none nv_none_pp nv_none_log; do TAG=$b bash $W/t207_run.sh $b $C & done
for b in nv_def nv_log; do
  nvidia-smi --query-gpu=memory.used --format=csv,noheader
  TAG=$b timeout 900 bash $W/t207_run.sh $b $C
done
wait
bash $W/t207_prectab.sh $C nv_def nv_log nv_none nv_none_pp nv_none_log | tee $W/logs/prectab_$C.txt
