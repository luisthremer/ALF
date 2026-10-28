#!/bin/bash
# Ticket 207, nividia: like t207_precrun.sh, but only the LOG pair (GPU nv_log, CPU nv_none_log) for several seed sets, to see
# the spread of the "max" statistic.  usage: t207_precrun2.sh CONFIG [CONFIG ...]    (each a prepared run directory, see t207_seeds.sh)
#   GPU runs one after the other; the CPU LOG run of a config goes in the background.  Table: logs/prectab2.txt
W=$HOME/alf-gpu/work/t207
: > $W/logs/prectab2.txt
for C in "$@"; do
  TAG=nv_none_log bash $W/t207_run.sh nv_none_log $C &
  nvidia-smi --query-gpu=memory.used --format=csv,noheader
  TAG=nv_log timeout 900 bash $W/t207_run.sh nv_log $C
  wait
  echo "== $C" >> $W/logs/prectab2.txt
  bash $W/t207_prectab.sh $C nv_log nv_none_log | tee -a $W/logs/prectab2.txt
done
