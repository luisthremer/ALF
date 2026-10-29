#!/bin/bash
# Ticket 207 (LOG tail follow-up), nividia: precision runs of the final GPU binaries, strictly one GPU process at a time.
#   usage: t207_finprec.sh BIN:CONFIG [BIN:CONFIG ...]     -> runs/fin_<BIN>/<CONFIG>
W=$HOME/alf-gpu/work/t207
for x in "$@"; do
  b=${x%%:*}; c=${x#*:}
  nvidia-smi --query-gpu=memory.used --format=csv,noheader
  TAG=fin_$b timeout 900 bash $W/t207_run.sh "$b" "$c" || { echo "FAIL: run fin_$b $c"; exit 1; }
  echo "== $c" >> $W/logs/finprec.txt
  bash $W/t207_prectab.sh $c fin_$b | tee -a $W/logs/finprec.txt || exit 1   # generated, never hand-assembled (2026-10-29)
done
