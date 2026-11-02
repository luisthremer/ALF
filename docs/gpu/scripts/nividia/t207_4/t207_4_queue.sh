#!/bin/bash
# Ticket 207-4 (development only; copy of t207_3_queue.sh), nividia: run chain steps one after another, after the chain currently running (if any) has ended.
#   usage: t207_4_queue.sh NAME 'PHASE ARGS' ['PHASE ARGS' ...]   -> logs/queue_<NAME>.txt
#   e.g.   t207_4_queue.sh q2 'main dev11' 'ctl dev14'
W=$HOME/alf-gpu/work/t207-4; cd $W || exit 1
N=$1; shift; out=$W/logs/queue_$N.txt
echo "queue $N start $(date -Is), waiting for running chains" >> $out
while pgrep -u "$USER" -f "t207_4_chain.sh" > /dev/null; do sleep 30; done
for step in "$@"; do
  echo "step [$step] $(date -Is)" >> $out
  bash $W/t207_4_chain.sh $step
  echo "step [$step] ended rc=$? $(date -Is)" >> $out
done
echo "queue $N done $(date -Is)" >> $out
