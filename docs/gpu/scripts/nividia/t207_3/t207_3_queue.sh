#!/bin/bash
# Ticket 207-3 (development only), nividia: run chain steps one after another, after the chain currently running (if any) has ended.
#   usage: t207_3_queue.sh NAME 'PHASE ARGS' ['PHASE ARGS' ...]   -> logs/queue_<NAME>.txt
#   e.g.   t207_3_queue.sh q2 'main dev2 dev3 dev4' 'ctl dev5'
W=${T207_3W:-$HOME/alf-gpu/work/t207-3};   # T207_3W: work dir of a follow-up agent (t207-3b)
cd $W || exit 1
N=$1; shift; out=$W/logs/queue_$N.txt
echo "queue $N start $(date -Is), waiting for running chains" >> $out
# fix round t207-3c: the chain runs as "bash $W/t207_3_chain.sh" (full path); the old pattern "bash t207_3_chain.sh" never matched it
while pgrep -u "$USER" -f "bash $W/t207_3_chain.sh" > /dev/null; do sleep 30; done
for step in "$@"; do
  echo "step [$step] $(date -Is)" >> $out
  bash $W/t207_3_chain.sh $step
  echo "step [$step] ended rc=$? $(date -Is)" >> $out
done
echo "queue $N done $(date -Is)" >> $out
