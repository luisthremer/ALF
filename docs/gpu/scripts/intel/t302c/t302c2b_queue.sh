#!/bin/bash
# Ticket 302c-2 round 2 (development only), on intel: the QREF run sequence, unattended (decision 36), one GPU run at a time, each
# processed before the next (t302c2b_chain.sh).  -> logs/queue.txt
#   seeds (committed dev seeds, rotation 37*N): main dev5, dev6 (fresh, not used in round 1); faults dev7 (near: 2e-15, 1e-14),
#   dev8 (far: 1e-12, 1e-10), dev9 (localized 2e-12 at (2,2)); round 1 used dev1..dev4 (beta 4) and dev1, dev2 (beta 10)
#   usage: t302c2b_queue.sh [4|10 ...]   (default: 4 10)
cd ~/alf-gpu/work/t302c-2b || exit 1
q() { echo "phase $* start $(date -Is)" >> logs/queue.txt; bash t302c2b_chain.sh "$@"; echo "phase $* end rc=$? $(date -Is)" >> logs/queue.txt; }
for b in ${@:-4 10}; do
  q main $b dev5 dev6
  for p in def log; do
    q flt $b $p dev7 2e-15 1e-14
    q flt $b $p dev8 1e-12 1e-10
    q fltat $b $p dev9 2:2 2e-12
  done
done
echo "QUEUE DONE $(date -Is)" >> logs/queue.txt
