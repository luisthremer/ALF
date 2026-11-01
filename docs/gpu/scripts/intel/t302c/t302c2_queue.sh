#!/bin/bash
# Ticket 302c-2 (development only), on intel: the whole run sequence, unattended (decision 36), one run at a time, each processed before
# the next (t302c2_chain.sh).  -> logs/queue.txt (start/end of each phase) and logs/chain_<phase>.txt
cd ~/alf-gpu/work/t302c-2 || exit 1
q() { echo "phase $* start $(date -Is)" >> logs/queue.txt; bash t302c2_chain.sh "$@"; echo "phase $* end $(date -Is)" >> logs/queue.txt; }
q main 4 dev1
q ctl 4 dev1
q main 4 dev2 dev3 dev4
q main 10 dev1 dev2
q ctl 10 dev1
echo "QUEUE DONE $(date -Is)" >> logs/queue.txt
