#!/bin/bash
# Ticket 302c-2 round 3 (development only), on intel: the QREF run sequence, unattended (decision 36), one GPU run at a time, each
# processed before the next (t302c2c_chain.sh).  -> logs/queue.txt
#   1. bit identity: reruns of round 2's m4_def_dev5 and m10_log_dev5 (same seeds, env) with the d3e82a36 binaries, compared with round 2
#      (t302c2c_bitid.sh -> logs/bitid.txt).  The queue stops unless both raw-dump manifests and both qref.txt (tq masked) are identical.
#   2. main beta = 10, dev10 and dev11, both paths (with round 2's dev5, dev6: 4 sets per path)
#   3. controls on fresh seed sets (dev12..dev18; rounds 1 and 2 used dev1..dev9; checked by seeds sha256 before the round):
#      beta 4:  dev12 def eps 4e-15 | dev13 def 4e-15 AT 0:2 | dev14 LOG 1e-14 AT 0:2
#      beta 10: dev15 eps 3e-12, 1e-11 (both paths) | dev16 localized 2e-11, 1e-10 AT 2:2 (both paths) | dev17 LOG 3e-14 |
#               dev18 LOG 3e-14 AT 0:2, def 1e-14 AT 0:2
#   usage: t302c2c_queue.sh [STEP...]   (default: bitid main ctl)
cd ~/alf-gpu/work/t302c-2c || exit 1; export TMPDIR=$PWD/tmp
q() { echo "phase $* start $(date -Is)" >> logs/queue.txt; bash t302c2c_chain.sh "$@"; local rc=$?; echo "phase $* end rc=$rc $(date -Is)" >> logs/queue.txt
      [ $rc = 0 ] || { echo "QUEUE STOPPED (rc $rc) $(date -Is)" >> logs/queue.txt; exit 1; }; }
for step in ${@:-bitid main ctl}; do case $step in
bitid)
  PATHS=def q main 4 dev5
  PATHS=log q main 10 dev5
  bash t302c2c_bitid.sh m4_def_dev5 m10_log_dev5 > logs/bitid.txt 2>&1
  ok=$(grep -cE 'qref_manifest.txt .*: identical|qref.txt with tq masked .*: IDENTICAL' logs/bitid.txt)
  echo "bitid: $ok of 4 checks (2 manifests, 2 qref.txt) identical $(date -Is)" >> logs/queue.txt
  [ "$ok" = 4 ] || { echo "QUEUE STOPPED: bit identity not shown (see logs/bitid.txt) $(date -Is)" >> logs/queue.txt; exit 1; } ;;
main)
  q main 10 dev10 dev11 ;;
ctl)
  q flt 4 def dev12 4e-15
  q fltat 4 def dev13 0:2 4e-15
  q fltat 4 log dev14 0:2 1e-14
  for p in def log; do q flt 10 $p dev15 3e-12 1e-11; done
  for p in def log; do q fltat 10 $p dev16 2:2 2e-11 1e-10; done
  q flt 10 log dev17 3e-14
  q fltat 10 log dev18 0:2 3e-14
  q fltat 10 def dev18 0:2 1e-14 ;;
*) echo "bad step $step" >> logs/queue.txt; exit 1 ;;
esac; done
echo "QUEUE DONE $(date -Is)" >> logs/queue.txt
