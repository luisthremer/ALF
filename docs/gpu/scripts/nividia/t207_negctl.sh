#!/bin/bash
# Ticket 207, nividia: negative controls of the tau-pass self-test (each must exit nonzero and say FAILED).
#   usage: t207_negctl.sh NAME [CONFIG]      (default config hub6x6_cb)   -> logs/negctl_NAME.txt
W=$HOME/alf-gpu/work/t207; N=$1; C=${2:-hub6x6_cb}; : > $W/logs/negctl_$N.txt
bad=0
for f in perm prop prec blocks; do
  case $f in prec) want='precision reduction';; *) want='tau';; esac
  TAG=nc_${N}_$f bash $W/t207_run.sh $N $C ALF_GPU_SELFTEST=4 ALF_GPU_EXIT=1 ALF_GPU_FAULT_TAU=$f >> $W/logs/negctl_$N.txt; rc=$?
  # a crash (rc 134/137/139) or a FAILED line about something else is not a detected fault (2026-10-29)
  if [ $rc = 1 ] && grep FAILED $W/runs/nc_${N}_$f/$C/run.log | grep -q "$want"; then echo "fault $f: detected" >> $W/logs/negctl_$N.txt
  else echo "fault $f: NOT DETECTED (rc $rc)" >> $W/logs/negctl_$N.txt; bad=$((bad+1)); fi
  grep -E "FAILED|PASSED|worst slice|precision reduction|CGR2_2 alone, natural" $W/runs/nc_${N}_$f/$C/run.log | sort -u | cut -c1-300 >> $W/logs/negctl_$N.txt
done
cat $W/logs/negctl_$N.txt
exit $bad
