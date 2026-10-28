#!/bin/bash
# Ticket 207, nividia: negative controls of the tau-pass self-test (each must exit nonzero and say FAILED).
#   usage: t207_negctl.sh NAME [CONFIG]      (default config hub6x6_cb)   -> logs/negctl_NAME.txt
W=$HOME/alf-gpu/work/t207; N=$1; C=${2:-hub6x6_cb}; : > $W/logs/negctl_$N.txt
for f in perm prop prec blocks; do
  TAG=nc_${N}_$f bash $W/t207_run.sh $N $C ALF_GPU_SELFTEST=4 ALF_GPU_EXIT=1 ALF_GPU_FAULT_TAU=$f >> $W/logs/negctl_$N.txt
  grep -E "FAILED|PASSED|worst slice|precision reduction|CGR2_2 alone, natural" $W/runs/nc_${N}_$f/$C/run.log | sort -u | cut -c1-300 >> $W/logs/negctl_$N.txt
done
cat $W/logs/negctl_$N.txt
