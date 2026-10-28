#!/bin/bash
# Ticket 205, nividia: negative controls of the update test (must exit nonzero and say FAILED).  usage: t205_negctl.sh NAME
W=$HOME/alf-gpu/work/t205; N=$1; : > $W/logs/negctl_$N.txt
for f in upd nsig; do
  TAG=nc_${N}_$f bash $W/t205_run.sh $N hub6x6_cb ALF_GPU_UPDTEST=1 ALF_GPU_EXIT=1 ALF_GPU_FAULT=$f >> $W/logs/negctl_$N.txt
  grep -E "FAILED|PASSED|trial" $W/runs/nc_${N}_$f/hub6x6_cb/run.log | sort -u | cut -c1-330 >> $W/logs/negctl_$N.txt
done
