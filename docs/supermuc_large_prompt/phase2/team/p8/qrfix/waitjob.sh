#!/bin/bash
# usage: waitjob.sh JOBID [maxseconds]  -> returns when the job has left the queue
for i in $(seq 1 $(( ${2:-1800} / 20 ))); do
  o=$(timeout 60 ssh supermuc2 "squeue -h -j $1 -o %T" 2>/dev/null | grep -E "RUNNING|PENDING|COMPLETING|CONFIGURING")
  [ -z "$o" ] && { echo "job $1 finished"; exit 0; }
  sleep 20
done; echo "timeout waiting"
