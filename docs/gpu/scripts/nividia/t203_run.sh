#!/bin/bash
# Ticket 203, nividia: run one binary on one config from runs/base.
#   usage: [TAG=dir] t203_run.sh NAME CONFIG [env assignments...]   -> runs/TAG(default NAME)/CONFIG (copy of base/CONFIG), log run.log, wall time
# NAME = binary suffix (bin/ALF_<NAME>.out).  Single process, 1 thread, GPU used only while this runs.
W=$HOME/alf-gpu/work/t203; export TMPDIR=$W/tmp
. ~/alf-gpu/env/nvhpc.sh
N=$1; C=$2; shift 2
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
T=${TAG:-$N}; safe "$N"; safe "$C"; safe "$T"; : "${HOME:?}"
D="$W/runs/$T/$C"; rm -rf -- "$D"; mkdir -p "$W/runs/$T"; cp -r "$W/runs/base/$C" "$D"; cd "$D" || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for kv in "$@"; do export "$kv"; done
s=$(date +%s.%N); $W/bin/ALF_$N.out > run.log 2>&1; rc=$?
e=$(date +%s.%N); printf '%s %s rc=%s wall=%.1fs\n' $N $C $rc $(echo "$e - $s" | bc)
