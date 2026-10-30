#!/bin/bash
# Ticket 212 (development only), nividia: run one binary on one config from dev/runs/base (single process, 1 thread).
#   usage: [TAG=dir] t212_run.sh NAME CONFIG [env assignments...]   -> dev/runs/TAG(default NAME)/CONFIG, run.log, rc + wall time
W=$HOME/alf-gpu/work/t212/dev; export TMPDIR=$W/tmp
. ~/alf-gpu/env/nvhpc.sh
N=$1; C=$2; shift 2
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
T=${TAG:-$N}; safe "$N"; safe "$C"; safe "$T"; : "${HOME:?}"
D="$W/runs/$T/$C"; rm -rf -- "$D"; mkdir -p "$W/runs/$T"; cp -r "$W/runs/base/$C" "$D"; cd "$D" || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_TARGET_OFFLOAD=MANDATORY
for kv in "$@"; do export "$kv"; done
s=$(date +%s.%N); $W/bin/ALF_$N.out > run.log 2>&1; rc=$?
e=$(date +%s.%N); printf '%s %s rc=%s wall=%.1fs\n' $N $C $rc $(echo "$e - $s" | bc)
