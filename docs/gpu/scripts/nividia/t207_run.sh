#!/bin/bash
# Ticket 207, nividia: run one binary on one config from runs/base.
#   usage: [TAG=dir] t207_run.sh NAME CONFIG [env assignments...]   -> runs/TAG(default NAME)/CONFIG (copy of base/CONFIG), log run.log, wall time
# NAME = binary suffix (bin/ALF_<NAME>.out).  Single process, 1 thread, GPU used only while this runs.
W=$HOME/alf-gpu/work/t207; export TMPDIR=$W/tmp
N=$1; C=$2; shift 2
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
T=${TAG:-$N}; safe "$N"; safe "$C"; safe "$T"; : "${HOME:?}"
D="$W/runs/$T/$C"; rm -rf -- "$D"; mkdir -p "$W/runs/$T"; cp -r "$W/runs/base/$C" "$D"; cd "$D" || exit 1
# env -i (2026-10-29): only the env file and the listed variables reach the program, nothing inherited from the caller
s=$(date +%s.%N)
env -i HOME="$HOME" USER="$USER" LANG=C PATH=/usr/bin:/bin "$@" bash --noprofile --norc -c '. "$HOME/alf-gpu/env/nvhpc.sh" > /dev/null 2>&1; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; env | sort > env.txt; exec "$0" > run.log 2>&1' "$W/bin/ALF_$N.out"; rc=$?
e=$(date +%s.%N); printf '%s %s rc=%s wall=%.1fs\n' $N $C $rc $(echo "$e - $s" | bc)
exit $rc
