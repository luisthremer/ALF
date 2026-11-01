#!/bin/bash
# Ticket 207-1 (development only), on intel: run one binary on one prepared config (1 thread).  usage: [SET=base] [TAG=dir] t207_1_irun.sh NAME CONFIG [env..]
# -> runs/TAG/CONFIG (refuses to overwrite).  GPU builds (names starting in_ but not in_none*) get the mandatory B580 environment.
W=$HOME/alf-gpu/work/t207-1; export TMPDIR=$W/tmp
N=$1; C=$2; shift 2
. ~/alf-gpu/env/intel.sh 2>/dev/null
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
T=${TAG:-$N}; S=${SET:-base}; safe "$N"; safe "$C"; safe "$T"; safe "$S"; : "${HOME:?}"
D="$W/runs/$T/$C"; [ -e "$D" ] && { echo "refusing: $D exists"; exit 1; }
mkdir -p "$W/runs/$T"; cp -r "$W/runs/$S/$C" "$D"; cd "$D" || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for kv in "$@"; do export "$kv"; done
case $N in in_none*) ;; *) export ONEAPI_DEVICE_SELECTOR=level_zero:0 OMP_TARGET_OFFLOAD=MANDATORY FORT_BUFFERED=false;; esac
s=$(date +%s.%N); timeout -k 30 3600 $W/bin/ALF_$N.out > run.log 2>&1; rc=$?; e=$(date +%s.%N)
printf '%s %s %s rc=%s wall=%.1fs bin %s env [%s]\n' "$T" $N $C $rc $(echo "$e - $s" | bc) "$(cut -c1-16 $W/bin/ALF_$N.out.sha256)" "$*"
