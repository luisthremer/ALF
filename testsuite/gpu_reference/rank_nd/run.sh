#!/bin/bash
# Ticket 401 (development): run one binary on one prepared config (single process, 1 thread), on nividia or intel.
#   usage: [TAG=dir] run.sh BIN CONFIG [env assignments...]   -> runs/TAG(default BIN)/CONFIG with run.log; prints rc and wall time
# GPU runs on nividia take ~/alf-gpu/gpu.lock (flock) when GPULOCK=1.  On intel: B580 only (level_zero:0).
W=$HOME/alf-gpu/work/t401b; export TMPDIR=$W/tmp
B=$1; C=$2; shift 2
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
T=${TAG:-$B}; safe "$B"; safe "$C"; safe "$T"; : "${HOME:?}"
if [ -d /opt/nvidia ]; then . ~/alf-gpu/env/nvhpc.sh; else . ~/alf-gpu/env/intel.sh; export ONEAPI_DEVICE_SELECTOR=level_zero:0; fi
D="$W/runs/$T/$C"; rm -rf -- "$D"; mkdir -p "$W/runs/$T"; cp -r "$W/runs/base/$C" "$D"; cd "$D" || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_TARGET_OFFLOAD=MANDATORY
case $B in cpu*) unset OMP_TARGET_OFFLOAD;; esac
for kv in "$@"; do export "$kv"; done
s=$(date +%s.%N)
if [ "$GPULOCK" = 1 ]; then flock ~/alf-gpu/gpu.lock $W/bin/ALF_$B.out > run.log 2>&1; rc=$?; else $W/bin/ALF_$B.out > run.log 2>&1; rc=$?; fi
e=$(date +%s.%N); printf '%s %s rc=%s wall=%.1fs\n' $B $C $rc $(echo "$e - $s" | bc)
