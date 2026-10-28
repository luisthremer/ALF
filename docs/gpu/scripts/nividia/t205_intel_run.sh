#!/bin/bash
# Ticket 205, on `intel`: run the GPU=intel binary on one prepared config: the one-slice update test (ALF_GPU_UPDTEST=1) + exit.
#   usage: t205_intel_run.sh CONFIG [env assignments...]
# Checks first that nobody else is busy; B580 only (level_zero:0); 1 CPU thread.
W=$HOME/alf-gpu/work/t205; export TMPDIR=$W/tmp
. ~/alf-gpu/env/intel.sh
C=$1; shift
case "$C" in ""|*..*|*/*) echo "unsafe name"; exit 1;; esac
: "${HOME:?}"
uptime; ps -eo user,pcpu,comm --sort=-pcpu | head -4
D="$W/runs/intel_$C"; rm -rf -- "$D"; mkdir -p "$D"; cp -r "$W/runs/base/$C/." "$D/"; cd "$D" || exit 1
export ONEAPI_DEVICE_SELECTOR=level_zero:0 OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1
export ALF_GPU_UPDTEST=1 ALF_GPU_EXIT=1
for kv in "$@"; do export "$kv"; done
timeout 600 $W/bin/ALF_intel.out > run.log 2>&1; echo "rc=$?"
grep -E "FAILED|PASSED|GPU update|ONEMKL|rror" run.log | cut -c1-300 | head -60
