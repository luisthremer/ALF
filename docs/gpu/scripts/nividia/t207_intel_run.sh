#!/bin/bash
# Ticket 207, on `intel`: run the GPU=intel binary on one prepared config (self-test ALF_GPU_SELFTEST=4; the stabilization and tau paths are inactive on GPU=intel until ticket 302 (alf_la backend), the run goes on the CPU path).
#   usage: t207_intel_run.sh CONFIG [env assignments...]
# Checks first that nobody else is busy; B580 only (level_zero:0); 1 CPU thread.
W=$HOME/alf-gpu/work/t207; export TMPDIR=$W/tmp
. ~/alf-gpu/env/intel.sh
# 2026-10-28: the owner's oneAPI maintenance on `intel` removed MKL 2025.3 (ifx 2025.3 still present); fall back to MKL 2026.2 (reported)
if [ ! -d /opt/intel/oneapi/mkl/2025.3/lib ] && [ -f /opt/intel/oneapi/mkl/2026.2/env/vars.sh ]; then . /opt/intel/oneapi/mkl/2026.2/env/vars.sh > /dev/null 2>&1; echo "MKL fallback: $MKLROOT"; fi
C=$1; shift
case "$C" in ""|*..*|*/*) echo "unsafe name"; exit 1;; esac
: "${HOME:?}"
uptime; ps -eo user,pcpu,comm --sort=-pcpu | head -4
D="$W/runs/intel_$C"; rm -rf -- "$D"; mkdir -p "$D"; cp -r "$W/runs/base/$C/." "$D/"; cd "$D" || exit 1
export ONEAPI_DEVICE_SELECTOR=level_zero:0 OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1
export ALF_GPU_SELFTEST=4 ALF_GPU_EXIT=1
for kv in "$@"; do export "$kv"; done
timeout 600 $W/bin/ALF_intel.out > run.log 2>&1; echo "rc=$?"
grep -E "FAILED|PASSED|GPU tau|ONEMKL|rror" run.log info | cut -c1-300 | head -60
