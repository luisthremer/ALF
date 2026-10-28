#!/bin/bash
# Ticket 207, run on `intel` (not nividia): GPU=intel build of work/t207/src_new (as shipped by 201/203, no flag surgery).
W=$HOME/alf-gpu/work/t207; export TMPDIR=$W/tmp; mkdir -p $TMPDIR $W/logs $W/bin
. ~/alf-gpu/env/intel.sh
# 2026-10-28: the owner's oneAPI maintenance on `intel` removed MKL 2025.3 (ifx 2025.3 still present); fall back to MKL 2026.2 (reported)
if [ ! -d /opt/intel/oneapi/mkl/2025.3/lib ] && [ -f /opt/intel/oneapi/mkl/2026.2/env/vars.sh ]; then . /opt/intel/oneapi/mkl/2026.2/env/vars.sh > /dev/null 2>&1; echo "MKL fallback: $MKLROOT"; fi
: "${HOME:?}"
rm -rf -- "$W/tree_intel"; mkdir -p "$W/tree_intel"; rsync -a "$W/src_new/" "$W/tree_intel/"; cd "$W/tree_intel" || exit 1
. ./configure.sh IntelLLVM noMPI GPU=intel NO-INTERACTIVE > $W/logs/conf_intel.log 2>&1 || { echo configure failed; tail $W/logs/conf_intel.log; exit 2; }
env | grep '^ALF_' | sort > $W/logs/env_intel.txt
/usr/bin/time -f "make wall %e s" make -j${JOBS:-6} program > $W/logs/make_intel.log 2>&1 || { echo make failed; grep -i -B2 -A6 error $W/logs/make_intel.log | head -60; exit 3; }
tail -1 $W/logs/make_intel.log; cp Prog/ALF.out $W/bin/ALF_intel.out; echo built
