#!/bin/bash
# usage: bld.sh HOST VARIANT   HOST=nividia|intel  VARIANT=src|srcT   (srcT = -DTIMING)
HOST=$1; V=$2
cd /Users/luis/Documents/_WORK/Coding/CLAUDE/ALF/.claude/worktrees/gpu-performance-bottlenecks-bc7e27 || exit 1
ssh $HOST "mkdir -p ~/alf-gpu/work/dlsub/$V"
rsync -a --exclude='.git' --exclude='*.o' --exclude='*.mod' --exclude='*.a' --exclude='*.out' --exclude='HDF5' --exclude='docs' --exclude='Documentation' --exclude='Mathematica' --exclude='testsuite' --exclude='.claude' --exclude='gpu_bench' --exclude='PARKING.md' ./ $HOST:alf-gpu/work/dlsub/$V/
if [ $HOST = nividia ]; then ENV='. ~/alf-gpu/env/nvhpc.sh; MACHINE="PGI noMPI GPU=nvidia"'; else ENV='. ~/alf-gpu/env/intel.sh 2>/dev/null; MACHINE="IntelLLVM noMPI GPU=intel"'; fi
FL=""; [ $V = srcT ] && FL='export ALF_FLAGS_PROG="$ALF_FLAGS_PROG -DTIMING"'
ssh $HOST bash -s <<R
cd ~/alf-gpu/work/dlsub/$V; export TMPDIR=~/alf-gpu/work/dlsub/tmp; $ENV
. ./configure.sh \$MACHINE NO-INTERACTIVE >/dev/null 2>&1
$FL
make -j6 all > ../make_$V.log 2>&1; echo "$HOST $V rc=\$?"; grep -i "NVFORTRAN-S\|error #" ../make_$V.log | head -3
R
