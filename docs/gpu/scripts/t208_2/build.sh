#!/bin/bash
# 207 close (decision 41): build one variant in a fresh tree.  usage: build.sh NAME TOOL(nvhpc|intel) configure-args...
# W=~/alf-gpu/work/t208-2, source in $W/src (git archive of gpu-port-207-final); binaries $W/bin/ALF_<NAME>.out; J = make jobs.
W=$HOME/alf-gpu/work/t208-2; export TMPDIR=$W/tmp; mkdir -p $TMPDIR $W/logs $W/bin
NAME=$1; TOOL=$2; shift 2
case "$NAME" in ""|*..*|*/*) echo "refusing unsafe name"; exit 1;; esac
[ -e "$W/bin/ALF_$NAME.out" ] && { echo "exists: $NAME"; exit 0; }
. ~/alf-gpu/env/$TOOL.sh 2>/dev/null
export ALF_HDF5_DIR=$HOME/alf-gpu/hdf5
rm -rf -- "$W/tree_$NAME"; mkdir -p "$W/tree_$NAME"; rsync -a "$W/src/" "$W/tree_$NAME/"
cd "$W/tree_$NAME" || exit 1
. ./configure.sh "$@" NO-INTERACTIVE > $W/logs/conf_$NAME.log 2>&1 || { echo "$NAME configure failed"; tail -20 $W/logs/conf_$NAME.log; exit 2; }
/usr/bin/time -f "make wall %e s" make -j${J:-4} program > $W/logs/make_$NAME.log 2>&1 || { echo "$NAME make failed"; grep -i -B2 -A8 "error" $W/logs/make_$NAME.log | head -60; exit 3; }
cp Prog/ALF.out $W/bin/ALF_$NAME.out && sha256sum $W/bin/ALF_$NAME.out > $W/bin/ALF_$NAME.out.sha256 && echo "built $NAME $(tail -1 $W/logs/make_$NAME.log)"
rm -rf -- "$W/tree_$NAME"
