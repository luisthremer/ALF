#!/bin/bash
# Ticket 205, run on nividia: build one variant in a fresh tree.
#   usage: t205_build.sh SRC NAME TOOL [configure args...]    e.g.  t205_build.sh src_new nv_def nvhpc PGI noMPI HDF5 GPU=nvidia
#   SRC  = src_new | src_head (under ~/alf-gpu/work/t205),  TOOL = gnu | nvhpc | intel (env file)
# Result: ~/alf-gpu/work/t205/bin/ALF_<NAME>.out, tree in tree_<NAME> (kept for ctest/inspection), logs in logs/.
# HDF5 builds are taken from a copy of ~/alf-gpu/src/ALF-*/HDF5 in work/t205/h5 (no rebuild).
W=$HOME/alf-gpu/work/t205; export TMPDIR=$W/tmp; mkdir -p $TMPDIR $W/logs $W/bin
SRC=$1; NAME=$2; TOOL=$3; shift 3
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
case $TOOL in gnu) . ~/alf-gpu/env/gnu.sh;; nvhpc) . ~/alf-gpu/env/nvhpc.sh;; intel) . ~/alf-gpu/env/intel.sh;; *) echo bad tool; exit 1;; esac
export ALF_HDF5_DIR=$W/h5
safe "$NAME"; safe "$SRC"; : "${HOME:?}"
rm -rf -- "$W/tree_$NAME"; mkdir -p "$W/tree_$NAME"; rsync -a "$W/$SRC/" "$W/tree_$NAME/"
cd "$W/tree_$NAME" || exit 1
. ./configure.sh "$@" NO-INTERACTIVE > $W/logs/conf_$NAME.log 2>&1 || { echo "configure failed"; tail -20 $W/logs/conf_$NAME.log; exit 2; }
env | grep '^ALF_' | sort > $W/logs/env_$NAME.txt
/usr/bin/time -f "make wall %e s" make -j${JOBS:-4} program > $W/logs/make_$NAME.log 2>&1 || { echo "make failed"; grep -i -B2 -A5 "error" $W/logs/make_$NAME.log | head -60; exit 3; }
tail -1 $W/logs/make_$NAME.log
cp Prog/ALF.out $W/bin/ALF_$NAME.out && echo "built $NAME"
