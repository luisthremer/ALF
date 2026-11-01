#!/bin/bash
# Ticket 207 (rebase on 212, development only), nividia: build one variant in a fresh tree.
#   usage: [ALF_FLAGS_EXT=-DQR_PREPIVOT] t207_d2_build.sh SRC NAME TOOL [configure args...]
#          e.g. t207_d2_build.sh src nv_def nvhpc PGI noMPI HDF5 GPU=nvidia ;  SRC = src | base_src (under dev2)
# Result: dev2/bin/ALF_<NAME>.out (+ .sha256); tree in dev2/tree_<NAME> (build tree, removable).
W=$HOME/alf-gpu/work/t207/dev2; export TMPDIR=$W/tmp; mkdir -p $TMPDIR $W/logs $W/bin
SRC=$1; NAME=$2; TOOL=$3; shift 3
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
case $TOOL in gnu) . ~/alf-gpu/env/gnu.sh;; nvhpc) . ~/alf-gpu/env/nvhpc.sh;; *) echo bad tool; exit 1;; esac
safe "$NAME"; safe "$SRC"; : "${HOME:?}"
[ -d $W/h5 ] || { mkdir -p $W/h5 && cp -a ~/alf-gpu/hdf5/. $W/h5/; }
export ALF_HDF5_DIR=$W/h5
rm -rf -- "$W/tree_$NAME"; mkdir -p "$W/tree_$NAME"; rsync -a "$W/$SRC/" "$W/tree_$NAME/"
cd "$W/tree_$NAME" || exit 1
. ./configure.sh "$@" NO-INTERACTIVE > $W/logs/conf_$NAME.log 2>&1 || { echo "configure failed"; tail -20 $W/logs/conf_$NAME.log; exit 2; }
env | grep '^ALF_' | sort > $W/logs/env_$NAME.txt
/usr/bin/time -f "make wall %e s" make -j4 program > $W/logs/make_$NAME.log 2>&1 || { echo "make failed"; grep -i -B2 -A5 "error" $W/logs/make_$NAME.log | head -60; exit 3; }
tail -1 $W/logs/make_$NAME.log
cp Prog/ALF.out $W/bin/ALF_$NAME.out && sha256sum $W/bin/ALF_$NAME.out > $W/bin/ALF_$NAME.out.sha256 && echo "built $NAME $(cut -c1-16 $W/bin/ALF_$NAME.out.sha256)"
