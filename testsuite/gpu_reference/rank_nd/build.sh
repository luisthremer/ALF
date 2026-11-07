#!/bin/bash
# Ticket 401 (development), on nividia or intel: build one variant in a fresh tree from work/t401b/src.
#   usage: build.sh NAME TOOL [configure args...]     TOOL = gnu | nvhpc | intel
#   e.g.   build.sh nv nvhpc PGI noMPI GPU=nvidia   |   build.sh cpu gnu GNU noMPI   |   build.sh ifx intel IntelLLVM noMPI GPU=intel
# ALF_FLAGS_EXT from the caller is kept (planted-defect builds: -DALF_PLANT_...).  Result: bin/ALF_<NAME>.out, tree_<NAME>.
W=$HOME/alf-gpu/work/${T401:-t401b}; export TMPDIR=$W/tmp; mkdir -p $TMPDIR $W/logs $W/bin
NAME=$1; TOOL=$2; shift 2
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
case $TOOL in gnu) . ~/alf-gpu/env/gnu.sh;; nvhpc) . ~/alf-gpu/env/nvhpc.sh;; intel) . ~/alf-gpu/env/intel.sh;; *) echo bad tool; exit 1;; esac
safe "$NAME"; : "${HOME:?}"
[ -d $W/h5 ] || { mkdir -p $W/h5 && cp -a ~/alf-gpu/hdf5/. $W/h5/; }
export ALF_HDF5_DIR=$W/h5
rm -rf -- "$W/tree_$NAME"; mkdir -p "$W/tree_$NAME"; rsync -a "$W/src/" "$W/tree_$NAME/"
cd "$W/tree_$NAME" || exit 1
. ./configure.sh "$@" NO-INTERACTIVE > $W/logs/conf_$NAME.log 2>&1 || { echo "configure failed"; tail -20 $W/logs/conf_$NAME.log; exit 2; }
env | grep '^ALF_' | sort > $W/logs/env_$NAME.txt
/usr/bin/time -f "make wall %e s" make -j${JOBS:-4} program > $W/logs/make_$NAME.log 2>&1 || { echo "make failed"; grep -i -B2 -A5 "error" $W/logs/make_$NAME.log | head -60; exit 3; }
tail -1 $W/logs/make_$NAME.log
cp Prog/ALF.out $W/bin/ALF_$NAME.out && echo "built $NAME"
