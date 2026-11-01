#!/bin/bash
# Ticket 207-1 (development only), on intel: build one variant in a fresh tree.  usage: [ALF_FLAGS_EXT=..] t207_1_ibuild.sh NAME [configure args...]
# Result: bin/ALF_<NAME>.out (+ .sha256); tree tree_<NAME> (build tree, removable); logs/{conf,make}_<NAME>.log.  J = make jobs (default 6, <= 8 rule).
W=$HOME/alf-gpu/work/t207-1; export TMPDIR=$W/tmp; mkdir -p $TMPDIR $W/logs $W/bin
NAME=$1; shift
case "$NAME" in ""|*..*|*/*) echo "refusing unsafe name"; exit 1;; esac; : "${HOME:?}"
[ -e "$W/bin/ALF_$NAME.out" ] && { echo "refusing: binary $NAME exists (binaries are kept)"; exit 1; }
. ~/alf-gpu/env/intel.sh 2>/dev/null
[ -d $W/h5 ] || { mkdir -p $W/h5 && cp -a ~/alf-gpu/hdf5/. $W/h5/; }
export ALF_HDF5_DIR=$W/h5
rm -rf -- "$W/tree_$NAME"; mkdir -p "$W/tree_$NAME"; rsync -a "$W/${SRC:-src}/" "$W/tree_$NAME/"
cd "$W/tree_$NAME" || exit 1
. ./configure.sh "$@" NO-INTERACTIVE > $W/logs/conf_$NAME.log 2>&1 || { echo "configure failed"; tail -20 $W/logs/conf_$NAME.log; exit 2; }
env | grep '^ALF_' | sort > $W/logs/env_$NAME.txt
/usr/bin/time -f "make wall %e s" make -j${J:-6} program > $W/logs/make_$NAME.log 2>&1 || { echo "make failed"; grep -i -B2 -A8 "error" $W/logs/make_$NAME.log | head -80; exit 3; }
tail -1 $W/logs/make_$NAME.log
cp Prog/ALF.out $W/bin/ALF_$NAME.out && sha256sum $W/bin/ALF_$NAME.out > $W/bin/ALF_$NAME.out.sha256 && echo "built $NAME $(cut -c1-16 $W/bin/ALF_$NAME.out.sha256)"
