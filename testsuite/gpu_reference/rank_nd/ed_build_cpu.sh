#!/bin/bash
# t401ed: build the CPU GNU HDF5 noMPI reference binary (GPU=none, upstream behaviour) from the synced tree.
# Runs ON nividia.  Writes only inside ~/alf-gpu/work/t401ed/.   usage: ed_build_cpu.sh
W=$HOME/alf-gpu/work/t401ed; export TMPDIR=$W/tmp; mkdir -p $TMPDIR $W/logs $W/bin
. ~/alf-gpu/env/gnu.sh
[ -d $W/h5 ] || { mkdir -p $W/h5 && cp -a ~/alf-gpu/hdf5/. $W/h5/; }
export ALF_HDF5_DIR=$W/h5
rm -rf -- "$W/tree_cpu"; mkdir -p "$W/tree_cpu"; rsync -a "$W/src/" "$W/tree_cpu/"
cd "$W/tree_cpu" || exit 1
. ./configure.sh GNU noMPI HDF5 NO-INTERACTIVE > $W/logs/conf_cpu.log 2>&1 || { echo "configure failed"; tail -20 $W/logs/conf_cpu.log; exit 2; }
env | grep '^ALF_' | sort > $W/logs/env_cpu.txt
/usr/bin/time -f "make wall %e s" make -j4 all > $W/logs/make_cpu.log 2>&1 || { grep -i -B2 -A5 "error" $W/logs/make_cpu.log | head -60; exit 3; }
tail -1 $W/logs/make_cpu.log
cp Prog/ALF.out $W/bin/ALF_cpu.out && cp Analysis/ana.out $W/bin/ana_cpu.out 2>/dev/null; echo "built cpu"
