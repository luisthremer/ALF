#!/bin/bash
# Ticket 212 (development only), nividia: ctest of src.  usage: t212_ctest.sh gnu|cuda [ctest args]
#   gnu : GNU devel noMPI (GPU=none), whole suite as in CI.   cuda: NVHPC devel noMPI GPU=nvidia, the alf_la tests only.
W=$HOME/alf-gpu/work/t212/dev; export TMPDIR=$W/tmp; mkdir -p $TMPDIR $W/logs
K=$1; shift; : "${HOME:?}"
[ -d $W/h5 ] || { mkdir -p $W/h5 && cp -a ~/alf-gpu/hdf5/. $W/h5/; }
export ALF_HDF5_DIR=$W/h5
if [ $K = gnu ]; then . ~/alf-gpu/env/gnu.sh; A="GNU devel noMPI"; else . ~/alf-gpu/env/nvhpc.sh; A="PGI devel noMPI GPU=nvidia"; export OMP_TARGET_OFFLOAD=MANDATORY; fi
rm -rf -- "$W/tree_ctest_$K"; mkdir -p "$W/tree_ctest_$K"; rsync -a "$W/src/" "$W/tree_ctest_$K/"; cd "$W/tree_ctest_$K" || exit 1
. ./configure.sh $A NO-INTERACTIVE > $W/logs/conf_ctest_$K.log 2>&1 || { echo configure failed; exit 2; }
make -j4 all > $W/logs/make_ctest_$K.log 2>&1 || { echo make failed; tail -20 $W/logs/make_ctest_$K.log; exit 3; }
cd testsuite && mkdir -p tests && cd tests
cmake .. > $W/logs/cmake_ctest_$K.log 2>&1 && cmake --build . -j4 >> $W/logs/cmake_ctest_$K.log 2>&1 || { echo cmake build failed; tail $W/logs/cmake_ctest_$K.log; exit 4; }
ctest -j1 --output-on-failure "$@" > $W/logs/ctest_$K.log 2>&1; tail -8 $W/logs/ctest_$K.log
