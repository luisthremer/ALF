#!/bin/bash
# Ticket 206, nividia: GNU devel noMPI (GPU=none) build of src_new + the whole ctest suite (as in CI), -j4.
W=$HOME/alf-gpu/work/t206; export TMPDIR=$W/tmp; mkdir -p $TMPDIR $W/logs
. ~/alf-gpu/env/gnu.sh
: "${HOME:?}"
rm -rf -- "$W/tree_ctest"; mkdir -p "$W/tree_ctest"; rsync -a "$W/src_new/" "$W/tree_ctest/"; cd "$W/tree_ctest" || exit 1
. ./configure.sh GNU devel noMPI NO-INTERACTIVE > $W/logs/conf_ctest.log 2>&1 || { echo configure failed; exit 2; }
make -j4 all > $W/logs/make_ctest.log 2>&1 || { echo make failed; tail -20 $W/logs/make_ctest.log; exit 3; }
cd testsuite && mkdir -p tests && cd tests
cmake .. > $W/logs/cmake_ctest.log 2>&1 && cmake --build . -j4 >> $W/logs/cmake_ctest.log 2>&1 || { echo cmake build failed; exit 4; }
ctest -j4 --output-on-failure > $W/logs/ctest.log 2>&1; tail -6 $W/logs/ctest.log
