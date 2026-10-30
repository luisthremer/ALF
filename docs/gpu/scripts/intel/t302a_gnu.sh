#!/bin/bash
# GNU devel noMPI build + ctest of the synced candidate tree on nividia (acceptance 4: ctest count 52): t302a_gnu.sh
W=$HOME/alf-gpu/work/t302a/dev; export TMPDIR=$W/tmp; mkdir -p $TMPDIR
rm -rf $W/build_gnu; mkdir $W/build_gnu; cp -a $W/src/. $W/build_gnu/; cd $W/build_gnu
. $HOME/alf-gpu/env/gnu.sh >/dev/null 2>&1
. ./configure.sh GNU devel noMPI NO-INTERACTIVE > configure.log 2>&1 || { tail configure.log; exit 1; }
make -j4 all > make.log 2>&1 || { tail -30 make.log; exit 1; }
cd testsuite && mkdir -p tests && cd tests && cmake .. > cmake.log 2>&1 && cmake --build . -j4 > cmake_build.log 2>&1 || exit 1
ctest -j4 > ctest.log 2>&1; tail -5 ctest.log
