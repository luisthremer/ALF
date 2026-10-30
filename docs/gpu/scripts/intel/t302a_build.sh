#!/bin/bash
# build GPU=intel from the synced tree into a fresh build dir on intel: t302a_build.sh TAG   (run on intel)
# tree: ~/alf-gpu/work/t302a/dev/build_TAG (build trees only are removable); log: dev/build_TAG.log
set -e
TAG=$1; W=$HOME/alf-gpu/work/t302a/dev; export TMPDIR=$W/tmp; mkdir -p $TMPDIR
rm -rf "$W/build_$TAG"; mkdir "$W/build_$TAG"; cp -a $W/src/. "$W/build_$TAG/"
cd "$W/build_$TAG"
. $HOME/alf-gpu/env/intel.sh 2>/dev/null
. ./configure.sh IntelLLVM devel noMPI GPU=intel NO-INTERACTIVE > configure.log 2>&1 || { tail -20 configure.log; exit 1; }
make -j6 all > make.log 2>&1 || { grep -n -B2 -A8 "rror" make.log | head -80; exit 1; }
cd testsuite && mkdir -p tests && cd tests
cmake .. > cmake.log 2>&1 && cmake --build . -j6 > cmake_build.log 2>&1 || { grep -n -B2 -A8 "rror" cmake_build.log | head -60; exit 1; }
echo BUILD_OK
