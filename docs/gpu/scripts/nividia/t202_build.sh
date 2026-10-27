#!/bin/bash
# Ticket 202, run on nividia: build the tree for one toolchain and run the whole ctest suite.
#   usage: t202_build.sh gnu|intel|pgi|cuda|cudahand [noctest]
#   gnu   : GNU devel noMPI            (HOST backend)
#   intel : IntelLLVM devel noMPI      (HOST backend)
#   pgi   : PGI devel noMPI            (HOST backend)
#   cuda  : PGI noMPI GPU=nvidia (ticket 201 switch)
#   cudahand: PGI noMPI + the flags hand-exported (identical result expected, compared via logs/env_*.txt):
#           -mp -> -mp=gpu -gpu=cc80, + -DALF_GPU -DALF_LA_CUDA, ALF_LIB + -mp=gpu -gpu=cc80 -cudalib=cublas,cusolver
# Source: ~/alf-gpu/work/t202/src (see t202_sync.sh); builds in tree_<tool>; logs in ~/alf-gpu/work/t202/logs.
# NB: configure.sh is sourced, never piped (the exports would be lost in a subshell).
T=$1; W=$HOME/alf-gpu/work/t202; export TMPDIR=$W/tmp; mkdir -p $TMPDIR $W/logs
L=$W/logs
case $T in
  gnu)   . ~/alf-gpu/env/gnu.sh;    M=GNU;       X="devel";;
  intel) . ~/alf-gpu/env/intel.sh;  M=IntelLLVM; X="devel";;
  pgi)   . ~/alf-gpu/env/nvhpc.sh;  M=PGI;       X="devel";;
  cuda)  . ~/alf-gpu/env/nvhpc.sh;  M=PGI;       X="GPU=nvidia";;
  cudahand) . ~/alf-gpu/env/nvhpc.sh; M=PGI;      X="";;
  *) echo "bad tool"; exit 1;;
esac
rm -rf $W/tree_$T; mkdir -p $W/tree_$T; rsync -a $W/src/ $W/tree_$T/
cd $W/tree_$T || exit 1
. ./configure.sh $M $X noMPI NO-INTERACTIVE > $L/conf_$T.log 2>&1 || { echo "configure failed"; exit 2; }
if [ $T = cudahand ]; then
  sub() { local o=""; for t in $1; do if [ "$t" = "-mp" ]; then o="$o -mp=gpu -gpu=cc80"; else o="$o $t"; fi; done; echo $o; }
  export ALF_FLAGS_MODULES="$(sub "$ALF_FLAGS_MODULES") -DALF_GPU -DALF_LA_CUDA"
  export ALF_FLAGS_PROG="$(sub "$ALF_FLAGS_PROG") -DALF_GPU -DALF_LA_CUDA"
  export ALF_LIB="$ALF_LIB -mp=gpu -gpu=cc80 -cudalib=cublas,cusolver"
fi
env | grep '^ALF_' | sort > $L/env_$T.txt
/usr/bin/time -f "make wall %e s" make -j4 all > $L/make_$T.log 2>&1 || { echo "make failed"; tail -30 $L/make_$T.log; exit 3; }
tail -1 $L/make_$T.log
[ "$2" = noctest ] && exit 0
cd testsuite && rm -rf tests && mkdir tests && cd tests
cmake .. > $L/cmake_$T.log 2>&1 && cmake --build . -j4 >> $L/cmake_$T.log 2>&1 || { echo "cmake build failed"; tail -30 $L/cmake_$T.log; exit 4; }
ctest -j4 --output-on-failure > $L/ctest_$T.log 2>&1
tail -6 $L/ctest_$T.log
