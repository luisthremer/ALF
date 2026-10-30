#!/bin/bash
# acceptance 5 (development version) on nividia: GPU=nvidia build of base (dev/base_src) and candidate (dev/src), alf_la.tests
# under CUDA, ctest -V output of both compared (timings dropped). usage: t302a_cuda.sh build | run   (run holds the GPU)
W=$HOME/alf-gpu/work/t302a/dev; export TMPDIR=$W/tmp; mkdir -p $TMPDIR
case $1 in
build)
  for v in base cand; do
    src=$W/base_src; [ $v = cand ] && src=$W/src
    rm -rf $W/cuda_$v; mkdir $W/cuda_$v; cp -a $src/. $W/cuda_$v/; cd $W/cuda_$v
    ( . $HOME/alf-gpu/env/nvhpc.sh >/dev/null 2>&1; . ./configure.sh PGI devel noMPI GPU=nvidia NO-INTERACTIVE > configure.log 2>&1 \
      && make -j${CUDA_JOBS:-4} all > make.log 2>&1 && cd testsuite && mkdir -p tests && cd tests && cmake .. > cmake.log 2>&1 && cmake --build . -j${CUDA_JOBS:-4} > cmake_build.log 2>&1 ) \
      || { echo "build $v failed"; exit 1; }
  done; echo BUILDS_OK ;;
run)
  nvidia-smi --query-compute-apps=pid --format=csv,noheader | grep -q . && { echo "GPU busy"; exit 3; }
  . $HOME/alf-gpu/env/nvhpc.sh >/dev/null 2>&1
  export OMP_TARGET_OFFLOAD=MANDATORY CUDA_VISIBLE_DEVICES=0
  for v in base cand; do
    ( cd $W/cuda_$v/testsuite/tests/alf_la.tests && ctest -V > $W/cuda_$v.ctest.log 2>&1; echo "$v ctest exit $?"; tail -4 $W/cuda_$v.ctest.log | head -2 )
  done
  diff <(grep -v -E "Test time|Total Test time|sec$|^\s*Start|Test project|tests passed|$W" $W/cuda_base.ctest.log | sed "s#$W/cuda_base#X#g") \
       <(grep -v -E "Test time|Total Test time|sec$|^\s*Start|Test project|tests passed|$W" $W/cuda_cand.ctest.log | sed "s#$W/cuda_cand#X#g") > $W/cuda_diff.txt
  echo "diff lines: $(wc -l < $W/cuda_diff.txt)"; grep -c "alf_la_mod.F90" $W/cuda_base.ctest.log ;;
esac
