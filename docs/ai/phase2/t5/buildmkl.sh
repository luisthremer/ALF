#!/bin/bash
n=$1; shift
M=/opt/intel/oneapi/mkl/latest/lib
cd ~/work/ALF_P2/t5 && rm -rf src_$n && mkdir src_$n && (cd ~/work/ALF_P2/wt-cpu-s1 && git archive HEAD) | tar -x -C src_$n
cd src_$n && . ./configure.sh GNU "$@" GPU=host NO-INTERACTIVE > ../cfg_$n.log 2>&1
export ALF_LIB="${ALF_LIB/-llapack -lblas/-L$M -Wl,-rpath,$M -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl}"
echo "ALF_LIB=$ALF_LIB" >> ../cfg_$n.log
make -j16 lib ana program > ../build_$n.log 2>&1 && echo "$n OK" || echo "$n FAIL"
