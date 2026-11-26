#!/bin/bash
# build.sh <name> <configure args...>
n=$1; shift
cd ~/work/ALF_P2/t5 && rm -rf src_$n && mkdir src_$n && (cd ~/work/ALF_P2/wt-cpu-s1 && git archive HEAD) | tar -x -C src_$n
cd src_$n && . ./configure.sh GNU "$@" GPU=host NO-INTERACTIVE > ../cfg_$n.log 2>&1 && make -j16 lib ana program > ../build_$n.log 2>&1 && echo "$n OK" || echo "$n FAIL"
