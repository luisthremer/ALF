#!/bin/bash -l
# build a candidate tree (default and LOG) on the Helma login node: lib, program, testsuite
T=$1
for v in "" _log; do
  cd ~/p2_nvidia/cand_$T$v || exit 1
  stab=""; [ "$v" = _log ] && stab=LOG
  . ./configure.sh HELMAGPU noMPI $stab GPU=nvidia NO-INTERACTIVE > cfg.log 2>&1
  module load cmake/3.30.5 >/dev/null 2>&1
  make -j16 all > build.log 2>&1 || { echo "build FAIL $v"; tail build.log; continue; }
  cmake -S testsuite -B testsuite/tests > cmake.log 2>&1 && cmake --build testsuite/tests -j16 > tbuild.log 2>&1 && echo "build ok $v" || { echo "tests build FAIL $v"; grep -i -B2 -A5 error tbuild.log | head -30; }
done
