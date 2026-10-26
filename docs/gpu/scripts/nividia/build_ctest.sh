#!/bin/bash
# usage: build_ctest.sh gnu|intel   -> devel noMPI build + unit tests, as in .github/workflows/tests.yml
set -o pipefail
t=$1; M=$([ $t = gnu ] && echo GNU || echo IntelLLVM)
. ~/alf-gpu/env/$t.sh
cd ~/alf-gpu/src/ALF-$t
{ . ./configure.sh $M devel noMPI NO-INTERACTIVE; echo "ALF_FC=$ALF_FC"; echo "ALF_FLAGS_PROG=$ALF_FLAGS_PROG"; echo "ALF_LIB=$ALF_LIB"
  /usr/bin/time -f "make wall %e s" make -j8 all || exit 2
  rm -rf testsuite/tests && mkdir testsuite/tests && cd testsuite/tests && cmake .. -DCMAKE_BUILD_TYPE=Release >/dev/null && cmake --build . -j8 >/dev/null || exit 3
  ctest -j8 --output-on-failure 2>&1 | tail -25; } > ~/alf-gpu/logs/ctest_$t.log 2>&1
echo "exit $?" >> ~/alf-gpu/logs/ctest_$t.log
