#!/bin/bash
# usage: ctest108.sh NAME "EXTRA_FLAGS"  -> devel noMPI GNU build + unit tests, log in work/t108/ctest_NAME.log
export TMPDIR=$HOME/alf-gpu/work/t108/tmp
. ~/alf-gpu/env/gnu.sh
cd ~/alf-gpu/work/t108/src
make clean >/dev/null 2>&1
if [ -n "$2" ]; then export ALF_FLAGS_EXT="$2"; else unset ALF_FLAGS_EXT; fi
{ . ./configure.sh GNU devel noMPI NO-INTERACTIVE; echo "FLAGS=$ALF_FLAGS_PROG"
  make -j6 all || exit 2
  rm -rf testsuite/tests && mkdir testsuite/tests && cd testsuite/tests && cmake .. -DCMAKE_BUILD_TYPE=Release >/dev/null && cmake --build . -j6 >/dev/null || exit 3
  ctest -j6 --output-on-failure 2>&1 | tail -8; } > ~/alf-gpu/work/t108/ctest_$1.log 2>&1
