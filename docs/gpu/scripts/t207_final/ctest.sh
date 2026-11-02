#!/bin/bash
# 207 close, check 5 (nividia): GNU devel noMPI build + ctest -j4 of $W/src (candidate) and $W/src_up (upstream master) -> logs/ctest_<name>.txt
W=$HOME/alf-gpu/work/t207-final; export TMPDIR=$W/tmp; . ~/alf-gpu/env/gnu.sh
export ALF_HDF5_DIR=$HOME/alf-gpu/hdf5
for n in cand:src up:src_up; do
  name=${n%%:*}; s=${n##*:}; T=$W/tree_ct_$name; rm -rf $T; mkdir -p $T; rsync -a $W/$s/ $T/; cd $T || exit 1
  { . ./configure.sh GNU devel noMPI NO-INTERACTIVE; make -j4 all || exit 2
    rm -rf testsuite/tests; mkdir testsuite/tests; cd testsuite/tests && cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null && cmake --build . -j4 > /dev/null || exit 3
    ctest -j4 2>&1 | tail -15; } > $W/logs/ctest_$name.txt 2>&1
  cd $W; rm -rf $T
done
echo CTEST_DONE >> $W/logs/ctest_up.txt
