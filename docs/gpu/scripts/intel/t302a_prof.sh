#!/bin/bash
# time single tests of build_TAG: t302a_prof.sh TAG test... (on intel); extra env (e.g. LIBOMPTARGET_PROFILE=T) is inherited
TAG=$1; shift
W=$HOME/alf-gpu/work/t302a/dev; export TMPDIR=$HOME/alf-gpu/work/t302a/tmp
. $HOME/alf-gpu/env/intel.sh 2>/dev/null
export ONEAPI_DEVICE_SELECTOR=level_zero:0 OMP_TARGET_OFFLOAD=MANDATORY
cd $W/build_$TAG/testsuite/tests/alf_la.tests
for t in "$@"; do
  s=$(date +%s.%N); ./$t > $TMPDIR/prof_$t.out 2>&1; rc=$?; e=$(date +%s.%N)
  echo "$t rc=$rc $(echo "$e - $s" | bc) s"
done
