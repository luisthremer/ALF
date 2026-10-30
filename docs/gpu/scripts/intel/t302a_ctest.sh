#!/bin/bash
# run ctest of a build on the B580: t302a_ctest.sh TAG [ctest args]   (run on intel; ONEAPI selector + MANDATORY set here)
TAG=$1; shift; ARGS=("$@")   # (sourcing intel.sh resets "$@")
W=$HOME/alf-gpu/work/t302a/dev
. $HOME/alf-gpu/env/intel.sh 2>/dev/null
export ONEAPI_DEVICE_SELECTOR=level_zero:0 OMP_TARGET_OFFLOAD=MANDATORY FORT_BUFFERED=false
cd $W/build_$TAG/testsuite/tests/alf_la.tests || exit 2
ctest --output-on-failure "${ARGS[@]}"
