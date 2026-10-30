#!/bin/bash
# build + run t302a_small against build_TAG: t302a_small.sh TAG  (on intel)
W=$HOME/alf-gpu/work/t302a/dev; B=$W/build_$1; export TMPDIR=$W/tmp
. $HOME/alf-gpu/env/intel.sh 2>/dev/null
cd $B/testsuite/tests/alf_la.tests
ifx -O2 -fiopenmp -fopenmp-targets=spir64 -I$B/Libraries/Modules -I$MKLROOT/include $W/src/docs/gpu/scripts/intel/t302a_small.F90 \
  $B/Libraries/Modules/modules_90.a -qmkl -lmkl_sycl_blas -lmkl_sycl_lapack -lsycl -lstdc++ -o t302a_small || exit 1
export ONEAPI_DEVICE_SELECTOR=level_zero:0 OMP_TARGET_OFFLOAD=MANDATORY
./t302a_small
