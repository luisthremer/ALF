#!/bin/bash
# build + run the singular-getrf probe against build_TAG (on intel): t302a_singular.sh TAG
W=$HOME/alf-gpu/work/t302a/dev; B=$W/build_$1; export TMPDIR=$W/tmp
. $HOME/alf-gpu/env/intel.sh 2>/dev/null
cd $B/testsuite/tests/alf_la.tests
[ -x t302a_singular ] || ifx -O2 -fiopenmp -fopenmp-targets=spir64 -I$B/Libraries/Modules -I$MKLROOT/include $W/src/docs/gpu/scripts/intel/t302a_singular.F90 \
  $B/Libraries/Modules/modules_90.a -qmkl -lmkl_sycl_blas -lmkl_sycl_lapack -lsycl -lstdc++ -o t302a_singular || exit 1
export ONEAPI_DEVICE_SELECTOR=level_zero:0 OMP_TARGET_OFFLOAD=MANDATORY
./t302a_singular; echo "exit code $?"
