#!/bin/bash
# build + run the timing copy of t_order against build_TAG: t302a_tprof.sh TAG [args]   (on intel)
W=$HOME/alf-gpu/work/t302a/dev; B=$W/build_$1; export TMPDIR=$HOME/alf-gpu/work/t302a/tmp; shift
. $HOME/alf-gpu/env/intel.sh 2>/dev/null
cd $B/testsuite/tests/alf_la.tests
ifx -O2 -fiopenmp -fopenmp-targets=spir64 -I$B/Libraries/Modules -I. -I$MKLROOT/include $W/src/docs/gpu/scripts/intel/t302a_tprof.F90 \
  libla_test_util.a $B/Libraries/Modules/modules_90.a -qmkl -lmkl_sycl_blas -lmkl_sycl_lapack -lsycl -lstdc++ -o t302a_tprof || exit 1
export ONEAPI_DEVICE_SELECTOR=level_zero:0 OMP_TARGET_OFFLOAD=MANDATORY
/usr/bin/time -f "wall %e s user %U sys %S" ./t302a_tprof "$@"
