#!/bin/bash
# build + run t302a_qrbench.cpp (on intel; in dev/qrbench)
W=$HOME/alf-gpu/work/t302a/dev; export TMPDIR=$HOME/alf-gpu/work/t302a/tmp; mkdir -p $W/qrbench; cd $W/qrbench
. $HOME/alf-gpu/env/intel.sh 2>/dev/null
icpx -fsycl -O2 -std=c++17 -I$MKLROOT/include $W/src/docs/gpu/scripts/intel/t302a_qrbench.cpp -qmkl -lmkl_sycl_lapack -lmkl_sycl_blas -lsycl -o qrbench || exit 1
export ONEAPI_DEVICE_SELECTOR=level_zero:0
./qrbench
