#!/bin/bash
# build + run the 302a probe on intel: t302a_probe.sh [own]   (run from ~/alf-gpu/work/t302a/dev/probe)
set -e
export TMPDIR=$HOME/alf-gpu/work/t302a/dev/tmp; mkdir -p $TMPDIR
. $HOME/alf-gpu/env/intel.sh 2>/dev/null
S=$HOME/alf-gpu/work/t302a/dev/src
icpx -fsycl -O2 -std=c++17 -I$MKLROOT/include -c $S/Libraries/Modules/alf_la_sycl_shim.cpp -o shim.o
ifx -O2 -g -fiopenmp -fopenmp-targets=spir64 -qmkl -I$MKLROOT/include -c $S/docs/gpu/scripts/intel/t302a_probe.F90 -o probe.o
ifx -O2 -g -fiopenmp -fopenmp-targets=spir64 -qmkl probe.o shim.o -lmkl_sycl_blas -lmkl_sycl_lapack -lsycl -lstdc++ -o probe
export ONEAPI_DEVICE_SELECTOR=level_zero:0 OMP_TARGET_OFFLOAD=MANDATORY
echo built
