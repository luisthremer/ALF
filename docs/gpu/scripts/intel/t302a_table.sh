#!/bin/bash
# acceptance 3 on intel: build the table tool against build_TAG, run timing (TABLE lines) and the MKL_VERBOSE trace.
# usage: t302a_table.sh TAG OUTDIR      (OUTDIR is created; raw outputs stay)
TAG=$1; OUT=$2; W=$HOME/alf-gpu/work/t302a/dev; B=$W/build_$TAG; export TMPDIR=$W/tmp
. $HOME/alf-gpu/env/intel.sh 2>/dev/null
mkdir -p "$OUT" && cd "$OUT" || exit 2
ifx -O2 -fiopenmp -fopenmp-targets=spir64 -I$B/Libraries/Modules -I$MKLROOT/include $W/src/docs/gpu/scripts/intel/t302a_table.F90 \
  $B/Libraries/Modules/modules_90.a -qmkl -lmkl_sycl_blas -lmkl_sycl_lapack -lsycl -lstdc++ -o t302a_table || exit 1
export ONEAPI_DEVICE_SELECTOR=level_zero:0 OMP_TARGET_OFFLOAD=MANDATORY FORT_BUFFERED=false
./t302a_table > table_timing.txt 2>&1; echo "timing exit $?"
MKL_VERBOSE=1 ./t302a_table trace > table_trace.txt 2>&1; echo "trace exit $?"
