#!/bin/bash -l
# Build the Step-1 micro-benchmarks on the login node (bash -l ib_build.sh). Toolchain = the verified ALF build toolchain.
set -e
D=$SCRATCH/large_prompt/batched
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
cd $D/src
icpx --version | head -1
icpx -fsycl -O3 -qmkl=sequential ib_bench.cpp -o $D/bin/ib_bench -lmkl_sycl_blas 2>&1 | grep -v warning || true
icpx -fsycl -O2 -qmkl=sequential -I$D/src ib_wtest.cpp -o $D/bin/ib_wtest -lmkl_sycl_blas 2>&1 | grep -v warning || true
ifx -O2 -fiopenmp -fopenmp-targets=spir64 ib_omp.F90 -o $D/bin/ib_omp
ls -la $D/bin
