#!/bin/bash
#SBATCH --job-name=p3_bb1
#SBATCH --output=p3_bb1.%j.log
#SBATCH --error=p3_bb1.%j.err
#SBATCH --time=00:15:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# Phase 3b step 2: block sequence benchmark
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P=$SCRATCH/large_prompt/team/p3; cd $P
# the SYCL device image of the user kernels is only registered if the C++ part is linked by the sycl driver: shared library
icpx -fsycl -fPIC -shared -O2 -std=c++17 -I$MKLROOT/include blkbench.cpp -o libblkbench.so -L$MKLROOT/lib -lmkl_sycl_blas -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core > b_sycl_cpp.log 2>&1 || tail b_sycl_cpp.log
AOT="-fopenmp-targets=spir64_gen -Xopenmp-target-backend -device -Xopenmp-target-backend pvc"
ifx -O2 -cpp -fiopenmp $AOT -qmkl -I$MKLROOT/include blkbench.F90 -L. -lblkbench -Wl,-rpath,$P -lmkl_sycl_blas -lsycl -lstdc++ -o bbench > b_sycl_f.log 2>&1 || tail b_sycl_f.log
ls -l bbench || exit 1
export OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 OMP_NUM_THREADS=8
timeout 600 ./bbench 2000 2>&1 | tail -30
for ce in none main all; do echo "== LIBOMPTARGET_LEVEL_ZERO_USE_COPY_ENGINE=$ce"; LIBOMPTARGET_LEVEL_ZERO_USE_COPY_ENGINE=$ce timeout 300 ./bbench 1000 2>&1 | grep -E "^E[0-9]|OMP sync" ; done
echo "== STAGING_BUFFER_SIZE=1MB"; LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 timeout 300 ./bbench 1000 2>&1 | grep -E "^E[0-9]|OMP sync"
echo "done $(date -Is)"
