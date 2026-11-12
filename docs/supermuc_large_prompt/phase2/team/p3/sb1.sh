#!/bin/bash
#SBATCH --job-name=p3_sb1
#SBATCH --output=p3_sb1.%j.log
#SBATCH --error=p3_sb1.%j.err
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
# Gate A': SYCL submission on the OpenMP interop queue vs OpenMP sync baseline
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P=$SCRATCH/large_prompt/team/p3; cd $P
# the SYCL device image of the user kernels is only registered if the C++ part is linked by the sycl driver: shared library
icpx -fsycl -fPIC -shared -O2 -std=c++17 -I$MKLROOT/include syclbench.cpp -o libsyclbench.so -L$MKLROOT/lib -lmkl_sycl_blas -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core > b_sycl_cpp.log 2>&1 || tail b_sycl_cpp.log
AOT="-fopenmp-targets=spir64_gen -Xopenmp-target-backend -device -Xopenmp-target-backend pvc"
ifx -O2 -cpp -fiopenmp $AOT -qmkl -I$MKLROOT/include syclbench.F90 -L. -lsyclbench -Wl,-rpath,$P -lmkl_sycl_blas -lsycl -lstdc++ -o sbench > b_sycl_f.log 2>&1 || tail b_sycl_f.log
ls -l sbench || exit 1
export OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 OMP_NUM_THREADS=8
for hw in 0 400; do echo "=== hw=$hw"; timeout 300 ./sbench 300 $hw 2>&1 | tail -20; done
echo "=== LIBOMPTARGET_INTEROP_INORDER_SYCL_QUEUE=1"
LIBOMPTARGET_INTEROP_INORDER_SYCL_QUEUE=1 timeout 300 ./sbench 300 0 2>&1 | tail -20
echo "done $(date -Is)"
