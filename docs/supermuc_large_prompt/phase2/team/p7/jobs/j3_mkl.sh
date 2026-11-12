#!/bin/bash
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --exclude=i20r02c06s02
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-node=8
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
source $SCRATCH/large_prompt/team/p7/jobs/p7lib.sh
export W=$P7/work_j3; mkdir -p $W $P7/lib $P7/data
cd $P7/lib && icpx -fsycl -fPIC -shared -O2 -std=c++17 -I$MKLROOT/include $P7/src/p7mkl.cpp -o libp7mkl.so -qmkl -lmkl_sycl_blas -ldl 2>&1 | tail -20
ls -la $P7/lib; nm -D --defined-only $P7/lib/libp7mkl.so | grep -c "column_major"
cd $PJ7
PRE=$P7/lib/libalfsb.so:$P7/lib/libp7mkl.so
for v in V7; do
  echo "=== $v serialised full profile with the MKL interposer, 1 and 3 sweeps"
  r1 ${v}_mf1 $v 1 LIBOMPTARGET_PLUGIN_PROFILE=T P7PRE=$PRE ALFP7_PROF=1 ALFDEVQR_PROF=1 ALFSH_TIMING=1
  r1 ${v}_mf3 $v 3 LIBOMPTARGET_PLUGIN_PROFILE=T P7PRE=$PRE ALFP7_PROF=1 ALFDEVQR_PROF=1 ALFSH_TIMING=1; tm ${v}_mf1; tm ${v}_mf3
  grep -h "ALFP7" $W/runs/${v}_mf3/rank0.log | head -40
done
echo "=== V7 interposer only (passes through, no ALFP7_PROF): overhead check"; r1 V7_ip3 V7 3 P7PRE=$PRE; tm V7_ip3
mkdir -p $P7/data/j3; for d in $W/runs/*; do t=$(basename $d); mkdir -p $P7/data/j3/$t; cp $d/rank0.log $d/info $P7/data/j3/$t/ 2>/dev/null; done
echo "end $(date -Is)"
