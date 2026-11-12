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
export W=$P7/work_j3b; mkdir -p $W $P7/lib $P7/data
cd $P7/lib && icpx -fsycl -fPIC -shared -O2 -std=c++17 -I$MKLROOT/include $P7/src/p7mkl.cpp -o libp7mkl.so -qmkl -lmkl_sycl_blas -ldl 2>&1 | tail -20
ls -la $P7/lib; nm -D --defined-only $P7/lib/libp7mkl.so | grep -c "column_major"
cd $PJ7
PRE=$P7/lib/libalfsb.so:$P7/lib/libp7mkl.so
echo "=== V7 MKL interposer: (a) with PLUGIN_PROFILE=T (interop queue may profile), (b) without; 1 and 3 sweeps"
for m in a b; do
  PP=""; [ $m = a ] && PP="LIBOMPTARGET_PLUGIN_PROFILE=T"
  r1 V7_m${m}1 V7 1 $PP P7PRE=$PRE ALFP7_PROF=1 ALFDEVQR_PROF=1 ALFSH_TIMING=1
  r1 V7_m${m}3 V7 3 $PP P7PRE=$PRE ALFP7_PROF=1 ALFDEVQR_PROF=1 ALFSH_TIMING=1; tm V7_m${m}1; tm V7_m${m}3
  grep -h "ALFP7_MKL" $W/runs/V7_m${m}3/rank0.log | head -30
done
mkdir -p $P7/data/j3b; for d in $W/runs/*; do t=$(basename $d); mkdir -p $P7/data/j3b/$t; cp $d/rank0.log $d/info $P7/data/j3b/$t/ 2>/dev/null; done
echo "end $(date -Is)"
