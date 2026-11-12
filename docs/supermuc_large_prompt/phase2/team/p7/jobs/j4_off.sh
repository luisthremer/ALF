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
export W=$P7/work_j4; mkdir -p $W $P7/lib $P7/data
cd $P7/lib
icpx -fsycl -fPIC -shared -O2 -std=c++17 -fsycl-targets=spir64_gen -Xs "-device pvc" -I$MKLROOT/include $P7/src/alf_la_sycl_blk_p7.cpp -o libalfsb.so -qmkl -lmkl_sycl_blas 2>&1 | tail -5
icpx -fsycl -fiopenmp -fPIC -shared -O2 -std=c++17 -I$MKLROOT/include $P7/src/p7off.cpp -o libp7off.so -qmkl -lmkl_sycl_blas -ldl 2>&1 | tail -20
ls -la $P7/lib
cd $PJ7
echo "=== V7 host submission cost per op (async mode, ALFP7_SUBMIT=1), 1 and 3 sweeps"
r1 V7_s1 V7 1 P7PRE=$P7/lib/libalfsb.so ALFP7_SUBMIT=1; r1 V7_s3 V7 3 P7PRE=$P7/lib/libalfsb.so ALFP7_SUBMIT=1; tm V7_s1; tm V7_s3; grep -h ALFP7_SUBMIT $W/runs/V7_s3/rank0.log
echo "=== V7 serialised full profile incl. dispatched BLAS (libp7off), 1 and 3 sweeps"
PRE=$P7/lib/libalfsb.so:$P7/lib/libp7off.so
r1 V7_of1 V7 1 LIBOMPTARGET_PLUGIN_PROFILE=T P7PRE=$PRE ALFP7_PROF=1 ALFDEVQR_PROF=1 ALFSH_TIMING=1
r1 V7_of3 V7 3 LIBOMPTARGET_PLUGIN_PROFILE=T P7PRE=$PRE ALFP7_PROF=1 ALFDEVQR_PROF=1 ALFSH_TIMING=1; tm V7_of1; tm V7_of3; pg V7_of3
grep -h "ALFP7_OFF" $W/runs/V7_of3/rank0.log | sort -k4 -n -r | head -40
echo "=== V7 with libp7off pass-through (no ALFP7_PROF): overhead and bitwise check"
r1 V7_pt3 V7 3 P7PRE=$PRE; tm V7_pt3; pg V7_pt3
mkdir -p $P7/data/j4; for d in $W/runs/*; do t=$(basename $d); mkdir -p $P7/data/j4/$t; cp $d/rank0.log $d/info $P7/data/j4/$t/ 2>/dev/null; done
echo "end $(date -Is)"
