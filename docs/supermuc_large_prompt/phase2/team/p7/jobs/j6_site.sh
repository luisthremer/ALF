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
#SBATCH --time=00:15:00
#SBATCH --ntasks-per-node=8
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
source $SCRATCH/large_prompt/team/p7/jobs/p7lib.sh
export W=$P7/work_j6; mkdir -p $W $P7/lib $P7/data
cd $P7/lib
icpx -fsycl -fiopenmp -fPIC -shared -O2 -std=c++17 -I$MKLROOT/include $P7/src/p7off.cpp -o libp7off.so -qmkl -lmkl_sycl_blas -ldl 2>&1 | tail -20
nm -n $BIN7 | grep -i " t " > $P7/data/ALF.out.nm.txt; wc -l $P7/data/ALF.out.nm.txt
cd $PJ7
PRE=$P7/lib/libalfsb.so:$P7/lib/libp7off.so
echo "=== single tile V7 serialised with call sites, 1 and 3 sweeps"
r1 V7_cs1 V7 1 P7PRE=$PRE ALFP7_PROF=1 ALFDEVQR_PROF=1 ALFSH_TIMING=1
r1 V7_cs3 V7 3 P7PRE=$PRE ALFP7_PROF=1 ALFDEVQR_PROF=1 ALFSH_TIMING=1; tm V7_cs1; tm V7_cs3
mkdir -p $P7/data/j6; for d in $W/runs/*; do t=$(basename $d); mkdir -p $P7/data/j6/$t; cp $d/rank*.log $d/info $P7/data/j6/$t/ 2>/dev/null; done
echo "end $(date -Is)"
