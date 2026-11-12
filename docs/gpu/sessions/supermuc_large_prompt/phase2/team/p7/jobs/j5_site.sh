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
#SBATCH --time=00:40:00
#SBATCH --ntasks-per-node=64
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
source $SCRATCH/large_prompt/team/p7/jobs/p7lib.sh
export W=$P7/work_j5; mkdir -p $W $P7/lib $P7/data
cd $P7/lib
icpx -fsycl -fiopenmp -fPIC -shared -O2 -std=c++17 -I$MKLROOT/include $P7/src/p7off.cpp -o libp7off.so -qmkl -lmkl_sycl_blas -ldl 2>&1 | tail -20
nm -n $BIN7 | grep -i " t " > $P7/data/ALF.out.nm.txt; wc -l $P7/data/ALF.out.nm.txt
cd $PJ7
PRE=$P7/lib/libalfsb.so:$P7/lib/libp7off.so
echo "=== single tile V7 serialised with call sites, 1 and 3 sweeps"
r1 V7_cs1 V7 1 P7PRE=$PRE ALFP7_PROF=1 ALFDEVQR_PROF=1 ALFSH_TIMING=1
r1 V7_cs3 V7 3 P7PRE=$PRE ALFP7_PROF=1 ALFDEVQR_PROF=1 ALFSH_TIMING=1; tm V7_cs1; tm V7_cs3
echo "=== node layout (32 ranks, 4/tile, CCS 4) V7: normal, then serialised profile (ops measured under real tile sharing), 1 and 3 sweeps"
r32 N_V7_n3 V7 3 ALFSH_TIMING=1; tm N_V7_n3
r32 N_V7_cs1 V7 1 P7PRE=$PRE ALFP7_PROF=1 ALFDEVQR_PROF=1 ALFSH_TIMING=1
r32 N_V7_cs3 V7 3 P7PRE=$PRE ALFP7_PROF=1 ALFDEVQR_PROF=1 ALFSH_TIMING=1; tm N_V7_cs1; tm N_V7_cs3; pg N_V7_cs3
mkdir -p $P7/data/j5; for d in $W/runs/*; do t=$(basename $d); mkdir -p $P7/data/j5/$t; cp $d/rank*.log $d/info $P7/data/j5/$t/ 2>/dev/null; done
echo "end $(date -Is)"
