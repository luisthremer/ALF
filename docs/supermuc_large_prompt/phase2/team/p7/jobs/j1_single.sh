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
#SBATCH --time=00:50:00
#SBATCH --ntasks-per-node=8
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
source $SCRATCH/large_prompt/team/p7/jobs/p7lib.sh
export W=$P7/work_j1; mkdir -p $W $P7/lib $P7/data
# --- build the instrumented block library (LD_PRELOAD) ---
cd $P7/lib && icpx -fsycl -fPIC -shared -O2 -std=c++17 -fsycl-targets=spir64_gen -Xs "-device pvc" -I$MKLROOT/include $P7/src/alf_la_sycl_blk_p7.cpp -o libalfsb.so -qmkl -lmkl_sycl_blas 2>&1 | tail -5
ls -la $P7/lib; ldd $BIN7 | grep -i alfsb
cd $PJ7
PRE=$P7/lib/libalfsb.so
for v in V7 V2; do
  echo "=== $v unperturbed timing (3 sweeps, ALFSH_TIMING)"; r1 ${v}_n3 $v 3 ALFSH_TIMING=1; tm ${v}_n3; grep -h "ALFSH_TIMING\|CPU Time" $W/runs/${v}_n3/rank0.log | head -14
  echo "=== $v PLUGIN_PROFILE 1 and 3 sweeps"; r1 ${v}_pp1 $v 1 LIBOMPTARGET_PLUGIN_PROFILE=T; r1 ${v}_pp3 $v 3 LIBOMPTARGET_PLUGIN_PROFILE=T; tm ${v}_pp1; tm ${v}_pp3
  echo "=== $v full profile (PLUGIN_PROFILE + serialised block path + devqr stage profile + shim timers) 1 and 3 sweeps"
  r1 ${v}_pf1 $v 1 LIBOMPTARGET_PLUGIN_PROFILE=T P7PRE=$PRE ALFP7_PROF=1 ALFDEVQR_PROF=1 ALFSH_TIMING=1
  r1 ${v}_pf3 $v 3 LIBOMPTARGET_PLUGIN_PROFILE=T P7PRE=$PRE ALFP7_PROF=1 ALFDEVQR_PROF=1 ALFSH_TIMING=1; tm ${v}_pf1; tm ${v}_pf3
  grep -h "ALFP7" $W/runs/${v}_pf3/rank0.log | head -20
done
echo "=== V7 MKL_VERBOSE 1 and 3 sweeps"; r1 V7_mk1 V7 1 MKL_VERBOSE=1; r1 V7_mk3 V7 3 MKL_VERBOSE=1
echo "=== V7 LIBOMPTARGET_INFO=-1 (launch sizes), 1 sweep"; r1 V7_info V7 1 LIBOMPTARGET_INFO=-1
python3 $PJ7/parse_info.py $W/runs/V7_info/rank0.log > $W/kinfo_V7.txt; wc -l $W/kinfo_V7.txt; rm -f $W/runs/V7_info/rank0.log
echo "=== standard / delayed / submatrix N_delay (V7, tile 0, 2 sweeps, same seeds)"
for c in "standard 32" "delayed 32" "submatrix 16" "submatrix 32" "submatrix 64"; do set -- $c
  SCHEME=$1 NDELAY=$2 timeout 900 bash -c "source $PJ7/p7lib.sh; export W=$W; r1 sch_$1_$2 V7 2" ; echo "SCHEME $1 N_delay $2:"; tm sch_$1_$2; grep -h "Precision Green" $W/runs/sch_$1_$2/info
done
mkdir -p $P7/data/j1; for d in $W/runs/*; do t=$(basename $d); mkdir -p $P7/data/j1/$t; cp $d/rank0.log $d/info $P7/data/j1/$t/ 2>/dev/null; done
echo "end $(date -Is)"
