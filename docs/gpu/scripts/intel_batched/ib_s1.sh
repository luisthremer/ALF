#!/bin/bash -l
#SBATCH -J ib_s1
#SBATCH -A pn93la
#SBATCH -p preempt
#SBATCH --qos=preempt
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH -t 00:25:00
#SBATCH -x i20r02c06s02
#SBATCH -o %x.%j.log
#SBATCH -e %x.%j.err
# Step 1 of the Intel batched study: bandwidth / gemm / batch / launch latency on one tile, alone and shared.
D=$SCRATCH/large_prompt/batched
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
B=$D/bin/ib_bench
O=$D/out/$SLURM_JOB_ID; mkdir -p $O
echo "start $(date -Is) $(hostname)"
export ZE_FLAT_DEVICE_HIERARCHY=FLAT
T=${TMIN:-0.4}
# 0. bitwise check of the multi-walker kernels (alf_la_sycl_blk.cpp *_w) + planted defect (must FAIL)
ZE_AFFINITY_MASK=0 $D/bin/ib_wtest > $O/wtest.txt 2>&1; echo "wtest rc=$?" >> $O/wtest.txt
ZE_AFFINITY_MASK=0 IB_PLANT=1 $D/bin/ib_wtest > $O/wtest_plant.txt 2>&1; echo "wtest_plant rc=$?" >> $O/wtest_plant.txt
# 1. one process on tile 0 (other tiles idle)
ZE_AFFINITY_MASK=0 $B all $T > $O/r1.txt 2>&1
# 2. launch latency with/without immediate command lists
for ic in 0 1; do
  ZE_AFFINITY_MASK=0 SYCL_PI_LEVEL_ZERO_USE_IMMEDIATE_COMMANDLISTS=$ic SYCL_UR_L0_USE_IMMEDIATE_COMMANDLISTS=$ic $B launch > $O/launch_ic$ic.txt 2>&1
done
# 3. OpenMP target latency (ifx)
ZE_AFFINITY_MASK=0 OMP_NUM_THREADS=4 $D/bin/ib_omp > $O/omp.txt 2>&1
ZE_AFFINITY_MASK=0 OMP_NUM_THREADS=4 LIBOMPTARGET_LEVEL_ZERO_COMMAND_MODE=async $D/bin/ib_omp > $O/omp_async.txt 2>&1
# 4. P processes sharing tile 0 (no CCS), and with CCS 4 (process gets one compute slice of the tile)
for P in 2 4; do
  for ccs in 1 4; do
    pids=""
    for p in $(seq 1 $P); do
      ZE_AFFINITY_MASK=0 ZEX_NUMBER_OF_CCS=0:$ccs taskset -c $((p*8))-$((p*8+7)) $B ${MODE2:-gemm} $T > $O/p${P}_ccs${ccs}_$p.txt 2>&1 &
      pids="$pids $!"
    done
    wait $pids
  done
done
# 5. block emulation with P = 4 processes per tile, CCS 4 (production layout) for comparison with batched W in one process
pids=""
for p in 1 2 3 4; do
  ZE_AFFINITY_MASK=0 ZEX_NUMBER_OF_CCS=0:4 taskset -c $((p*8))-$((p*8+7)) $B block $T > $O/blk4_ccs4_$p.txt 2>&1 &
  pids="$pids $!"
done
wait $pids
# 6. all 8 tiles busy, one process each (node loaded; tile 0 reported)
pids=""
for t in 0 1 2 3 4 5 6 7; do
  ZE_AFFINITY_MASK=$t taskset -c $((t*14))-$((t*14+13)) $B gemm $T > $O/n8_t$t.txt 2>&1 &
  pids="$pids $!"
done
wait $pids
echo "end $(date -Is)"
