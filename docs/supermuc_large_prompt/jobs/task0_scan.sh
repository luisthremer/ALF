#!/bin/bash
#SBATCH --job-name=lp_scan
#SBATCH --output=lp_scan.%j.log
#SBATCH --error=lp_scan.%j.err
#SBATCH --time=00:55:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# Task 0: oversubscription scan, L=24 reference parameters, ONE node. usage: sbatch task0_scan.sh "8 16 32" [ccs]
# Device of rank r = r mod 8 (ALF: omp_set_default_device(rank mod #devices); ZE_AFFINITY_MASK unset, ZE_FLAT_DEVICE_HIERARCHY=FLAT)
# Aggregate throughput from the marginal cost of 3 extra sweeps (runs of 1 and 4 sweeps per rank, one bin each), so start-up cost cancels.
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
L=$SCRATCH/large_prompt; J=$L/jobs
export W=$L/work_scan SEEDS=$L/ref_L24/seeds PARAMS=$L/ref_L24/parameters; mkdir -p $W; cd $J
G=${GPUBIN:-$L/build_mpi/Prog/ALF.out}   # HDF5 build: the text-output builds cost ~70 s per bin (measured), HDF5 ~1 s
export OMP_TARGET_OFFLOAD=MANDATORY
echo "host $(hostname) $(date -Is)  ZE_FLAT_DEVICE_HIERARCHY=${ZE_FLAT_DEVICE_HIERARCHY:-unset} ZE_AFFINITY_MASK=${ZE_AFFINITY_MASK:-unset} bin=$G"
if [ "$2" = ccs ]; then export ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4; echo "ZEX_NUMBER_OF_CCS=$ZEX_NUMBER_OF_CCS"; fi
cpu_time() { grep -h "CPU Time" $W/runs/$1/info | awk '{print $NF}'; }
for R in $1; do
  echo "=== $R ranks (rank mod 8 -> $((R/8)) chain(s) per tile)"
  ( sleep 45; echo "-- xpu-smi ps"; xpu-smi ps 2>&1 | grep -v lspci | head -50 ) > $W/ps_$R.txt 2>&1 &
  python3 run_ref_sng.py scan_${R}_s1 $G $R 1 1 | tail -1     # 1 sweep, 1 bin
  python3 run_ref_sng.py scan_${R}_s4 $G $R 4 1 | tail -1     # 4 sweeps in the same single bin
  wait
  T1=$(cpu_time scan_${R}_s1); T4=$(cpu_time scan_${R}_s4)
  python3 - <<PY
t1=float("$T1"); t4=float("$T4"); R=$R
m=(t4-t1)/3
print(f"RESULT ranks={R} T(1 sweep)={t1:.1f}s T(4 sweeps)={t4:.1f}s  marginal={m:.2f} s/sweep/chain  aggregate={R/m:.3f} sweeps/s  per-tile={R/m/8:.3f} sweeps/s")
PY
done
echo "done $(date -Is)"
