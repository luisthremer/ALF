#!/bin/bash
#SBATCH --job-name=lp_scanL
#SBATCH --output=lp_scanL.%j.log
#SBATCH --error=lp_scanL.%j.err
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
# System-size scan: s/sweep/chain and aggregate throughput for L=28 and 32 (beta=L, U=5.3, N_Phi=L^2/2, ZEX_NUMBER_OF_CCS=4 per tile)
# usage: sbatch scanL.sh "28 32" "8 32"
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
L=$SCRATCH/large_prompt; J=$L/jobs
export W=$L/work_scanL SEEDS=$L/ref_L24/seeds PARAMS=$L/ref_L24/parameters; mkdir -p $W; cd $J
G=${GPUBIN:-$L/build_mpi/Prog/ALF.out}
export OMP_TARGET_OFFLOAD=MANDATORY ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4
echo "host $(hostname) $(date -Is) bin=$G"
for S in $1; do
  ND=32; [ $S -ge 32 ] && ND=64
  for R in $2; do
    echo "=== L=$S N_delay=$ND ranks=$R"
    LSIZE=$S NDELAY=$ND python3 run_ref_sng.py sl_${S}_${R}_s1 $G $R 1 1 | tail -1
    LSIZE=$S NDELAY=$ND python3 run_ref_sng.py sl_${S}_${R}_s3 $G $R 3 1 | tail -1
    T1=$(grep -h "CPU Time" $W/runs/sl_${S}_${R}_s1/info | awk '{print $NF}'); T3=$(grep -h "CPU Time" $W/runs/sl_${S}_${R}_s3/info | awk '{print $NF}')
    python3 - <<PY
t1=float("$T1"); t3=float("$T3"); R=$R
m=(t3-t1)/2
print(f"RESULT L=$S ranks={R} T(1)={t1:.1f}s T(3)={t3:.1f}s marginal={m:.2f} s/sweep/chain aggregate={R/m:.3f} sweeps/s")
PY
  done
done
echo "done $(date -Is)"
