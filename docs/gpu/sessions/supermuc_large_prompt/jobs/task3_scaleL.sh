#!/bin/bash
#SBATCH --job-name=lp_sizes
#SBATCH --output=lp_sizes.%j.log
#SBATCH --error=lp_sizes.%j.err
#SBATCH --time=00:55:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# Task 3 (test only): time per sweep at L=28 and L=32 (Beta=L, N_Phi=L^2/2, other parameters = reference), 32 ranks = 4 per tile,
# ZEX_NUMBER_OF_CCS=:4 per tile. usage: sbatch task3_scaleL.sh "28 32" [nranks]
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt; J=$L0/jobs; cd $J
NR=${2:-32}
export SEEDS=$L0/ref_L24/seeds OMP_TARGET_OFFLOAD=MANDATORY ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4
G=${GPUBIN:-$L0/build_mpi/Prog/ALF.out}
echo "host $(hostname) $(date -Is) ranks=$NR bin=$G"
for L in $1; do
  export W=$L0/work_L$L; mkdir -p $W
  sed -e "s/L1 *= *24/L1 = $L/" -e "s/L2 *= *24/L2 = $L/" -e "s/N_Phi *= *288/N_Phi = $((L*L/2))/" -e "s/Beta *= *24.d0/Beta = $L.d0/" $L0/ref_L24/parameters > $W/params_L$L
  export PARAMS=$W/params_L$L
  grep -E "L1|N_Phi|Beta" $PARAMS | tr '\n' ' '; echo
  python3 run_ref_sng.py L${L}_s1 $G $NR 1 1 | tail -1
  python3 run_ref_sng.py L${L}_s3 $G $NR 3 1 | tail -1
  T1=$(grep -h "CPU Time" $W/runs/L${L}_s1/info | awk '{print $NF}'); T3=$(grep -h "CPU Time" $W/runs/L${L}_s3/info | awk '{print $NF}')
  python3 - <<PY
t1=float("$T1"); t3=float("$T3"); R=$NR
m=(t3-t1)/2
print(f"RESULT L=$L ranks={R} T(1)={t1:.1f}s T(3)={t3:.1f}s marginal={m:.2f} s/sweep/chain aggregate={R/m:.3f} sweeps/s")
PY
  grep -h -E "Precision Green|Acceptance" $W/runs/L${L}_s3/info | head -3
  (xpu-smi dump -d -1 -m 18 -n 1 2>&1 | grep -v lspci | tail -5) 
done
echo "done $(date -Is)"
