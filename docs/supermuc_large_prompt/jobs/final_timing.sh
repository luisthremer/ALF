#!/bin/bash
#SBATCH --job-name=lp_final
#SBATCH --output=lp_final.%j.log
#SBATCH --error=lp_final.%j.err
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# Clean marginal s/sweep/chain of the TUNED production configuration (32 ranks = 4 per tile, ZEX CCS, 4 host threads, Nwrap 20, N_delay 64,
# no pin domain), L = $1, from runs of 1 and 4 sweeps (HDF5 PP-independent binary build_mpi).
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt; J=$L0/jobs; cd $J; LL=$1
export SEEDS=$L0/ref_L24/seeds OMP_TARGET_OFFLOAD=MANDATORY ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4 OMPT=4
export W=$L0/work_final_L$LL; mkdir -p $W
sed -e "s/L1 *= *24/L1 = $LL/" -e "s/L2 *= *24/L2 = $LL/" -e "s/N_Phi *= *288/N_Phi = $((LL*LL/2))/" -e "s/Beta *= *24.d0/Beta = $LL.d0/" -e "s/Nwrap *= *10/Nwrap = 20/" -e "s/N_delay *= *32/N_delay = 64/" $L0/ref_L24/parameters > $W/params; export PARAMS=$W/params
echo "host $(hostname) $(date -Is)"
python3 run_ref_sng.py f1 $L0/build_mpi/Prog/ALF.out 32 1 1 | tail -1
python3 run_ref_sng.py f4 $L0/build_mpi/Prog/ALF.out 32 4 1 | tail -1
python3 - <<PY
import re,subprocess
def cpu(t): return float(re.search(r"CPU Time\s*:\s*([\d.]+)", open("$W/runs/%s/info"%t).read()).group(1))
w=lambda l: float(re.search(r"wall=([\d.]+)s", l).group(1))
print("RESULT L=$LL (CPU-time based) T1=%.1f T4=%.1f"%(cpu("f1"),cpu("f4")))
PY
echo "done $(date -Is)"
