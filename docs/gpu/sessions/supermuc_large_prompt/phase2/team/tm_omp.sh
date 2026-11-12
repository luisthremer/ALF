#!/bin/bash
#SBATCH --job-name=tm_omp
#SBATCH --output=tm_omp.%j.log
#SBATCH --error=tm_omp.%j.err
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# Device-side profile of the non-LAPACK part: Intel libomptarget LIBOMPTARGET_PROFILE (kernel counts, device time, transfers),
# L=24 one tile, 8 host threads, 2 and 6 sweeps (marginal = difference), shim timers on.
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; cd $J
export W=$L0/work_omp SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 ALFSH_TIMING=1
export LIBOMPTARGET_PROFILE=T,usec
B=$L0/build_aottim/Prog/ALF.out
for ns in 2 6; do OMPT=8 python3 run_ref_sng.py p_s$ns $B 1 $ns 1 | tail -1; done
for ns in 2 6; do echo "=== s$ns"; grep -c . $W/runs/p_s$ns/run.log; done
echo "done $(date -Is)"
