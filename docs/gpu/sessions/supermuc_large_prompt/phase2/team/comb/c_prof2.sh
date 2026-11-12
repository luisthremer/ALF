#!/bin/bash
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:59:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; C=$L0/comb; cd $C/jobs
# Phase 5 profile, clean per-call LAPACK times: V3 (devqr, NO async block path -> every shim call starts on an idle queue)
# and V1 (oneMKL geqrf) at L24/L32 one tile, ALFSH_TIMING=1. Runs: NSWEEP 4 vs 2 to subtract setup per call count.
source c_lib.sh; export W=$C/work_prof2
for v in V3 V1; do for L in 24 32; do for ns in 2 4; do
  t=P_L${L}_${v}_s$ns; if [ $L = 32 ]; then one $t $v $ns LSIZE=32 ALFSH_TIMING=1; else one $t $v $ns ALFSH_TIMING=1; fi
  grep -h -E "ALFSH_TIMING (geqrf|geqrf_dev|unmqr|getrf) " $W/runs/$t/run.log
done; done; done
echo "done $(date -Is)"
