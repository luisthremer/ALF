#!/bin/bash
#SBATCH --job-name=lp_vrep
#SBATCH --output=lp_vrep.%j.log
#SBATCH --error=lp_vrep.%j.err
#SBATCH --time=00:10:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
module load intel-oneapi-vtune/2024.2.1
W=$SCRATCH/large_prompt/work_vtune2
vtune -report top-down -result-dir $W/r_hot -call-stack-mode all -column="CPU Time:Self,CPU Time:Total" -limit 3000 > $W/topdown.txt 2>&1; wc -l $W/topdown.txt
