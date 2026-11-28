#!/bin/bash
#SBATCH --job-name=tune_fin
#SBATCH --output=/hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27/auto/tune_fin.%x.%j.log
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=3
#SBATCH --time=01:00:00
# beta = L finals for L = 48, 52, 56 with the extended layout fallback (2 and 1 walkers per tile). usage: sbatch tune_final_job.sh SIDE
module load slurm_setup >/dev/null 2>&1
cd /hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27/jobs/nt
python3 tune.py final --beta L --Ls 48,52,56 --sides $1 --inalloc
