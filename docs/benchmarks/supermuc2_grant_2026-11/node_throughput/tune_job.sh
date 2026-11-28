#!/bin/bash
#SBATCH --job-name=tune_all
#SBATCH --output=/hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27/auto/tune_all.%j.log
#SBATCH --error=/hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27/auto/tune_all.%j.err
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=general
#SBATCH --nodes=32
#SBATCH --time=08:00:00
#SBATCH --exclude=i20r02c06s02,i20r02c06s06,i20r03c06s06
# Everything in one allocation: the automatic search (tune.py search) and then the beta = L runs with the tuned settings
# (tune.py final --beta L). Each run is one srun step on one full node; steps wait for free nodes. Finished runs are reused (cache).
module load slurm_setup >/dev/null 2>&1
cd /hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27/jobs/nt
python3 tune.py search --inalloc
python3 tune.py final --beta L --inalloc
