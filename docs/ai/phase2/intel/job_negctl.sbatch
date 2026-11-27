#!/bin/bash
# Negative controls of the gate only (gate_cases.sh lines negctl_*): both must FAIL.
#SBATCH -J p2i_neg
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --time=00:15:00
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH -o %x.%j.out
module load slurm_setup
S=$HOME/p2_intel/scripts; B=$HOME/p2_intel/build; R=$SCRATCH/p2_intel/neg_$SLURM_JOB_ID
. $S/env_smuc2.sh
export OMP_NUM_THREADS=14
mkdir -p $R
bash $S/gate_cases.sh $B | grep ^negctl > $R/cases.txt
python3 $S/gate_driver.py $R/cases.txt $R
