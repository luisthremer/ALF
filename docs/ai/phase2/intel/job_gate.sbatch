#!/bin/bash
# Device-path gate on one node (8 tiles): gate_cases.sh -> gate_driver.py -> gate.py (see gate.py for the criteria).
#SBATCH -J p2i_gate
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --time=00:59:00
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH -o %x.%j.out
module load slurm_setup
S=$HOME/p2_intel/scripts; B=$HOME/p2_intel/build; R=$SCRATCH/p2_intel/gate_$SLURM_JOB_ID
. $S/env_smuc2.sh
export OMP_NUM_THREADS=14
mkdir -p $R
bash $S/gate_cases.sh $B > $R/cases.txt
python3 $S/gate_driver.py $R/cases.txt $R | sort -k2 | tee $R/gate.txt
echo "PASS $(grep -c ^PASS $R/gate.txt) FAIL $(grep -c ^FAIL $R/gate.txt) of $(wc -l < $R/cases.txt) (negctl_* must FAIL)"
