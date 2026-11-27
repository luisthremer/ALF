#!/bin/bash
# Gate (gate.py) with the Intel switches ON in the GPU lanes: ALFSH_ASYNCBLK=1 ALFSH_BLKFUSE=5 ALFSH_BLKPF=1 ALFSH_OBSCORR=1;
# the default build (intel_* cases of gate_cases.sh, incl. delayed: the block switches act on submatrix only) and STABLOG (log_*),
# plus the negative controls. Same criteria as job_gate.sbatch.
#SBATCH -J p2i_gon
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --time=00:59:00
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH -o %x.%j.out
module load slurm_setup
S=$HOME/p2_intel/scripts; B=$HOME/p2_intel/build; R=$SCRATCH/p2_intel/gon_$SLURM_JOB_ID
. $S/env_smuc2.sh
export OMP_NUM_THREADS=14
mkdir -p $R
ON="env:ALFSH_ASYNCBLK=1 env:ALFSH_BLKFUSE=5 env:ALFSH_BLKPF=1 env:ALFSH_OBSCORR=1"
bash $S/gate_cases.sh $B | grep -E "^(intel_|log_|negctl_)" | sed "s|\$| $ON|" > $R/cases.txt
python3 $S/gate_driver.py $R/cases.txt $R | sort -k2 > $R/gate.txt
cut -c1-200 $R/gate.txt
echo "PASS $(grep -c ^PASS $R/gate.txt) FAIL $(grep -c ^FAIL $R/gate.txt) of $(wc -l < $R/cases.txt) (negctl_* must FAIL)"
grep -h "Intel block\|OBSCORR\|prefetch" $R/intel_mzF_L16_lt1_submatrix/gpu/info
