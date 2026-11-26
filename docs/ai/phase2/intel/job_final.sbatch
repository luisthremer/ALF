#!/bin/bash
# Re-check of the final (rebased) canvas on one node: device ctest (JIT test tree) and the gate with the default build
# (gate_cases.sh lines intel_* and negctl_*). Builds: build_smuc2.sh + the JIT test tree, from the same source.
#SBATCH -J p2i_fin
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --time=00:59:00
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH -o %x.%j.out
module load slurm_setup
S=$HOME/p2_intel/scripts; B=$HOME/p2_intel/build; R=$SCRATCH/p2_intel/fin_$SLURM_JOB_ID
. $S/env_smuc2.sh
module load cmake/3.30.0 >/dev/null 2>&1
export OMP_NUM_THREADS=14
mkdir -p $R
( export ZE_AFFINITY_MASK=0; cd $B/src_test/testsuite/tests && ctest -R "gpu|linalg|det-phase" -j2 --timeout 2400 --output-on-failure > $R/ctest.log 2>&1; tail -3 $R/ctest.log ) &
bash $S/gate_cases.sh $B | grep -E "^(intel_|negctl_)" > $R/cases.txt
python3 $S/gate_driver.py $R/cases.txt $R | sort -k2 > $R/gate.txt
wait
cut -c1-220 $R/gate.txt
echo "PASS $(grep -c ^PASS $R/gate.txt) FAIL $(grep -c ^FAIL $R/gate.txt) of $(wc -l < $R/cases.txt) (negctl_* must FAIL)"
