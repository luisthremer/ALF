#!/bin/bash
# S-a on SuperMUC-NG Phase 2: first trajectory checks on one tile (device ctest: job_ctest.sbatch).
# usage: sbatch job_sa.sbatch   (from $HOME/p2_intel/scripts; builds in $HOME/p2_intel/build by build_smuc2.sh)
#SBATCH -J p2i_sa
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --time=00:59:00
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH -o %x.%j.out
module load slurm_setup
S=$HOME/p2_intel/scripts; B=$HOME/p2_intel/build; R=$SCRATCH/p2_intel/sa_$SLURM_JOB_ID
. $S/env_smuc2.sh
export ZE_AFFINITY_MASK=0
mkdir -p $R
ST=$B/src_none/Scripts_and_Parameters_files/Start
for L in 6 12; do for mz in .F. .T.; do
  echo "== L=$L Mz=$mz"
  python3 $S/traj.py $R/traj_L${L}_mz$mz $B/ALF-intel.out $B/ALF-none.out $ST 1e-8 L1=$L L2=$L Mz=$mz NSweep=20 NBin=5
done; done
