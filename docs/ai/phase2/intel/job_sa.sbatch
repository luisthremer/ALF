#!/bin/bash
# S-a on SuperMUC-NG Phase 2: device ctest (13-linalg, GPU tests, e2e) and first trajectory checks on one tile.
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
( cd $B/src_intel && rm -rf testsuite/tests && mkdir -p testsuite/tests && cd testsuite/tests && cmake .. > $R/cmake.log 2>&1 \
  && make -j16 > $R/tmake.log 2>&1 && ctest -R "gpu|linalg|det-phase" -j2 --timeout 3000 --output-on-failure > $R/ctest.log 2>&1; \
  tail -3 $R/ctest.log )
ST=$B/src_none/Scripts_and_Parameters_files/Start
for L in 6 12; do for mz in .F. .T.; do
  echo "== L=$L Mz=$mz"
  python3 $S/traj.py $R/traj_L${L}_mz$mz $B/ALF-intel.out $B/ALF-none.out $ST 1e-8 L1=$L L2=$L Mz=$mz NSweep=20 NBin=5
done; done
