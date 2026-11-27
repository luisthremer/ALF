#!/bin/bash
# Device ctest on SuperMUC-NG Phase 2 (one tile): 13-linalg, the GPU tests and the e2e comparison.
# The test tree is a JIT build ($HOME/p2_intel/build/src_test, configure.sh SUPERMUCNG2 noMPI GPU=intel without
# ALF_INTEL_AOT): the testsuite CMake drops the repeated -Xopenmp-target-backend of the AOT flags.
#SBATCH -J p2i_ct
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --time=00:45:00
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH -o %x.%j.out
module load slurm_setup
S=$HOME/p2_intel/scripts; R=$SCRATCH/p2_intel/ct_$SLURM_JOB_ID
. $S/env_smuc2.sh
module load cmake/3.30.0 >/dev/null 2>&1
export ZE_AFFINITY_MASK=0
mkdir -p $R
cd $HOME/p2_intel/build/src_test/testsuite/tests && ctest -R "gpu|linalg|det-phase" -j2 --timeout 2400 --output-on-failure > $R/ctest.log 2>&1
tail -4 $R/ctest.log
ctest -R gpu-state-e2e -V 2>&1 | grep -E "OK|FAIL|vs" > $R/e2e.txt; head -40 $R/e2e.txt
