#!/bin/bash -l
# One-node Slurm job running testsuite/lowrank/cluster/validate.sh (all groups concurrently).
# Submit from the ALF directory, e.g. on SuperMUC-NG:
#   sbatch --account=<project> --partition=test --export=ALL,ALF_MACHINE=SuperMUC-NG,ALF_ENV_SH=<workdir>/env_alf.sh \
#          testsuite/lowrank/cluster/validate.sbatch [GROUP ...]
# The binaries must have been built on the login node first (build_variants.sh, parallel_tier.sh --build-only).
#SBATCH -J alf_validate
#SBATCH -o runs/cluster_val_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=00:30:00
#SBATCH --get-user-env
module load slurm_setup >/dev/null 2>&1
echo "job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ), groups: ${*:-all}"
testsuite/lowrank/cluster/validate.sh "$@"
echo "validate.sh exit $?"
