#!/bin/bash
#SBATCH --job-name=lp_tune
#SBATCH --output=lp_tune.%j.log
#SBATCH --error=lp_tune.%j.err
#SBATCH --time=00:55:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
# Task 1: single rank on one tile, L=24 reference parameters, HDF5 binary; runtime-knob sweep + offload plugin profile.
# usage: sbatch tune1.sh [binary]   (CPU time of 1 bin with 4 sweeps, start-up ~ 10 s)
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
L=$SCRATCH/large_prompt; J=$L/jobs
export W=$L/work_tune SEEDS=$L/ref_L24/seeds PARAMS=$L/ref_L24/parameters; mkdir -p $W; cd $J
G=${1:-$L/build_mpi/Prog/ALF.out}
export OMP_TARGET_OFFLOAD=MANDATORY
echo "host $(hostname) $(date -Is) bin=$G"
run() { # tag, env assignments...
  tag=$1; shift
  env "$@" python3 run_ref_sng.py tune_$tag $G 1 4 1 | tail -1
  echo "   CPU Time $(grep -h 'CPU Time' $W/runs/tune_$tag/info | awk '{print $NF}') s   [$*]"
}
run base            X=1
run immcl_all       LIBOMPTARGET_LEVEL_ZERO_USE_IMMEDIATE_COMMAND_LIST=all
run immcl_compute   LIBOMPTARGET_LEVEL_ZERO_USE_IMMEDIATE_COMMAND_LIST=compute
run sycl_imm        SYCL_PI_LEVEL_ZERO_USE_IMMEDIATE_COMMANDLISTS=1 UR_L0_USE_IMMEDIATE_COMMANDLISTS=1
run profile         LIBOMPTARGET_PLUGIN_PROFILE=T
echo "=== plugin profile (rank 0), head"
grep -v "^ *$" $W/runs/tune_profile/run.log | head -90 | cut -c1-170
echo "done $(date -Is)"
