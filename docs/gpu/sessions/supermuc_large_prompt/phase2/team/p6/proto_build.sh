#!/bin/bash
#SBATCH --job-name=p6_pbuild
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:59:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH -x i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# Phase 6 prototype builds from $P6/src_alf (= p5/src_alf + the three prototype files): aottt (AOT, phase timers, HDF5: timing) and txtnew (text: same-seed / bitwise gates)
echo "start $(date -Is) node $(hostname)"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P6=$SCRATCH/large_prompt/team/p6; P5=$SCRATCH/large_prompt/team/p5
cd $P6; export SRC=$P6/src_alf
BLD=$P6/build_aottt bash $P5/jobs/build_comb.sh aottt > b_aottt.out 2>&1 &
BLD=$P6/build_txtnew bash $P5/jobs/build_comb.sh txtnew > b_txtnew.out 2>&1 &
wait
ls -l $P6/build_aottt/Prog/ALF.out $P6/build_txtnew/Prog/ALF.out $P6/build_*/Libraries/Modules/lib*.so
tail -3 b_aottt.out b_txtnew.out
echo "end $(date -Is)"
