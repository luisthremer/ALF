#!/bin/bash
#SBATCH --job-name=qf_build
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH -x i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# ALF builds from the git-archive tree of the worktree HEAD (p8/qrfix/src_head), same recipe as p5/c_build5.sh:
# aottt (AOT for PVC + phase timers; PG gates) and txtnew (text observables: same-seed gates)
echo "start $(date -Is) node $(hostname) HEAD $(cat $SCRATCH/large_prompt/team/p8/qrfix/src_head/GIT_HEAD)"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P=$SCRATCH/large_prompt/team/p8/qrfix; cd $P; export SRC=$P/src_head
BLD=$P/build_aottt bash $SCRATCH/large_prompt/team/p5/jobs/build_comb.sh aottt > b_aottt.out 2>&1 &
BLD=$P/build_txtnew bash $SCRATCH/large_prompt/team/p5/jobs/build_comb.sh txtnew > b_txtnew.out 2>&1 &
wait
ls -l $P/build_aottt/Prog/ALF.out $P/build_txtnew/Prog/ALF.out $P/build_*/Libraries/Modules/lib*.so
tail -3 b_aottt.out b_txtnew.out
md5sum $P/build_aottt/Libraries/Modules/libalfdevqr.so $SCRATCH/large_prompt/team/p5/build_aottt/Libraries/Modules/libalfdevqr.so
echo "end $(date -Is)"
