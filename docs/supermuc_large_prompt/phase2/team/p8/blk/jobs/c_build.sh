#!/bin/bash
#SBATCH --job-name=p8b_build
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# p8/blk: aottxt (AOT for PVC + phase timers + text observables: gates G1-G4) and txtnew (JIT text, same flags as the p5 txtnew: default-vs-previous gate)
# sources: $SCRATCH/large_prompt/team/p8/blk/src_alf = git archive HEAD + the 2 changed files
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P=$SCRATCH/large_prompt/team/p8/blk; cd $P/jobs; export SRC=$P/src_alf
echo "start $(date -Is) $(hostname)"
BLD=$P/build_aottxt bash build_comb.sh aottxt > b_aottxt.out 2>&1 &
BLD=$P/build_txtnew bash build_comb.sh txtnew > b_txtnew.out 2>&1 &
wait
ls -l $P/build_aottxt/Prog/ALF.out $P/build_txtnew/Prog/ALF.out $P/build_*/Libraries/Modules/lib*.so
tail -3 b_aottxt.out b_txtnew.out
echo "end $(date -Is)"
