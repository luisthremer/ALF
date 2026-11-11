#!/bin/bash
#SBATCH --job-name=tm_ball
#SBATCH --output=tm_ball.%j.log
#SBATCH --error=tm_ball.%j.err
#SBATCH --time=00:50:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=224
# Team session: the four builds in parallel on one node (aot = timing binary, gpu-txt/cpu-txt = same-seed comparison,
# gpu-ser = repo ctest on the GPU). usage: sbatch build_all.sh [variants...]
J=$SCRATCH/large_prompt/team/jobs; cd $J
VS=${@:-aot gpu-txt cpu-txt gpu-ser}
for v in $VS; do
  if [ $v = gpu-ser ]; then bash $J/build_intel.sh $v ctest > b_$v.out 2>&1 & else bash $J/build_intel.sh $v > b_$v.out 2>&1 & fi
done
wait; date -Is
for v in $VS; do echo "== $v"; ls -l $SCRATCH/large_prompt/team/build_$v/Prog/ALF.out; tail -4 b_$v.out; done
