#!/bin/bash
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:59:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; C=$L0/comb; cd $C/jobs
# combined branch: build the timer (AOT) and text builds in parallel
bash build_comb.sh aottt > b_aottt.out 2>&1 & bash build_comb.sh txtnew > b_txtnew.out 2>&1 & wait
ls -l $C/build_aottt/Prog/ALF.out $C/build_txtnew/Prog/ALF.out $C/build_*/Libraries/Modules/lib*.so || exit 1
