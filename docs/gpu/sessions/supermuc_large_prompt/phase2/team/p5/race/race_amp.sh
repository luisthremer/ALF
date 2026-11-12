#!/bin/bash
#SBATCH --job-name=rc_amp
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH -x i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:20:00
# real kernels with the G=1 amplifier (ALFDEVQR_G=1, valid for N<=640 only): the amplifier exposed the race_v plant at 530 BAD / 713k iterations
# (default panel shape: 1 / 851k). usage: sbatch race_amp.sh [sec=600]
SEC=${1:-600}
source $SCRATCH/large_prompt/team/p5/race/race_lib.sh
runcase real_N576_G1 576 $SEC ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1 ALFSH_DEVQR_LIB=$LIB ALFDEVQR_G=1
echo "end $(date -Is)"
