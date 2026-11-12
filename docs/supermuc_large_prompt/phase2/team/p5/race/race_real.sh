#!/bin/bash
#SBATCH --job-name=rc_real
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
#SBATCH --time=01:15:00
# usage: sbatch race_real.sh <N> [sec_real=1800] [sec_ctl=600]: real kernels (DEVQR=2 + DEVUNMQR=1), then controls oneMKL and DEVQR=1, same node/job
N=${1:-1024}; SEC=${2:-1800}; SC=${3:-600}
source $SCRATCH/large_prompt/team/p5/race/race_lib.sh
runcase real_N$N $N $SEC ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1 ALFSH_DEVQR_LIB=$LIB
runcase mkl_N$N $N $SC ALFSH_DEVQR=0 ALFSH_DEVUNMQR=0
runcase devqr1_N$N $N $SC ALFSH_DEVQR=1 ALFSH_DEVUNMQR=0 ALFSH_DEVQR_LIB=$LIB
echo "end $(date -Is)"
