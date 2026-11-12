#!/bin/bash
#SBATCH --job-name=rc_plant
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
#SBATCH --time=00:42:00
# usage: sbatch race_plant.sh <plant: race_tau|race_ex|race_v|onebuf|nobar|tau|conj> [N=1024] [sec=1500]
PL=${1:-race_tau}; N=${2:-1024}; SEC=${3:-1500}
source $SCRATCH/large_prompt/team/p5/race/race_lib.sh
runcase plant_${PL}_N$N $N $SEC ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1 ALFSH_DEVQR_LIB=$P/plant/$PL/libalfdevqr.so
echo "end $(date -Is)"
