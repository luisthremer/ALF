#!/bin/bash
#SBATCH --job-name=rc_plant2
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
# usage: sbatch race_plant2.sh <plant> [N=576] [sec per case=900]: default panel shape, then the amplifier ALFDEVQR_G=1 (existing runtime knob of the
# library: every panel uses 1 sub-group per column, i.e. the barrier-free-handoff variant of the panel kernel for all M)
PL=${1:-race_tau}; N=${2:-576}; SEC=${3:-900}
source $SCRATCH/large_prompt/team/p5/race/race_lib.sh
runcase plant_${PL}_N$N $N $SEC ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1 ALFSH_DEVQR_LIB=$P/plant/$PL/libalfdevqr.so
runcase plant_${PL}_N${N}_G1 $N $SEC ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1 ALFSH_DEVQR_LIB=$P/plant/$PL/libalfdevqr.so ALFDEVQR_G=1
echo "end $(date -Is)"
