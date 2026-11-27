#!/bin/bash
# Layout question (2.5x vs gpu-port's documented L = 24 numbers): host threads and pinning. L = 24 reference set
# (Beta 24, N_FL 2, U 5.3, h0 1, submatrix N_delay 32), s/sweep = CPU(NBin 2) - CPU(NBin 1) over 3 sweeps.
#SBATCH -J p2i_lay
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --time=00:45:00
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH -o %x.%j.out
module load slurm_setup
S=$HOME/p2_intel/scripts; B=$HOME/p2_intel/build; R=$SCRATCH/p2_intel/lay_$SLURM_JOB_ID
. $S/env_smuc2.sh
mkdir -p $R
ST=$B/src_none/Scripts_and_Parameters_files/Start
P="L1=24 L2=24 Beta=24.d0 N_Phi=288 N_FL=2 Mz=.F. ham_h0=1.d0 ham_U=5.3d0 Nwrap=10 Ltau=0 Sequential_scheme=\"submatrix\" N_delay=32 NSweep=3"
one() {  # label tile threads cpulist(or -) ; runs NBin 1 and 2
  local lab=$1 t=$2 th=$3 cpus=$4 pre=""
  [ "$cpus" != - ] && pre="numactl -C $cpus"
  ( export ZE_AFFINITY_MASK=$t OMP_NUM_THREADS=$th
    [ "$cpus" != - ] && export OMP_PROC_BIND=close OMP_PLACES=cores
    for nb in 1 2; do $pre python3 $S/run1.py $R/$lab/nb$nb $B/ALF-intel.out $ST $P NBin=$nb > /dev/null; done
    c1=$(grep "CPU Time" $R/$lab/nb1/info | awk '{print $NF}'); c2=$(grep "CPU Time" $R/$lab/nb2/info | awk '{print $NF}')
    echo "$lab threads=$th cpus=$cpus s/sweep $(python3 -c "print('%.2f' % ((float('$c2')-float('$c1'))/3))")" )
}
one alone_t1   0 1  -
one alone_t12p 0 12 0-11
one alone_t14  0 14 -
for t in 0 1 2 3 4 5 6 7; do one c8_t12p_$t $t 12 $((12*t))-$((12*t+11)) & done; wait
for t in 0 1 2 3 4 5 6 7; do one c8_t1_$t $t 1 - & done; wait
for t in 0 1 2 3 4 5 6 7; do one c8_t14_$t $t 14 - & done; wait
