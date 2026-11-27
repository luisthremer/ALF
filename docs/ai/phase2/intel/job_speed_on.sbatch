#!/bin/bash
# Speed with the Intel switches: this line off / on (ASYNCBLK=1 BLKFUSE=5 BLKPF=1 OBSCORR=1), intel-batched tip with the same
# switches, gpu-port tip; L = 24 and 32, reference set (Beta = L, N_FL 2, U 5.3, h0 1, Nwrap 10, submatrix N_delay 32), one run per
# tile, 8 at a time (node full). s/sweep = (CPU time NBin 2 - NBin 1) / 3 (NSweep 3).
#SBATCH -J p2i_spd
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --time=00:59:00
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH -o %x.%j.out
module load slurm_setup
S=$HOME/p2_intel/scripts; B=$HOME/p2_intel/build; R=$SCRATCH/p2_intel/spd_$SLURM_JOB_ID
. $S/env_smuc2.sh
export OMP_NUM_THREADS=14
mkdir -p $R
STB=$B/src_none/Scripts_and_Parameters_files/Start; STG=$B/gp/Scripts_and_Parameters_files/Start; STI=$B/ib/Scripts_and_Parameters_files/Start
ON="ALFSH_ASYNCBLK=1 ALFSH_BLKFUSE=5 ALFSH_BLKPF=1 ALFSH_OBSCORR=1"
par() { echo "L1=$1 L2=$1 Beta=$1.d0 N_Phi=$(( $1*$1/2 )) N_FL=2 Mz=.F. ham_h0=1.d0 ham_U=5.3d0 Nwrap=10 Ltau=0 Sequential_scheme=\"submatrix\" N_delay=32 NSweep=3"; }
task() {  # tile label exe start L [env...]
  local t=$1 lab=$2 exe=$3 st=$4 L=$5; shift 5
  ( export ZE_AFFINITY_MASK=$t "$@"
    for nb in 1 2; do python3 $S/run1.py $R/$lab/nb$nb $exe $st $(par $L) NBin=$nb > /dev/null; done
    c1=$(grep "CPU Time" $R/$lab/nb1/info | awk '{print $NF}'); c2=$(grep "CPU Time" $R/$lab/nb2/info | awk '{print $NF}')
    echo "$lab s/sweep $(python3 -c "print('%.2f' % ((float('$c2')-float('$c1'))/3))") | $(grep -h 'Intel block\|OBSCORR\|prefetch' $R/$lab/nb2/info | tr -s ' ' | tr '\n' ';')" ) > $R/$lab.txt 2>&1 &
}
task 0 B_off_L24  $B/ALF-intel.out   $STB 24
task 1 B_on_L24   $B/ALF-intel.out   $STB 24 $ON
task 2 IB_on_L24  $B/ib/Prog/ALF.out $STI 24 $ON
task 3 GP_L24     $B/gp/Prog/ALF.out $STG 24
task 4 B_off_L32  $B/ALF-intel.out   $STB 32
task 5 B_on_L32   $B/ALF-intel.out   $STB 32 $ON
task 6 IB_on_L32  $B/ib/Prog/ALF.out $STI 32 $ON
task 7 GP_L32     $B/gp/Prog/ALF.out $STG 32
wait
cat $R/*.txt
