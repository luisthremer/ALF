#!/bin/bash
# S-b: speed of this line (B) vs gpu-port (GP) at L = 24 / 32 and Precision Green at L = 24 / 28 / 32 on one node,
# 8 runs in parallel, one per tile (node fully busy, as in production). Parameters: the reference production set
# (jakob_params: Beta = L, N_FL = 2, Mz = F, ham_h0 = 1, U = 5.3, Nwrap 10, submatrix) with N_delay = 128 as in
# gpu-port's Intel numbers; h0 = 0 once more for the like-for-like case (GP measures on the device only with h0 = 0).
# s/sweep = (CPU time of NBin = 3 - CPU time of NBin = 1) / 10 (NSweep = 5), so start-up does not count.
#SBATCH -J p2i_sb
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --time=00:59:00
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH -o %x.%j.out
module load slurm_setup
S=$HOME/p2_intel/scripts; B=$HOME/p2_intel/build; R=$SCRATCH/p2_intel/sb_$SLURM_JOB_ID
. $S/env_smuc2.sh
export OMP_NUM_THREADS=14
mkdir -p $R
STB=$B/src_none/Scripts_and_Parameters_files/Start; STG=$B/gp/Scripts_and_Parameters_files/Start
par() { echo "L1=$1 L2=$1 Beta=$1.d0 N_Phi=$(( $1*$1/2 )) N_FL=2 Mz=.F. ham_h0=$2 ham_U=5.3d0 Nwrap=10 Ltau=0 Sequential_scheme=\"submatrix\" N_delay=128 NSweep=5"; }
task() {  # tile label exe start L h0 [env...]
  local t=$1 lab=$2 exe=$3 st=$4 L=$5 h=$6; shift 6
  ( export ZE_AFFINITY_MASK=$t "$@"
    for nb in 1 3; do python3 $S/run1.py $R/$lab/nb$nb $exe $st $(par $L $h) NBin=$nb; done ) > $R/$lab.txt 2>&1 &
}
task 0 B_L24_h1  $B/ALF-intel.out  $STB 24 1.d0
task 1 GP_L24_h1 $B/gp/Prog/ALF.out $STG 24 1.d0
task 2 B_L24_h0  $B/ALF-intel.out  $STB 24 0.d0
task 3 GP_L24_h0 $B/gp/Prog/ALF.out $STG 24 0.d0
task 4 B_L32_h1  $B/ALF-intel.out  $STB 32 1.d0
task 5 GP_L32_h1 $B/gp/Prog/ALF.out $STG 32 1.d0
task 6 B_L24_h1_LU $B/ALF-intel.out $STB 24 1.d0 ALF_GPU_DEBUG=1 ALF_GPU_CGR_LU=1
task 7 B_L28_h1  $B/ALF-intel.out  $STB 28 1.d0
wait
cd $R && for f in *.txt; do echo "== $f"; cat $f; done
for lab in B_L24_h1 GP_L24_h1 B_L24_h0 GP_L24_h0 B_L32_h1 GP_L32_h1 B_L24_h1_LU B_L28_h1; do
  c1=$(grep -h "CPU Time" $R/$lab/nb1/info | awk '{print $NF}'); c3=$(grep -h "CPU Time" $R/$lab/nb3/info | awk '{print $NF}')
  echo "$lab s/sweep $(python3 -c "print('%.2f' % ((float('$c3')-float('$c1'))/10))" 2>/dev/null)"
done
