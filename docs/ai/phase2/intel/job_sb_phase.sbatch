#!/bin/bash
# S-b: determinant phase at N >~ 2950 (KNOWN_ISSUES "GPU: NaN determinant phase"): L = 56 / 60 / 64, Beta = 6.4
# (Nt = 64), U = 5, h0 = 0, submatrix N_delay 32, 1 sweep x 2 bins, Ltau = 0 (as gpu-port's check_large_N_phase.sh);
# tracked det(U) (default) and the LU path (ALF_GPU_CGR_LU=1, LU_PHASE_SAFE rule). Pass: acceptance >= 0.30, no NaN.
# Symm = .F.: with Symm = .T. the CPU start-up check Op_V_is_symmetric takes hours at this N (KNOWN_ISSUES P1); one bin.
#SBATCH -J p2i_ph
#SBATCH --account=pn93la
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --time=01:30:00
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH -o %x.%j.out
module load slurm_setup
S=$HOME/p2_intel/scripts; B=$HOME/p2_intel/build; R=$SCRATCH/p2_intel/ph_$SLURM_JOB_ID
. $S/env_smuc2.sh
export OMP_NUM_THREADS=14
mkdir -p $R
ST=$B/src_none/Scripts_and_Parameters_files/Start
t=0
for L in 56 60 64; do for mode in tracked lu; do
  ( export ZE_AFFINITY_MASK=$t; [ $mode = lu ] && export ALF_GPU_DEBUG=1 ALF_GPU_CGR_LU=1
    python3 $S/run1.py $R/L${L}_$mode $B/ALF-intel.out $ST L1=$L L2=$L Beta=6.4d0 N_FL=2 Mz=.F. ham_h0=0.d0 ham_U=5.0d0 \
      Symm=.F. Nwrap=10 Ltau=0 Sequential_scheme=\"submatrix\" N_delay=32 NSweep=1 NBin=1 ) > $R/L${L}_$mode.txt 2>&1 &
  t=$((t+1))
done; done
wait
for f in $R/*.txt; do cat $f; acc=$(grep -m1 "Acceptance  " ${f%.txt}/info | awk '{print $NF}')
  awk -v a="$acc" -v f="$(basename $f .txt)" 'BEGIN{ ok = !(a == "" || a+0 != a+0 || a+0 < 0.30); print f, "acceptance", a, (ok ? "OK" : "FAIL") }'; done
