#!/bin/bash
# Writes the case list of the gate (B = build directory): Mz and SU(2); L = 6, 8, 12, 16; Ltau 0/1; delayed/submatrix
# (N_delay 16) with the default build; the same with STABLOG and the -O1 build at Ltau = 1; L = 24 once per model.
B=$1; ST=$B/src_none/Scripts_and_Parameters_files/Start
c() { echo "$@"; }
for b in intel log o1; do for mz in T F; do for L in 6 8 12 16; do for lt in 0 1; do for sc in delayed submatrix; do
  [ $b != intel ] && [ $lt = 0 ] && continue
  c "${b}_mz${mz}_L${L}_lt${lt}_$sc $B/ALF-$b.out $ST L1=$L L2=$L Mz=.$mz. Ltau=$lt Sequential_scheme=\"$sc\" N_delay=16 NSweep=10 NBin=2"
done; done; done; done; done
for mz in T F; do c "intel_mz${mz}_L24_lt0_submatrix $B/ALF-intel.out $ST L1=24 L2=24 Mz=.$mz. Ltau=0 Sequential_scheme=\"submatrix\" N_delay=16 NSweep=5 NBin=2"; done
# negative controls (must FAIL): device update fault (acts on the standard update path) and device observable fault
c "negctl_faultupd_mzT_L8 $B/ALF-intel.out $ST env:ALF_GPU_DEBUG=1 env:ALF_GPU_FAULT=upd L1=8 L2=8 Mz=.T. Ltau=0 Sequential_scheme=\"standard\" NSweep=5 NBin=2"
c "negctl_faultobs_mzF_L8 $B/ALF-intel.out $ST env:ALF_GPU_DEBUG=1 env:ALF_GPU_FAULT=obs L1=8 L2=8 Mz=.F. Ltau=1 Sequential_scheme=\"submatrix\" N_delay=16 NSweep=5 NBin=2"
