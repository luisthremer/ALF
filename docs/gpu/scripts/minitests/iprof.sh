#!/bin/bash
# usage: iprof.sh LABEL L NSWEEP NBIN scheme [Ltau] ; runs ALF on the B580
LABEL=$1; L=$2; NS=$3; NB=$4; sch=$5; LT=${6:-1}
. ~/alf-gpu/env/intel.sh 2>/dev/null
export ONEAPI_DEVICE_SELECTOR=level_zero:0 OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 TMPDIR=~/alf-gpu/work/dlsub/tmp
SRC=${SRC:-~/alf-gpu/work/dlsub/src}
d=~/alf-gpu/work/dlsub/q_${LABEL}_$sch; rm -rf $d; mkdir -p $d
cp $SRC/Scripts_and_Parameters_files/Start/{parameters,seeds} $d/
sed -i "s/^Sequential_scheme .*/Sequential_scheme   = \"$sch\"/; s/^NSweep .*/NSweep = $NS/; s/^NBin .*/NBin = $NB/; s/^L1 .*/L1 = $L/; s/^L2 .*/L2 = $L/; s/^Ltau .*/Ltau = $LT/" $d/parameters
cd $d; s=$(date +%s.%N); $SRC/Prog/ALF.out > out.log 2>&1; rc=$?; e=$(date +%s.%N)
echo "$sch L=$L rc=$rc wall=$(echo "$e - $s" | bc) $(grep -h 'GPU time in seq' info 2>/dev/null | cut -c40-70) | $(grep -h 'GPU path (gpu_state)' info | cut -c1-90)"
