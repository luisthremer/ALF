#!/bin/bash
# usage: runall.sh LABEL NSWEEP L [schemes...]  ; runs ALF.out per scheme in ~/alf-gpu/work/dlsub/run_LABEL_scheme
LABEL=$1; NS=$2; L=$3; shift 3
. ~/alf-gpu/env/nvhpc.sh
export TMPDIR=~/alf-gpu/work/dlsub/tmp OMP_TARGET_OFFLOAD=MANDATORY
SRC=~/alf-gpu/work/dlsub/src
for sch in "$@"; do
  d=~/alf-gpu/work/dlsub/run_${LABEL}_$sch; rm -rf $d; mkdir -p $d
  cp $SRC/Scripts_and_Parameters_files/Start/{parameters,seeds} $d/
  sed -i "s/^Sequential_scheme .*/Sequential_scheme   = \"$sch\"/; s/^NSweep .*/NSweep = $NS/; s/^NBin .*/NBin = 10/; s/^L1 .*/L1 = $L/; s/^L2 .*/L2 = $L/" $d/parameters
  ( cd $d; s=$(date +%s.%N); $SRC/Prog/ALF.out > out.log 2>&1; rc=$?; e=$(date +%s.%N); echo "$sch rc=$rc wall=$(echo "$e - $s" | bc)" )
done
