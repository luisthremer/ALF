#!/bin/bash
# usage: rc.sh MACH LABEL SEDFILE [HAM] ; env SCHEME=standard|delayed|submatrix  VAR=src|srcT|srcC|srcP  TMO=seconds
MACH=$1; LABEL=$2; SEDF=$3; HAM=${4:-Hubbard}
W=~/alf-gpu/work/dlsub; VAR=${VAR:-src}
case $MACH in
 nv)    . ~/alf-gpu/env/nvhpc.sh; export OMP_TARGET_OFFLOAD=MANDATORY;;
 intel) . ~/alf-gpu/env/intel.sh 2>/dev/null; export ONEAPI_DEVICE_SELECTOR=level_zero:0 OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1 MKL_NUM_THREADS=1;;
 cpu)   . ~/alf-gpu/env/gnu.sh; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1;;
esac
BIN=$W/$VAR/Prog/ALF.out; export TMPDIR=$W/tmp
d=$W/r_${LABEL}; rm -rf "$d"; mkdir -p "$d"
cp $W/${PSRC:-src}/Scripts_and_Parameters_files/Start/{parameters,seeds} $d/
sed -i "s/^ham_name .*/ham_name = \"$HAM\"/" $d/parameters
[ -n "$SEDF" ] && [ -f "$SEDF" ] && sed -i -f "$SEDF" $d/parameters
[ -n "$SCHEME" ] && sed -i "s/^Sequential_scheme .*/Sequential_scheme = \"$SCHEME\"/" $d/parameters
cd $d; s=$(date +%s.%N); timeout ${TMO:-300} $BIN > out.log 2>&1; rc=$?; e=$(date +%s.%N)
act=$(grep -h "GPU path (gpu_state)" info 2>/dev/null | cut -c1-100)
echo "$LABEL rc=$rc wall=$(printf %.1f $(echo "$e - $s" | bc)) $act"
[ $rc -ne 0 ] && tail -3 out.log | cut -c1-140
exit 0
