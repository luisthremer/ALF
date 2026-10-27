#!/bin/sh
# Build and run both backends; write gpu_bench/results/nividia_A30_<date>.txt
# usage (on nividia, from gpu_bench/):  . ~/alf-gpu/env/nvhpc.sh; ./run.sh [N-list, e.g. 64,144,256]
set -e
cd "$(dirname "$0")"
NLIST=${1:-}
DATE=$(date +%Y-%m-%d)
OUT=results/nividia_A30_${DATE}.txt
mkdir -p results
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=false
{
  echo "# gpu_bench run $DATE on $(hostname)"
  echo "# compiler: $(nvfortran --version | sed -n 2p)"
  nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
  echo "# N list: ${NLIST:-default}"
  echo "# host load at start: $(cut -d' ' -f1-3 /proc/loadavg)"
} > $OUT
for B in cuda host; do
  make BACKEND=$B > build_$B.log 2>&1 || { cat build_$B.log; exit 1; }
  echo "" >> $OUT
  echo "################ BACKEND=$B ################" >> $OUT
  if [ $B = cuda ]; then export OMP_TARGET_OFFLOAD=MANDATORY; else unset OMP_TARGET_OFFLOAD; fi
  # exit status 1 = at least one check FAILED (printed in the file); keep going
  # stdout -> table file; stderr (nvfortran runtime 'ieee_inexact' warning, time) kept apart so tables stay intact
  /usr/bin/time -o build/time_$B.txt -f "# wall %e s (BACKEND=$B)" ./build/$B/bench $NLIST >> $OUT 2> build/stderr_$B.txt || echo "# bench exit status $? (at least one check FAILED, see FAIL lines)" >> $OUT
  cat build/time_$B.txt >> $OUT
  sed 's/^/# stderr: /' build/stderr_$B.txt >> $OUT
done
rm -f build_*.log
echo "results: $OUT"
