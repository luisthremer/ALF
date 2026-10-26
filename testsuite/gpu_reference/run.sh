#!/bin/bash
# usage: run.sh ALF.out NRANKS RUNDIR [RUNDIR ...]
# Runs each prepared directory with NRANKS MPI ranks (1 OpenMP/BLAS thread each), records wall time.
set -euo pipefail
bin=$(readlink -f "$1"); n=$2; shift 2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for d in "$@"; do
  ( cd "$d"
    start=$(date +%s.%N)
    mpirun -n "$n" "$bin" > run.log 2>&1
    end=$(date +%s.%N)
    printf 'binary=%s\nranks=%s\nwall_s=%.1f\n' "$bin" "$n" "$(echo "$end - $start" | bc)" > timing.txt
    echo "$d: $(grep wall_s timing.txt)" )
done
