#!/bin/bash
# usage: run_all.sh PREPDIR RUNDIR JOBS case:build [case:build ...]
# build in {def-up,def-pp,log-up,log-pp,def-orig}; binary ~/alf-gpu/work/t104/bin/ALF-<build>.out
# Copies PREPDIR/<case>/{parameters,seeds} to RUNDIR/<case>/<build>/ and runs 1 rank, JOBS in parallel.
. ~/alf-gpu/env/gnu.sh
prep=$(readlink -f $1); run=$(readlink -f $2); jobs=$3; shift 3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
one() {
  c=${1%%:*}; b=${1##*:}; d=$RUN/$c/$b
  [ -f $d/timing.txt ] && return
  mkdir -p $d && cp $PREP/$c/parameters $PREP/$c/seeds $d/ && cd $d || return
  s=$(date +%s.%N)
  mpirun -n 1 $HOME/alf-gpu/work/t104/bin/ALF-$b.out > run.log 2>&1
  e=$(date +%s.%N)
  printf 'wall_s=%.1f\n' "$(echo "$e - $s" | bc)" > timing.txt
  rm -f confout_*
}
export -f one; export PREP=$prep RUN=$run
printf '%s\n' "$@" | xargs -P $jobs -I{} bash -c 'one {}'
