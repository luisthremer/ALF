#!/bin/bash
# usage: run108.sh PREPDIR RUNDIR JOBS case:build [case:build ...]
# build = name of ~/alf-gpu/work/t108/bin/ALF-<build>.out. Same parameters+seeds for every build of a case.
export TMPDIR=$HOME/alf-gpu/work/t108/tmp
. ~/alf-gpu/env/gnu.sh
prep=$(readlink -f $1); run=$(readlink -f $2); jobs=$3; shift 3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
one() {
  c=${1%%:*}; b=${1##*:}; d=$RUN/$c/$b
  [ -f $d/timing.txt ] && return
  mkdir -p $d && cp $PREP/$c/parameters $PREP/$c/seeds $d/ && cd $d || return
  s=$(date +%s.%N)
  timeout 900 mpirun -n 1 $HOME/alf-gpu/work/t108/bin/ALF-$b.out > run.log 2>&1
  e=$(date +%s.%N)
  printf 'wall_s=%.1f\n' "$(echo "$e - $s" | bc)" > timing.txt
  rm -f confout_*
}
export -f one; export PREP=$prep RUN=$run
printf '%s\n' "$@" | xargs -P $jobs -I{} bash -c 'one {}'
