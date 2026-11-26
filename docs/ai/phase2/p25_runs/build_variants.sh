#!/bin/bash
# Builds the GPU=host variants needed for the P2.5 evidence (one detached worktree per variant, -j16).
# usage: build_variants.sh <base repo (a checkout of the commit to test)> <dir for the variant trees>
BASE=$(cd "$1" && pwd); OUT=$2; mkdir -p "$OUT"; REV=$(git -C "$BASE" rev-parse HEAD)
build() {   # name, configure args...
  name=$1; shift
  d=$OUT/$name
  [ -d "$d" ] || git -C "$BASE" worktree add --detach "$d" "$REV" >/dev/null 2>&1
  ( cd "$d" && { [ "$name" = qrref ] && sed -i 's/^# STABCONFIGURATION="\${STABCONFIGURATION} -DQRREF"/STABCONFIGURATION="${STABCONFIGURATION} -DQRREF"/' configure.sh; true; } && source ./configure.sh "$@" GPU=host > "$OUT/$name.cfg.log" 2>&1 \
      && make -j16 > "$OUT/$name.make.log" 2>&1 && echo "$name: built $(ls -la Prog/ALF.out | awk '{print $5}')" || echo "$name: BUILD FAILED" )
}
build log     GNU noMPI LOG
build stab1   GNU noMPI STAB1
build stab2   GNU noMPI STAB2
build stab3   GNU noMPI STAB3
build mpi     GNU
build tempering GNU Tempering
build pp      GNU PP
build qrref   GNU noMPI
build hdf5    GNU noMPI HDF5 NO-INTERACTIVE
