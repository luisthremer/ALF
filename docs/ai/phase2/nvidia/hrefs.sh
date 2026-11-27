#!/bin/bash -l
# GPU=none (nvfortran) and MPI GPU=nvidia builds of a candidate, then submit ctest, device suites, MPI
T=$1; P=~/p2_nvidia
for v in none mpi; do
  rm -rf $P/cand_${T}_$v; mkdir $P/cand_${T}_$v
  (cd $P/cand_$T && tar cf - --exclude=testsuite/tests --exclude='*.o' --exclude='*.mod' --exclude='*.a' --exclude=ALF.out .) | tar xf - -C $P/cand_${T}_$v
  cd $P/cand_${T}_$v
  if [ $v = none ]; then a="noMPI GPU=none"; else a="MPI GPU=nvidia"; fi
  (. ./configure.sh HELMAGPU $a NO-INTERACTIVE > cfg.log 2>&1; make clean >/dev/null 2>&1; make -j16 lib program > build.log 2>&1) && echo "$v ok" || { echo "$v FAIL"; tail -3 build.log; }
done
cd $P
sbatch ctest.sbatch $T
sbatch dev.sbatch $T traj phase g2 speed speedgp trajlog
sbatch mpi.sbatch $T
