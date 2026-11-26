#!/bin/bash
# Runs every P2.5 evidence case.  usage: run_all.sh <repo checkout with Prog/ALF.out (GPU=host, noMPI)> <variants dir (build_variants.sh)> <work dir>
REPO=$(cd "$1" && pwd); V=$(cd "$2" && pwd); W=$3; R=$REPO/docs/ai/phase2/p25_runs; S=$REPO/Scripts_and_Parameters_files/Start
P="python3 $R/p25_run.py --start $S --work $W"
rm -rf "$W"; mkdir -p "$W"
$P --build serial --exe $REPO/Prog/ALF.out $R/cases_main.py
$P --build serial --exe $REPO/Prog/ALF.out --restart $R/cases_restart.py
# stabilisation variants: LOG runs the whole table, STAB1/2/3 are refused
$P --build LOG --exe $V/log/Prog/ALF.out --work $W/log $R/cases_main.py
for s in stab1 stab2 stab3; do $P --build $s --exe $V/$s/Prog/ALF.out --work $W/$s $R/cases_stab.py; done
$P --build QRREF --exe $V/qrref/Prog/ALF.out --work $W/qrref --only hub-mz-4x4,hub-su2-4x4,sch-delayed-N4 $R/cases_main.py
# MPI: two independent chains (mpirun -np 2), then PARALLEL_PARAMS and TEMPERING with two groups of one rank
$P --build MPI-np2 --exe $V/mpi/Prog/ALF.out --work $W/mpi --mpi 2 --only hub-mz-4x4,hub-su2-4x4,sch-delayed-N4,sch-submatrix-N4,mv-global,mod-kondo $R/cases_main.py
$P --build PP-np2 --exe $V/pp/Prog/ALF.out --work $W/pp --mpi 2 --groups 2 $R/cases_group.py
$P --build TEMPERING-np2 --exe $V/tempering/Prog/ALF.out --work $W/tempering --mpi 2 --groups 2 $R/cases_tempering.py
# HDF5 (host-side output; variant built with `configure.sh GNU noMPI HDF5`), skipped when the variant was not built
if [ -x $V/hdf5/Prog/ALF.out ]; then
  export PATH=$(ls -d $V/hdf5/HDF5/*/bin | head -1):$PATH
  $P --build HDF5 --exe $V/hdf5/Prog/ALF.out --work $W/hdf5 --only hub-mz-4x4,hub-su2-4x4,sch-delayed-N4,mod-kondo $R/cases_main.py
fi
