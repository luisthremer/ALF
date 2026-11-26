#!/bin/bash
# Login-node check of the SuperMUC-NG Phase 2 machine block with MPI and HDF5 (production-style builds):
#   source configure.sh SUPERMUCNG2 MPI GPU=intel HDF5 NO-INTERACTIVE && make   (and the same with GPU=none)
# usage: build_mpi_check.sh SRC_DIR OUT_DIR
src=$1; out=$2
mkdir -p "$out"
for t in intel none; do
  rm -rf "$out/mpi_$t"; cp -r "$src" "$out/mpi_$t"
  ( cd "$out/mpi_$t"
    source configure.sh SUPERMUCNG2 MPI GPU=$t HDF5 NO-INTERACTIVE > "$out/configure_mpi_$t.log" 2>&1
    make -j16 > "$out/make_mpi_$t.log" 2>&1
    echo "mpi $t: make rc=$? $(ls -la Prog/ALF.out 2>&1)" )
done
