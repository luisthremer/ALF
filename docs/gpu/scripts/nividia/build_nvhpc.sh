#!/bin/bash
# usage: build_nvhpc.sh [SRC_DIR]   (default ~/alf-gpu/work/t105)
# nvfortran 26.9 CPU builds with the unchanged upstream PGI profile of configure.sh (pgfortran = nvfortran).
#  1. devel noMPI unit-test build + ctest      -> logs/t105_ctest.log
#  2. production MPI HDF5 (HPC-X mpifort)      -> ~/alf-gpu/bin/ALF-nvhpc-DEFAULT.out
#  3. production noMPI HDF5                    -> ~/alf-gpu/bin/ALF-nvhpc-noMPI.out
# NB: do not pipe `. configure.sh` into tail/tee, the exports would be lost in the subshell.
set -o pipefail
SRC=${1:-$HOME/alf-gpu/work/t105}
. ~/alf-gpu/env/nvhpc.sh
cd "$SRC" || exit 1
L=~/alf-gpu/logs; mkdir -p $L ~/alf-gpu/bin

make clean >/dev/null 2>&1
. ./configure.sh PGI devel noMPI NO-INTERACTIVE > $L/t105_conf_ctest.log 2>&1
make -j4 all > $L/t105_make_ctest.log 2>&1 || { echo "ctest build failed"; exit 2; }
( cd testsuite && rm -rf tests && mkdir tests && cd tests \
  && cmake .. >/dev/null && cmake --build . -j4 > $L/t105_cmake.log 2>&1 \
  && ctest -j4 --output-on-failure > $L/t105_ctest.log 2>&1; tail -4 $L/t105_ctest.log )

make clean >/dev/null 2>&1
. ./configure.sh PGI MPI HDF5 NO-INTERACTIVE > $L/t105_conf_mpi.log 2>&1
make -j4 all > $L/t105_make_mpi.log 2>&1 && cp Prog/ALF.out ~/alf-gpu/bin/ALF-nvhpc-DEFAULT.out && echo "MPI HDF5 OK"

make clean >/dev/null 2>&1
. ./configure.sh PGI noMPI HDF5 NO-INTERACTIVE > $L/t105_conf_nompi.log 2>&1
make -j4 all > $L/t105_make_nompi.log 2>&1 && cp Prog/ALF.out ~/alf-gpu/bin/ALF-nvhpc-noMPI.out && echo "noMPI HDF5 OK"
