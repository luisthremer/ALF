#!/bin/bash
# Ticket 108: builds MPI+HDF5 STABLOG binaries into ~/alf-gpu/work/t108/bin/ALF-<name>.out
# up = modified sources, macro off; orig = unmodified (pre-104) cgr sources; v1/v2/v3nN = -DQR_PREPIVOT=1/2/3 [-DQR_PPNB=N]
export TMPDIR=$HOME/alf-gpu/work/t108/tmp
. ~/alf-gpu/env/gnu.sh
W=$HOME/alf-gpu/work/t108
cd $W/src
build() { # name ext
  make clean >/dev/null 2>&1
  if [ -n "$2" ]; then export ALF_FLAGS_EXT="$2"; else unset ALF_FLAGS_EXT; fi
  { . ./configure.sh GNU MPI LOG HDF5 NO-INTERACTIVE; echo "FLAGS=$ALF_FLAGS_PROG"
    make -j6 all && cp Prog/ALF.out ../bin/ALF-$1.out && echo BUILD_OK; } > ../build_$1.log 2>&1
  echo "$1 $(grep -c BUILD_OK ../build_$1.log)"
}
for spec in "$@"; do
  case $spec in
    up) build up "" ;;
    v1) build v1 "-DQR_PREPIVOT=1" ;;
    v2) build v2 "-DQR_PREPIVOT=2" ;;
    v3n*) n=${spec#v3n}; build $spec "-DQR_PREPIVOT=3 -DQR_PPNB=$n" ;;
    orig) cp Prog/cgr1_mod.F90 ../mod_cgr1.F90; cp Prog/cgr2_2_mod.F90 ../mod_cgr2_2.F90
          cp ../orig/*.F90 Prog/; build orig ""
          cp ../mod_cgr1.F90 Prog/cgr1_mod.F90; cp ../mod_cgr2_2.F90 Prog/cgr2_2_mod.F90 ;;
  esac
done
