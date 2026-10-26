#!/bin/bash
# usage: build_prod.sh gnu|intel  -> MPI builds, default and LOG stabilization, binaries in ~/alf-gpu/bin/
t=$1; M=$([ $t = gnu ] && echo GNU || echo IntelLLVM)
. ~/alf-gpu/env/$t.sh
cd ~/alf-gpu/src/ALF-$t
for stab in DEFAULT LOG; do
  arg=$([ $stab = LOG ] && echo LOG || echo "")
  { make clean >/dev/null 2>&1
    . ./configure.sh $M MPI $arg HDF5 NO-INTERACTIVE
    echo "ALF_FC=$ALF_FC"; echo "ALF_FLAGS_PROG=$ALF_FLAGS_PROG"; echo "ALF_LIB=$ALF_LIB"
    /usr/bin/time -f "make wall %e s" make -j8 all && cp Prog/ALF.out ~/alf-gpu/bin/ALF-$t-$stab.out && echo BUILD_OK
  } > ~/alf-gpu/logs/build_${t}_$stab.log 2>&1
  echo "$t $stab: $(grep -cE BUILD_OK ~/alf-gpu/logs/build_${t}_$stab.log) $(grep 'make wall' ~/alf-gpu/logs/build_${t}_$stab.log)"
done
