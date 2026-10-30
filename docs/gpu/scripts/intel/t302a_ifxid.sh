#!/bin/bash
# acceptance 4, ifx part (development version): ifx 2025.3.3 CPU build of base and candidate (IntelLLVM noMPI HDF5), one run of
# hub6x6_cb with the same seeds, data.h5 / info compared. usage on nividia: t302a_ifxid.sh   (<= 2 cores)
W=$HOME/alf-gpu/work/t302a/dev; export TMPDIR=$W/tmp; mkdir -p $TMPDIR
export ALF_HDF5_DIR=$HOME/alf-gpu/hdf5
for v in base cand; do
  src=$W/base_src; [ $v = cand ] && src=$W/src
  rm -rf $W/ifx_$v; mkdir $W/ifx_$v; cp -a $src/. $W/ifx_$v/; cd $W/ifx_$v
  ( . $HOME/alf-gpu/env/intel.sh >/dev/null 2>&1; . ./configure.sh IntelLLVM noMPI HDF5 NO-INTERACTIVE > configure.log 2>&1 && make -j2 all > make.log 2>&1 ) \
    || { echo "build $v failed"; tail -5 make.log; exit 1; }
  ls Prog/ALF.out > /dev/null || exit 1
done
echo BUILDS_OK
