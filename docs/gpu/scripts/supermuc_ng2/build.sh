#!/bin/bash
# usage: [AOT=1] bash -l build_generic.sh NAME configure-args...   Builds from ~/alf-gpu/srcm into ~/alf-gpu/bld_NAME, installs ~/alf-gpu/bin/ALF-NAME.out
# AOT=1: ahead-of-time compilation for PVC through the ifx wrapper (~/alf-gpu/wrap/ifx). HDF5 in args: LRZ HDF5 module via a shim dir.
NAME=$1; shift
source ~/alf-gpu/env251.sh >/dev/null 2>&1 || exit 1
if [ -n "$AOT" ]; then export REAL_IFX=$(which ifx); export PATH=$HOME/alf-gpu/wrap:$PATH; echo "AOT wrapper: $(which ifx)"; fi
B=~/alf-gpu/bld_$NAME; rm -rf $B; mkdir -p $B ~/alf-gpu/bin
rsync -a --exclude='*.o' --exclude='*.mod' --exclude='*.a' --exclude='ALF.out' ~/alf-gpu/srcm/ $B/
cd $B; export ALF_DIR=$PWD
if echo "$@" | grep -q HDF5; then
  ifx -o /tmp/gcv.$$ get_compiler_version.F90 && CV=$(/tmp/gcv.$$ | sed 's/[ ,()]/_/g')
  export ALF_HDF5_DIR=$HOME/alf-gpu/h5shim_$NAME; S=$ALF_HDF5_DIR/$CV
  rm -rf $ALF_HDF5_DIR; mkdir -p $S/bin $S/lib; ln -sfn $HDF5_BASE/include $S/include
  for f in $HDF5_BASE/lib/*; do ln -sfn $f $S/lib/$(basename $f); done
  ln -sfn $HDF5_BASE/lib/libhdf5_hl_fortran.a $S/lib/libhdf5hl_fortran.a
  printf '#!/bin/bash\n%s/bin/h5fc "$@" | sed -E "s/(Extra libraries:).*/\\1 -lhdf5_f90cstub -lhdf5_hl_f90cstub -lhdf5_hl -lhdf5 -lz -lm -ldl/"\n' "$HDF5_BASE" > $S/bin/h5fc; chmod +x $S/bin/h5fc
fi
ifx --version | head -1
source ./configure.sh "$@" || exit 2
make -j8 lib program 2>&1 | tail -6
cp Prog/ALF.out ~/alf-gpu/bin/ALF-$NAME.out && ls -la ~/alf-gpu/bin/ALF-$NAME.out
