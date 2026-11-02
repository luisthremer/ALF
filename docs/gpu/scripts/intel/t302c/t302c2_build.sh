#!/bin/bash
# Ticket 302c-2 (development only), on intel: build one variant in a fresh tree from src/ (as t207_d5_build.sh on nividia).
#   usage: [ALF_FLAGS_EXT=-DQR_PREPIVOT] [J=6] t302c2_build.sh NAME intel|gnu [configure args...]
#     intel: oneAPI of ~/alf-gpu/env/intel.sh (ifx 2026.2, oneMKL 2026.2); gnu: system gfortran 15.2 + OpenBLAS (-llapack -lblas)
# Result: bin/ALF_<NAME>.out (+ .sha256, .srctree = hash of the Prog/Libraries sources, .srcstate = HEAD/diff of the laptop worktree);
#         tree_<NAME> (build tree, removable; its objects are the CPU oracle of the offline program for gnu_*); logs/{conf,make,env}_<NAME>.
W=${T302C2_W:-$HOME/alf-gpu/work/t302c-2}   # round 2: T302C2_W=$HOME/alf-gpu/work/t302c-2b, round 3: .../t302c-2c (only these three)
case "$W" in "$HOME/alf-gpu/work/t302c-2"|"$HOME/alf-gpu/work/t302c-2b"|"$HOME/alf-gpu/work/t302c-2c") ;; *) echo "refusing work dir $W"; exit 1;; esac
:; export TMPDIR=$W/tmp; mkdir -p $TMPDIR $W/logs $W/bin
NAME=$1; TOOL=$2; shift 2
case "$NAME" in ""|*..*|*/*) echo "refusing unsafe name: [$NAME]"; exit 1;; esac; : "${HOME:?}"
[ -e "$W/bin/ALF_$NAME.out" ] && { echo "refusing: binary $NAME exists (binaries are kept)"; exit 1; }
case $TOOL in
  intel) . ~/alf-gpu/env/intel.sh 2>/dev/null
         [ -d $W/h5 ] || { mkdir -p $W/h5 && cp -a ~/alf-gpu/hdf5/. $W/h5/; }
         export ALF_HDF5_DIR=$W/h5 ;;
  gnu) ;;
  *) echo "bad tool $TOOL"; exit 1 ;;
esac
rm -rf -- "$W/tree_$NAME"; mkdir -p "$W/tree_$NAME"; rsync -a "$W/src/" "$W/tree_$NAME/"
cd "$W/tree_$NAME" || exit 1
. ./configure.sh "$@" NO-INTERACTIVE > $W/logs/conf_$NAME.log 2>&1 || { echo "configure failed"; tail -20 $W/logs/conf_$NAME.log; exit 2; }
env | grep '^ALF_' | sort > $W/logs/env_$NAME.txt
/usr/bin/time -f "make wall %e s" make -j${J:-6} program > $W/logs/make_$NAME.log 2>&1 || { echo "make failed"; grep -i -B2 -A8 "error" $W/logs/make_$NAME.log | head -80; exit 3; }
tail -1 $W/logs/make_$NAME.log
(cd "$W/src" && find Prog Libraries -type f \( -name "*.F90" -o -name "*.f90" -o -name Makefile \) | LC_ALL=C sort | xargs sha256sum | sha256sum | cut -c1-16) > $W/bin/ALF_$NAME.srctree
cp "$W/src/.srcstate" $W/bin/ALF_$NAME.srcstate
{ $ALF_FC --version 2>&1 | head -1; echo "ALF_FLAGS_PROG=$ALF_FLAGS_PROG"; echo "ALF_LIB=$ALF_LIB"; } > $W/bin/ALF_$NAME.buildinfo
cp Prog/ALF.out $W/bin/ALF_$NAME.out && sha256sum $W/bin/ALF_$NAME.out > $W/bin/ALF_$NAME.out.sha256 && echo "built $NAME $(cut -c1-16 $W/bin/ALF_$NAME.out.sha256)"
