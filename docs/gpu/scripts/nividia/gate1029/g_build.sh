#!/bin/bash
# Gate 2026-10-29, nividia: build one variant from a read-only source tree that came from `git archive <sha>`.
#   usage: g_build.sh SHA12 NAME TOOL [configure args...]     (ALF_FLAGS_EXT may be set by the caller, it is recorded)
#   source: $G/src_<SHA12> (with SOURCE_SHA, TREE_SHA, TAR_SHA256); result: $G/bin/<SHA12>_<NAME>.out + .build.txt stamp
# Exit code: 0 only if configure and make succeed and the binary is new.
set -o pipefail
G=$HOME/alf-gpu/work/gate_1029; export TMPDIR=$G/tmp; mkdir -p $TMPDIR $G/bin $G/logs
S=$1; N=$2; TOOL=$3; shift 3
case "$S$N" in ""|*..*|*/*) echo "unsafe name"; exit 2;; esac
SRC=$G/src_$S; [ -f $SRC/SOURCE_SHA ] || { echo "no source $SRC"; exit 2; }
case $TOOL in gnu) . ~/alf-gpu/env/gnu.sh;; nvhpc) . ~/alf-gpu/env/nvhpc.sh;; *) echo "bad tool"; exit 2;; esac
export ALF_HDF5_DIR=$HOME/alf-gpu/work/t207/h5        # read only (prebuilt HDF5 1.14.6 per compiler)
T=$G/tree_${S}_$N; B=$G/bin/${S}_$N.out
rm -rf -- "$T" "$B"; mkdir -p "$T"; cp -a $SRC/. "$T/"; chmod -R u+w "$T"; cd "$T" || exit 2
. ./configure.sh "$@" NO-INTERACTIVE > $G/logs/conf_${S}_$N.log 2>&1 || { echo "configure failed"; exit 3; }
make -j4 program > $G/logs/make_${S}_$N.log 2>&1 || { echo "make failed"; exit 4; }
cp Prog/ALF.out "$B" || exit 5
{ echo "source_sha $(cat $SRC/SOURCE_SHA)"; echo "tree_sha $(cat $SRC/TREE_SHA)"; echo "tar_sha256 $(cat $SRC/TAR_SHA256)"
  echo "binary_sha256 $(sha256sum $B | cut -d' ' -f1)"; echo "configure $* ALF_FLAGS_EXT=${ALF_FLAGS_EXT:-}"
  echo "flags_prog $ALF_FLAGS_PROG"; echo "lib $ALF_LIB"
  echo "nvfortran $(nvfortran --version 2>/dev/null | sed -n 2p)"; echo "gfortran $(gfortran --version 2>/dev/null | head -1)"
  echo "built $(date -Is)"; } > $B.build.txt
echo "built ${S}_$N $(sha256sum $B | cut -c1-16)"
