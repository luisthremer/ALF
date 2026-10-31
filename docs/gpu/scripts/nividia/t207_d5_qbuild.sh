#!/bin/bash
# Ticket 207 dev5 (accuracy test against a quad reference, development only), nividia: build the offline quad-reference program
# t207_d5_qref.F90 (gfortran, real(16)) against the CPU-oracle objects of a GNU build tree of t207_d5_build.sh.
#   usage: t207_d5_qbuild.sh TREE OUT CPPFLAGS     e.g. t207_d5_qbuild.sh gnu_pp qref_pp -DQR_PREPIVOT ; ... gnu_log qref_log -DSTABLOG
# Parts: 1 module q5 (quad reference) and 3 the program with -O2 (no fast-math); 2 module q5nul (null copy of CGR2_2) with the tree's
# ALF_FLAGS_PROG, as cgr2_2_mod.  Result: dev5/bin/<OUT>.out (+ .sha256); objects in dev5/qbuild_<OUT> (build tree, removable).
W=$HOME/alf-gpu/work/t207/dev5; export TMPDIR=$W/tmp; . ~/alf-gpu/env/gnu.sh > /dev/null
TR=$1; OUT=$2; FL=$3
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
safe "$TR"; safe "$OUT"
[ -e "$W/bin/$OUT.out" ] && { echo "refusing: $OUT exists (binaries are kept)"; exit 1; }
T=$W/tree_$TR; mkdir -p $W/qbuild_$OUT; cd $W/qbuild_$OUT || exit 1
LIBS=$(grep '^ALF_LIB=' $W/logs/env_$TR.txt | cut -d= -f2-)
PF=$(grep '^ALF_FLAGS_PROG=' $W/logs/env_$TR.txt | cut -d= -f2-)
{ gfortran -O2 -fopenmp -cpp -DQ5_PART=1 $FL -c $W/t207_d5_qref.F90 -o q5.o &&
  gfortran $PF -I$T/Prog -I. -DQ5_PART=2 -c $W/t207_d5_qref.F90 -o q5nul.o &&
  gfortran -O2 -fopenmp -cpp -DQ5_PART=3 $FL -I$T/Prog -I$T/Libraries/Modules -I. -c $W/t207_d5_qref.F90 -o q5main.o &&
  gfortran -fopenmp q5main.o q5nul.o q5.o $T/Prog/cgr2_2_mod.o $T/Prog/udv_state_mod.o $T/Prog/QDRP_decompose_mod.o $LIBS -o $W/bin/$OUT.out
} > $W/logs/qbuild_$OUT.log 2>&1 || { echo "qbuild failed"; tail -30 $W/logs/qbuild_$OUT.log; exit 2; }
sha256sum $W/t207_d5_qref.F90 > $W/bin/$OUT.src.sha256
sha256sum $W/bin/$OUT.out > $W/bin/$OUT.out.sha256 && echo "built $OUT $(cut -c1-16 $W/bin/$OUT.out.sha256) (part 2 flags: $PF; libs: $LIBS)"
