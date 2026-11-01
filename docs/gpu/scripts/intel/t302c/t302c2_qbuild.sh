#!/bin/bash
# Ticket 302c-2 (development only), on intel: build the offline binary128 program t207_d5_qref.F90 (207-2 pipeline, commit 5d876291)
# against the CPU-oracle objects of a GNU build tree of t302c2_build.sh; the same three parts and flags as t207_d5_qbuild.sh on nividia
# (only the work dir and the compiler environment differ: system gfortran 15.2 + OpenBLAS on intel).
#   usage: t302c2_qbuild.sh TREE OUT CPPFLAGS     t302c2_qbuild.sh gnu_pp qref_pp -DQR_PREPIVOT ; t302c2_qbuild.sh gnu_log qref_log -DSTABLOG
# Result: bin/<OUT>.out (+ .sha256, .src.sha256 = sha256 of the t207_d5_qref.F90 next to t207_d5_proc.sh, which proc.sh verifies).
W=$HOME/alf-gpu/work/t302c-2; export TMPDIR=$W/tmp
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
