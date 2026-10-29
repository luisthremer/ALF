#!/bin/bash
# Acceptance harness, MACHINE side: GPU=none guarantee of a shipped candidate against its base (both from git archive).
#   usage: acc_gpunone.sh TICKET SHA12 OUTDIR
#   (1) preprocessed sources: every Fortran file that differs between base_src and src (computed, not typed), for 6
#       macro sets, gfortran -E -cpp -P, diff -B after masking the __LINE__ numbers of runtime-error calls. FAIL on
#       any preprocessor error output, an empty preprocessed file, or any remaining difference.
#   (2) ctest (GNU devel noMPI, as in CI): exit code and the "0 tests failed out of N" line are written for verdict.py.
# Writes OUTDIR/{ppcmp.txt,ppcmp.exit,ctest.log,ctest.exit}. The data.h5 identity runs (GNU noMPI HDF5, base vs
# candidate, same seeds) are ordinary acceptance runs (builds gnu_cand / base_gnu) judged by the replay gate.
# Exit code: 0 only if ppcmp and ctest both pass.
set -o pipefail
TICKET=$1; S12=$2; OUT=$3
C=$HOME/alf-gpu/work/t$TICKET/cand/$S12; [ -d $C/src ] && [ -d $C/base_src ] || { echo "need src and base_src"; exit 2; }
[ -e "$OUT" ] && { echo "refusing: $OUT exists"; exit 2; }
mkdir -p "$OUT"; export TMPDIR=$HOME/alf-gpu/work/t$TICKET/tmp; mkdir -p $TMPDIR
. ~/alf-gpu/env/gnu.sh > /dev/null 2>&1
W=$(mktemp -d $TMPDIR/ppcmp.XXXXXX)
# generated headers come from a GNU devel build of the candidate (configure + make all)
cp -a $C/src/. $W/tree/ 2>/dev/null || { mkdir -p $W/tree; cp -a $C/src/. $W/tree/; }; chmod -R u+w $W/tree
( cd $W/tree && . ./configure.sh GNU devel noMPI NO-INTERACTIVE > $OUT/ctest_configure.log 2>&1 && make -j4 all > $OUT/ctest_make.log 2>&1 )
mk=$?
files=$(cd $C && diff -rq base_src src 2>/dev/null | sed -n 's/^Files base_src\/\(.*\) and src\/.* differ$/\1/p' | grep -E '\.(F90|f90|F)$' | grep -v '^Prog/gpu_\|^gpu_bench/\|^testsuite/' )
tot=0; err=0; n=0
{
  echo "changed Fortran files (GPU=none compiled): ${files:-none}"
  for f in $files; do
    for m in "-DMPI -DHDF5 -DHDF5_ZLIB" "-DHDF5 -DHDF5_ZLIB" "-DMPI -DHDF5 -DHDF5_ZLIB -DSTABLOG" "-DMPI -DHDF5 -DHDF5_ZLIB -DQR_PREPIVOT" "-DMPI -DTEMPERING" "-DMPI -DHDF5 -DHDF5_ZLIB -DSTAB3"; do
      for v in base_src src; do
        d=$W/$v; mkdir -p $d/$(dirname $f); cp $C/$v/$f $d/$f; cp $W/tree/Prog/*.h $d/$(dirname $f)/ 2>/dev/null
        ( cd $d/$(dirname $f) && gfortran -E -cpp -P $m $(basename $f) > $W/$v.i 2> $W/$v.err )
        [ -s $W/$v.err ] && { echo "PREPROCESSOR ERROR $v $f [$m]: $(head -2 $W/$v.err)"; err=$((err+1)); }
        [ -s $W/$v.i ] || { echo "EMPTY $v $f [$m]"; err=$((err+1)); }
        sed -E 's/"([A-Za-z_0-9]+\.F90)",[0-9]+/"\1",N/' $W/$v.i > $W/$v.n
      done
      k=$(diff -B $W/base_src.n $W/src.n | wc -l); n=$((n+1)); tot=$((tot+k))
      echo "$f [$m]: $(wc -l < $W/src.i) lines, differing diff lines after __LINE__ masking: $k"
    done
  done
  echo "comparisons $n preprocessor_errors $err differing_lines $tot"
} > $OUT/ppcmp.txt
pp=$(( (err > 0 || tot > 0) ? 1 : 0 )); echo $pp > $OUT/ppcmp.exit
if [ $mk = 0 ]; then
  ( cd $W/tree/testsuite && mkdir -p tests && cd tests && cmake .. > $OUT/ctest_cmake.log 2>&1 && cmake --build . -j4 >> $OUT/ctest_cmake.log 2>&1 && ctest -j4 --output-on-failure > $OUT/ctest.log 2>&1 )
  echo $? > $OUT/ctest.exit
else
  echo "build failed" > $OUT/ctest.log; echo 99 > $OUT/ctest.exit
fi
rm -rf -- "$W"
echo "ppcmp exit $pp ($(tail -1 $OUT/ppcmp.txt)); ctest exit $(cat $OUT/ctest.exit) ($(grep 'tests failed out of' $OUT/ctest.log))"
[ $pp = 0 ] && [ "$(cat $OUT/ctest.exit)" = 0 ]
