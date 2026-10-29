#!/bin/bash
# Acceptance harness, MACHINE side: GPU=none guarantee of a shipped candidate against its base (both from git archive).
#   usage: acc_gpunone.sh TICKET WHO SHA12 OUTDIR [--plant]
#   (1) preprocessed sources: every Fortran file that git reports as modified between base and candidate
#       (cand/<sha12>/changed_files.txt from acc_ship.sh), except Prog/gpu_*.F90 (compiled only with ALF_GPU),
#       gpu_bench/ and testsuite/; 6 macro sets; gfortran -E -cpp -P with each side's OWN generated headers;
#       diff -B after masking the __LINE__ numbers of runtime-error calls. Macro sets: plain (noMPI, no HDF5, as in CI),
#       MPI+HDF5, HDF5, STABLOG, QR_PREPIVOT, TEMPERING, STAB3, PARALLEL_PARAMS, OBS_LEGACY (gfortran -E; ifx-only
#       branches such as __INTEL_COMPILER_BUILD_DATE are not covered by this check).
#       Counted against the verdict: preprocessor errors, empty outputs, differing lines, added or deleted CPU-path
#       Fortran files, changed build files (Makefile*, configure.sh, CMakeLists.txt, *.mk), and other changed files
#       outside docs/, testsuite/, gpu_bench/, .claude/ and Prog/gpu_* (e.g. Hamiltonians.list, scripts, C/C++, headers);
#       these are counts that the acceptance file must expect exactly (reviewed changes only).
#   In --plant mode a failed plant or any copy/preprocessor error makes the run INVALID (exit value 2, no summary
#   line), so it can never count as a detection.
#   (2) ctest (GNU devel noMPI, as in CI) of the candidate: exit code and the "tests failed out of N" line.
#   --plant: positive control - one token is appended to the candidate copy of the first compared file; the
#       comparison must then report differing lines (ctest is skipped).
# Writes OUTDIR/{stamp,ppcmp.txt,ppcmp.exit,ctest.log,ctest.exit}. Exit 0 only if ppcmp and ctest both pass.
set -o pipefail
export LC_ALL=C
TICKET=$1; WHO=$2; S12=$3; OUT=$4; PLANT=${5:-}
C=$HOME/alf-gpu/work/t$TICKET/cand/$WHO/$S12
[ -d $C/src ] && [ -d $C/base_src ] && [ -f $C/changed_files.txt ] || { echo "need src, base_src, changed_files.txt"; exit 2; }
mkdir "$OUT" || { echo "refusing: $OUT exists"; exit 2; }
grep -E '^(sha|base_sha) ' $C/.alf_candidate > "$OUT/stamp"; echo "plant ${PLANT:-no}" >> "$OUT/stamp"
export TMPDIR=$HOME/alf-gpu/work/t$TICKET/tmp; mkdir -p $TMPDIR
W=$(mktemp -d $TMPDIR/gpunone.XXXXXX) || exit 2
# GNU devel builds of both sides: generated headers (and the candidate's ctest)
for v in base_src src; do   # from the tarballs after a fresh hash check, not from the extracted trees
  key=tar_sha256; [ $v = base_src ] && key=base_tar_sha256
  want=$(awk -v k=$key '$1==k {print $2}' $C/.alf_candidate); have=$(sha256sum $C/$v.tar | cut -d' ' -f1)
  [ -n "$want" ] && [ "$want" = "$have" ] || { echo "TAR HASH MISMATCH $v"; exit 2; }
  mkdir $W/tree_$v $W/raw_$v && tar -xf $C/$v.tar -C $W/tree_$v && tar -xf $C/$v.tar -C $W/raw_$v || exit 2
  env -i HOME="$HOME" USER="$USER" LANG=C LC_ALL=C PATH=/usr/bin:/bin TMPDIR="$TMPDIR" bash --noprofile --norc -c '
    . "$HOME/alf-gpu/env/gnu.sh" > /dev/null 2>&1 || exit 2; cd "$0" || exit 2
    . ./configure.sh GNU devel noMPI NO-INTERACTIVE > "$1/configure_$2.log" 2>&1 && make -j4 all > "$1/make_$2.log" 2>&1' \
    $W/tree_$v "$OUT" $v
  echo $? > "$OUT/build_$v.exit"
done
mod=$(awk '$1=="M" || $1=="T" {print $2}' $C/changed_files.txt | grep -E '\.(F90|f90|F|f)$' | grep -v -E '^(Prog/gpu_|gpu_bench/|testsuite/)')
addel=$(awk '$1=="A" || $1=="D" {print $1" "$2}' $C/changed_files.txt | grep -E '\.(F90|f90|F|f)$' | grep -v -E ' (Prog/gpu_|gpu_bench/|testsuite/)' | wc -l)
build=$(awk '{print $2}' $C/changed_files.txt | grep -E '(^|/)(Makefile[^/]*|configure\.sh|CMakeLists\.txt|[^/]*\.mk)$' | wc -l)
otherl=$(awk '{print $2}' $C/changed_files.txt | grep -v -E '^(docs/|testsuite/|gpu_bench/|\.claude/|Prog/gpu_)' \
        | grep -v -E '\.(F90|f90|F|f)$' | grep -v -E '(^|/)(Makefile[^/]*|configure\.sh|CMakeLists\.txt|[^/]*\.mk)$')
other=$(echo -n "$otherl" | grep -c . || true)
# the plant control needs a compared file: without changed CPU-path Fortran it compares Prog/main.F90 (identical on both sides)
[ "$PLANT" = --plant ] && [ -z "$mod" ] && mod=Prog/main.F90
tot=0; err=0; n=0; first=1
{
  echo "compared files: ${mod:-none}"
  echo "added_or_deleted_cpu_fortran $addel"
  echo "changed_build_files $build: $(awk '{print $2}' $C/changed_files.txt | grep -E '(^|/)(Makefile[^/]*|configure\.sh|CMakeLists\.txt|[^/]*\.mk)$' | tr '\n' ' ')"
  for f in $mod; do
    for m in "-UNONE" "-DMPI -DHDF5 -DHDF5_ZLIB" "-DHDF5 -DHDF5_ZLIB" "-DMPI -DHDF5 -DHDF5_ZLIB -DSTABLOG" "-DMPI -DHDF5 -DHDF5_ZLIB -DQR_PREPIVOT" "-DMPI -DTEMPERING" "-DMPI -DHDF5 -DHDF5_ZLIB -DSTAB3" "-DMPI -DPARALLEL_PARAMS" "-DMPI -DOBS_LEGACY"; do
      for v in base_src src; do
        d=$W/pp_$v; rm -rf -- "$W/pp_$v"; mkdir -p $d/$(dirname $f)
        cp $W/raw_$v/$f $d/$f && chmod u+w $d/$f || { echo "COPY ERROR $v $f"; err=$((err+1)); }
        if [ "$PLANT" = --plant ] && [ $v = src ] && [ $first = 1 ]; then
          echo "integer, parameter :: planted_control_token = 1" >> $d/$f && grep -q planted_control_token $d/$f \
            || { echo "PLANT FAILED $f"; err=$((err+1)); }   # a plant that did not happen must not look like a pass
        fi
        cp $W/tree_$v/Prog/*.h $d/$(dirname $f)/ || { echo "HEADER COPY ERROR $v"; err=$((err+1)); }
        ( cd $d/$(dirname $f) && gfortran -E -cpp -P $m $(basename $f) > $W/$v.i 2> $W/$v.err )
        [ -s $W/$v.err ] && { echo "PREPROCESSOR ERROR $v $f [$m]: $(head -2 $W/$v.err)"; err=$((err+1)); }
        [ -s $W/$v.i ] || { echo "EMPTY $v $f [$m]"; err=$((err+1)); }
        sed -E 's/"([A-Za-z_0-9]+\.F90)",[0-9]+/"\1",N/' $W/$v.i > $W/$v.n
      done
      k=$(diff -B $W/base_src.n $W/src.n | wc -l); n=$((n+1)); tot=$((tot+k))
      echo "$f [$m]: $(wc -l < $W/src.i) lines, differing diff lines after __LINE__ masking: $k"
    done
    first=0
  done
  echo "other_changed_files $other: $(echo $otherl)"
  if [ "$PLANT" = --plant ] && [ $err -gt 0 ]; then
    echo "PLANT RUN INVALID: $err copy/plant/preprocessor errors (no summary line)"
  else
    echo "comparisons $n preprocessor_errors $err differing_lines $tot added_or_deleted $addel build_files $build other_files $other"
  fi
} > "$OUT/ppcmp.txt"
# added/deleted CPU-path Fortran files and build-file changes are reported as counts; the acceptance file states the
# expected counts (verdict: exact match), so an unexpected one is a FAIL
pp=$(( (err > 0 || tot > 0) ? 1 : 0 ))
if [ "$PLANT" = --plant ] && [ $err -gt 0 ]; then pp=2; fi
echo $pp > "$OUT/ppcmp.exit"
if [ "$PLANT" = --plant ]; then echo "skipped (plant control)" > "$OUT/ctest.log"; echo 2 > "$OUT/ctest.exit"
elif [ "$(cat $OUT/build_src.exit)" = 0 ]; then
  ( cd $W/tree_src/testsuite && mkdir -p tests && cd tests && cmake .. > "$OUT/ctest_cmake.log" 2>&1 \
    && cmake --build . -j4 >> "$OUT/ctest_cmake.log" 2>&1 && ctest -j4 --output-on-failure > "$OUT/ctest.log" 2>&1 )
  rc=$?   # normalised like acc_ctest.sh: 0 passed, 1 a test failed its check, 2 invalid (crash, timeout, build)
  if grep -q -E '\*\*\*(Exception|Timeout|Not Run)|SEGFAULT|Child aborted' "$OUT/ctest.log"; then echo 2 > "$OUT/ctest.exit"
  elif [ $rc = 0 ]; then echo 0 > "$OUT/ctest.exit"
  elif grep -q '\*\*\*Failed' "$OUT/ctest.log"; then echo 1 > "$OUT/ctest.exit"
  else echo 2 > "$OUT/ctest.exit"; fi
else echo "build failed" > "$OUT/ctest.log"; echo 2 > "$OUT/ctest.exit"; fi
rm -rf -- "$W"
echo "ppcmp exit $pp ($(tail -1 "$OUT/ppcmp.txt")); ctest exit $(cat "$OUT/ctest.exit") ($(grep 'tests failed out of' "$OUT/ctest.log"))"
[ $pp = 0 ] && [ "$(cat "$OUT/ctest.exit")" = 0 ]
