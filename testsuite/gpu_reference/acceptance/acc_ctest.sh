#!/bin/bash
# Acceptance harness, MACHINE side: build the unit tests of a shipped candidate from its verified tarball and run ctest.
#   usage: acc_ctest.sh TICKET WHO SHA12 OUTDIR TOOL [--base] [--regex R] [--repeat N] [--env VAR=value ...] -- CONFIGURE ARGS...
#     --base     : build the base (base_src.tar) instead of the candidate (for identity comparisons);
#     TOOL: env file (gnu | nvhpc | intel); CONFIGURE ARGS e.g. "GNU devel noMPI" or "IntelLLVM devel noMPI GPU=intel".
#     --regex R  : ctest -R R (default: all tests);   --repeat N : run the selected tests N times (e.g. t_order);
#     --env      : variables for the test runs (e.g. a declared fault mode of a control); ALF_GPU_* need the
#                  acceptance file to declare the run as a control (the driver passes them only then).
# Writes OUTDIR/{stamp, build.log, ctest.log, ctest.exit, summary}: ctest.exit is 0 (all passed) or 1 (a test failed);
# a build or setup failure writes 2 (INVALID, never a pass and never a detection). summary: "<regex>: P of N repetitions
# passed" and the last "tests failed out of" line. Exit code = the value in ctest.exit.
set -o pipefail
export LC_ALL=C
TICKET=$1; WHO=$2; S12=$3; OUT=$4; TOOL=$5; shift 5
RX=""; REP=1; ENVX=(); SRCT=src; KEY=tar_sha256
while [ $# -gt 0 ] && [ "$1" != -- ]; do
  case "$1" in --base) SRCT=base_src; KEY=base_tar_sha256; shift;; --regex) RX=$2; shift 2;; --repeat) REP=$2; shift 2;; --env) ENVX+=("$2"); shift 2;; *) echo "bad option $1"; exit 2;; esac
done
[ "${1:-}" = -- ] && shift
C=$HOME/alf-gpu/work/t$TICKET/cand/$WHO/$S12; [ -f $C/.alf_candidate ] || { echo "no shipped candidate $C"; exit 2; }
want=$(awk -v k=$KEY '$1==k {print $2}' $C/.alf_candidate); have=$(sha256sum $C/$SRCT.tar | cut -d' ' -f1)
[ -n "$want" ] && [ "$want" = "$have" ] || { echo "TAR HASH MISMATCH"; exit 2; }
ENVF=$HOME/alf-gpu/env/$TOOL.sh; [ -f "$ENVF" ] || { echo "no env file $ENVF"; exit 2; }
mkdir "$OUT" || { echo "refusing: $OUT exists"; exit 2; }
{ grep -E '^(sha|tar_sha256|base_sha|who) ' $C/.alf_candidate; echo "source $SRCT"; echo "tool $TOOL"; echo "configure $*"; echo "regex ${RX:-all}"
  echo "repeat $REP"; echo "env ${ENVX[*]:-none}"; echo "start $(date -Is)"; } > "$OUT/stamp"
export TMPDIR=$HOME/alf-gpu/work/t$TICKET/tmp; mkdir -p $TMPDIR
W=$(mktemp -d $TMPDIR/ctest.XXXXXX) || exit 2
tar -xf $C/$SRCT.tar -C $W || { echo 2 > "$OUT/ctest.exit"; exit 2; }
env -i HOME="$HOME" USER="$USER" LANG=C LC_ALL=C PATH=/usr/bin:/bin TMPDIR="$TMPDIR" bash --noprofile --norc -c '
  . "$1" > /dev/null 2>&1 || exit 2; shift; cd "$0" || exit 2
  . ./configure.sh "$@" NO-INTERACTIVE && make -j4 all && cd testsuite && mkdir -p tests && cd tests && cmake .. && cmake --build . -j4' \
  "$W" "$ENVF" "$@" > "$OUT/build.log" 2>&1
if [ $? != 0 ]; then echo 2 > "$OUT/ctest.exit"; echo "BUILD FAILED (INVALID)" > "$OUT/summary"; rm -rf -- "$W"; exit 2; fi
pass=0; fails=0
for i in $(seq 1 "$REP"); do
  env -i HOME="$HOME" USER="$USER" LANG=C LC_ALL=C PATH=/usr/bin:/bin TMPDIR="$TMPDIR" "${ENVX[@]}" bash --noprofile --norc -c '
    . "$1" > /dev/null 2>&1 || exit 2; cd "$0" && ctest -j4 --output-on-failure ${2:+-R "$2"}' "$W/testsuite/tests" "$ENVF" "$RX" >> "$OUT/ctest.log" 2>&1
  rc=$?
  case $rc in 0) pass=$((pass+1));; 8) fails=$((fails+1));; *) [ $rc -ne 0 ] && fails=$((fails+1));; esac
done
ntests=$(grep -c "tests passed" "$OUT/ctest.log" || true)
{ echo "${RX:-all}: $pass of $REP repetitions passed"; grep "tests failed out of" "$OUT/ctest.log" | tail -1; } > "$OUT/summary"
if [ "$ntests" != "$REP" ]; then echo 2 > "$OUT/ctest.exit"; echo "INVALID: $ntests ctest summaries for $REP repetitions" >> "$OUT/summary"
elif [ $pass = "$REP" ]; then echo 0 > "$OUT/ctest.exit"; else echo 1 > "$OUT/ctest.exit"; fi
echo "end $(date -Is)" >> "$OUT/stamp"
rm -rf -- "$W"
cat "$OUT/summary"; exit "$(cat "$OUT/ctest.exit")"
