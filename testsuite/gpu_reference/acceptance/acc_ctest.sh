#!/bin/bash
# Acceptance harness, MACHINE side: build the unit tests of a shipped candidate from its verified tarball and run ctest.
#   usage: acc_ctest.sh TICKET WHO SHA12 OUTDIR TOOL [--base] [--regex R] [--repeat N] [--env VAR=value ...] -- CONFIGURE ARGS...
#     --base     : build the base (base_src.tar) instead of the candidate (for identity comparisons);
#     TOOL: env file (gnu | nvhpc | intel); CONFIGURE ARGS e.g. "GNU devel noMPI" or "IntelLLVM devel noMPI GPU=intel".
#     --regex R  : ctest -R R (default: all tests);   --repeat N : run the selected tests N times (e.g. t_order);
#     --env      : variables for the test runs (e.g. a declared fault mode of a control); ALF_GPU_* and ALF_LA_TEST_*
#                  need the acceptance file to declare the run as a control (the driver passes them only then);
#                  OMP_TARGET_OFFLOAD, ONEAPI_DEVICE_SELECTOR, ZE_AFFINITY_MASK, CUDA_VISIBLE_DEVICES are refused.
# Writes OUTDIR/{stamp, build.log, ctest.log, ctest.exit, summary}: ctest.exit is 0 (all passed) or 1 (a test failed
# its own check: "***Failed"); a build/setup failure or any crash, exception, timeout, "Not Run" or abort writes 2
# (INVALID, never a pass and never a detection). Tests run with ctest -V (their output is in ctest.log).
# GPU configurations (GPU=nvidia|intel) run with OMP_TARGET_OFFLOAD=MANDATORY (GPU=intel, which needs TOOL intel:
# also ONEAPI_DEVICE_SELECTOR=level_zero:0, exported after the env file and asserted; an env file that sets
# ZE_AFFINITY_MASK, ZE_ENABLE_PCI_ID_DEVICE_ORDER, ZE_FLAT_DEVICE_HIERARCHY or SYCL_DEVICE_FILTER makes the repetition
# INVALID; the stamp line "selector" records it, env_rep_N.txt the effective environment, checked after each
# repetition: offload, selector and threads 1 as stamped, else INVALID), one test at a time, all tests with OMP/MKL/OpenBLAS
# threads 1, under flock ~/alf-gpu/gpu.lock. summary: "<regex>: P of N repetitions
# passed" and the last "tests failed out of" line. Exit code = the value in ctest.exit.
set -o pipefail
export LC_ALL=C
TICKET=$1; WHO=$2; S12=$3; OUT=$4; TOOL=$5; shift 5
RX=""; REP=1; ENVX=(); SRCT=src; KEY=tar_sha256
while [ $# -gt 0 ] && [ "$1" != -- ]; do
  case "$1" in --base) SRCT=base_src; KEY=base_tar_sha256; shift;; --regex) RX=$2; shift 2;; --repeat) REP=$2; shift 2;;
    --env) case "$2" in OMP_TARGET_OFFLOAD=*|ONEAPI_DEVICE_SELECTOR=*|ZE_AFFINITY_MASK=*|ZE_ENABLE_PCI_ID_DEVICE_ORDER=*|ZE_FLAT_DEVICE_HIERARCHY=*|SYCL_DEVICE_FILTER=*|CUDA_VISIBLE_DEVICES=*)
             echo "refusing to override $2"; exit 2;; esac
           ENVX+=("$2"); shift 2;;
    *) echo "bad option $1"; exit 2;; esac
done
[ "${1:-}" = -- ] && shift
C=$HOME/alf-gpu/work/t$TICKET/cand/$WHO/$S12; [ -f $C/.alf_candidate ] || { echo "no shipped candidate $C"; exit 2; }
want=$(awk -v k=$KEY '$1==k {print $2}' $C/.alf_candidate); have=$(sha256sum $C/$SRCT.tar | cut -d' ' -f1)
[ -n "$want" ] && [ "$want" = "$have" ] || { echo "TAR HASH MISMATCH"; exit 2; }
ENVF=$HOME/alf-gpu/env/$TOOL.sh; [ -f "$ENVF" ] || { echo "no env file $ENVF"; exit 2; }
case " $* " in *" GPU=intel "*) [ "$TOOL" = intel ] || { echo "refusing: GPU=intel needs TOOL intel (got $TOOL)"; exit 2; };; esac
mkdir "$OUT" || { echo "refusing: $OUT exists"; exit 2; }
GPUENV=(); OFF=none; JOBS=1; LOCK=(); SELV=""   # one test at a time: deterministic output order for identity comparisons
case " $* " in *" GPU=nvidia "*|*" GPU=intel "*) GPUENV=(OMP_TARGET_OFFLOAD=MANDATORY); OFF=MANDATORY; JOBS=1
  LOCK=(flock "$HOME/alf-gpu/gpu.lock");; esac
case " $* " in *" GPU=intel "*) SELV=ONEAPI_DEVICE_SELECTOR=level_zero:0;; esac
{ grep -E '^(sha|who|base_sha) ' $C/.alf_candidate; echo "$KEY $want"; echo "source $SRCT"; echo "offload $OFF"; S2=${SELV#ONEAPI_DEVICE_SELECTOR=}; echo "selector ${S2:-none}"; echo "tool $TOOL"; echo "configure $*"; echo "regex ${RX:-all}"
  echo "repeat $REP"; echo "env ${ENVX[*]:-none}"; echo "start $(date -Is)"; } > "$OUT/stamp"
export TMPDIR=$HOME/alf-gpu/work/t$TICKET/tmp; mkdir -p $TMPDIR
W=$(mktemp -d $TMPDIR/ctest.XXXXXX) || exit 2
tar -xf $C/$SRCT.tar -C $W || { echo 2 > "$OUT/ctest.exit"; exit 2; }
env -i HOME="$HOME" USER="$USER" LANG=C LC_ALL=C PATH=/usr/bin:/bin TMPDIR="$TMPDIR" bash --noprofile --norc -c '
  # readonly copies (an env file may reset the positionals); NO-FALLBACK: never a silent default (CPU) machine
  readonly __acc_w="$0" __acc_envf="$1"; readonly -a __acc_args=("${@:2}")
  . "$__acc_envf" > /dev/null 2>&1 || exit 2
  builtin cd "$__acc_w" || exit 2
  . ./configure.sh "${__acc_args[@]}" NO-FALLBACK NO-INTERACTIVE && make -j4 all && cd testsuite && mkdir -p tests && cd tests && cmake .. && cmake --build . -j4' \
  "$W" "$ENVF" "$@" > "$OUT/build.log" 2>&1
if [ $? != 0 ]; then echo 2 > "$OUT/ctest.exit"; echo "BUILD FAILED (INVALID)" > "$OUT/summary"; rm -rf -- "$W"; exit 2; fi
pass=0; fails=0; crash=0
for i in $(seq 1 "$REP"); do
  "${LOCK[@]}" env -i HOME="$HOME" USER="$USER" LANG=C LC_ALL=C PATH=/usr/bin:/bin TMPDIR="$TMPDIR" "${GPUENV[@]}" "${ENVX[@]}" \
    bash --noprofile --norc -c '
      # readonly copies (an env file may reset the positionals or reuse short names); builtin/command against functions
      readonly __acc_d="$0" __acc_envf="$1" __acc_rx="$2" __acc_j="$3" __acc_sel="$4" __acc_envout="$5" __acc_off="$6" __acc_elog="$7"
      . "$__acc_envf" > "$__acc_elog" 2>&1 || builtin exit 2
      command grep -q "readonly variable\|readonly function\|cannot unset" "$__acc_elog" && { builtin echo "env file tried to change the harness variables: refusing"; builtin exit 2; }
      if [ -n "$__acc_off" ]; then
        builtin export "$__acc_off"; [ "$OMP_TARGET_OFFLOAD" = MANDATORY ] || { builtin echo "offload not in effect: refusing"; builtin exit 2; }
      fi
      if [ -n "$__acc_sel" ]; then
        [ -n "${ZE_AFFINITY_MASK+x}${ZE_ENABLE_PCI_ID_DEVICE_ORDER+x}${ZE_FLAT_DEVICE_HIERARCHY+x}${SYCL_DEVICE_FILTER+x}" ] && { builtin echo "env file sets a device-ordering variable: refusing"; builtin exit 2; }
        builtin export "$__acc_sel"; [ "$ONEAPI_DEVICE_SELECTOR" = level_zero:0 ] || { builtin echo "device selector not in effect: refusing"; builtin exit 2; }
      fi
      builtin export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
      command env | command sort > "$__acc_envout"
      builtin cd "$__acc_d" && command ctest -j "$__acc_j" -V ${__acc_rx:+-R "$__acc_rx"}' \
    "$W/testsuite/tests" "$ENVF" "$RX" "$JOBS" "$SELV" "$OUT/env_rep_$i.txt" "${GPUENV[0]:-}" "$OUT/env_file_rep_$i.log" > "$OUT/rep_$i.log" 2>&1
  rc=$?
  # the pinning as it took effect (checked outside the inner shell): a missing or wrong env makes the result INVALID
  E="$OUT/env_rep_$i.txt"; envok=1
  [ -f "$E" ] || envok=0
  for v in OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; do [ $envok = 1 ] && { grep -qx "$v" "$E" || envok=0; }; done
  [ $envok = 1 ] && [ "$OFF" = MANDATORY ] && { grep -qx 'OMP_TARGET_OFFLOAD=MANDATORY' "$E" || envok=0; }
  [ $envok = 1 ] && [ -n "$SELV" ] && { grep -qx 'ONEAPI_DEVICE_SELECTOR=level_zero:0' "$E" && \
    ! grep -q -E '^(ZE_AFFINITY_MASK|ZE_ENABLE_PCI_ID_DEVICE_ORDER|ZE_FLAT_DEVICE_HIERARCHY|SYCL_DEVICE_FILTER)=' "$E" || envok=0; }
  [ $envok = 1 ] || { crash=1; echo "repetition $i: effective environment not as stamped (offload/selector)" >> "$OUT/env_check.txt"; }
  # the temporary build path and backtrace addresses (ASLR, e.g. from the expected-failure tests) differ between runs
  sed -i -E "s#$W#<BUILD>#g; s/0x[0-9a-f]{6,}/0xADDR/g; s/^[0-9]+: //" "$OUT/rep_$i.log"   # also the ctest -V test-index prefix
  cat "$OUT/rep_$i.log" >> "$OUT/ctest.log"
  # a check failure = every failed test printed its own "FAILED <name>" line (la_test_util convention) and nothing
  # looks like a crash; anything else that exits nonzero is a crash (never a detection)
  nf=$(grep -c '\*\*\*Failed' "$OUT/rep_$i.log" || true); nc=$(grep -c -E '(^|: )FAILED [A-Za-z_]' "$OUT/rep_$i.log" || true)
  if grep -q -E '\*\*\*(Exception|Timeout|Not Run)|SEGFAULT|Child aborted|Subprocess aborted|Illegal|forrtl: severe|Fortran runtime error|Program received signal|Segmentation fault|omptarget|terminate called' "$OUT/rep_$i.log"; then crash=1
  elif [ $rc = 0 ]; then pass=$((pass+1))
  elif [ "$nf" -gt 0 ] && [ "$nc" -ge "$nf" ]; then fails=$((fails+1))
  else crash=1; fi
done
ntests=$(grep -c -E "^[0-9]+% tests passed, [0-9]+ tests failed out of" "$OUT/ctest.log" || true)
{ echo "${RX:-all}: $pass of $REP repetitions passed"; grep "tests failed out of" "$OUT/ctest.log" | tail -1; } > "$OUT/summary"
if [ "$ntests" != "$REP" ] || [ $crash = 1 ]; then echo 2 > "$OUT/ctest.exit"; echo "INVALID: $ntests ctest summaries for $REP repetitions, crash=$crash" >> "$OUT/summary"; [ -f "$OUT/env_check.txt" ] && cat "$OUT/env_check.txt" >> "$OUT/summary"
elif [ $pass = "$REP" ]; then echo 0 > "$OUT/ctest.exit"; else echo 1 > "$OUT/ctest.exit"; fi
echo "end $(date -Is)" >> "$OUT/stamp"
rm -rf -- "$W"
cat "$OUT/summary"; exit "$(cat "$OUT/ctest.exit")"
