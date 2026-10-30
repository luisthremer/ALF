#!/bin/bash
# repeat t_order of build_TAG: t302a_reps.sh TAG NREP OUTDIR [HOOKVAR]  (on intel; per-rep logs, exit codes and times stay in OUTDIR)
TAG=$1; NREP=$2; OUT=$3; HOOK=$4
W=$HOME/alf-gpu/work/t302a/dev; mkdir "$OUT" || exit 2
. $HOME/alf-gpu/env/intel.sh 2>/dev/null
export ONEAPI_DEVICE_SELECTOR=level_zero:0 OMP_TARGET_OFFLOAD=MANDATORY
[ -n "$HOOK" ] && export $HOOK=1
cd $W/build_$TAG/testsuite/tests/alf_la.tests
pass=0
for i in $(seq 1 $NREP); do
  s=$(date +%s.%N)
  ./t_order > $OUT/rep_$i.out 2>&1; rc=$?
  e=$(date +%s.%N)
  echo "$i $rc $(echo "$e - $s" | bc)" >> $OUT/reps.txt
  [ $rc = 0 ] && grep -q "^PASSED t_order" $OUT/rep_$i.out && pass=$((pass+1))
done
echo "^t_order\$: $pass of $NREP repetitions passed" > $OUT/summary
echo done >> $OUT/summary
