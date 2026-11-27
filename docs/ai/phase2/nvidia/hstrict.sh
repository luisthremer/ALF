#!/bin/bash -l
# variant of a candidate without -Mfprelaxed -fast (accuracy diagnosis), GPU=nvidia noMPI
T=$1; P=~/p2_nvidia; D=$P/cand_${T}_strict
rm -rf $D; mkdir $D
(cd $P/cand_$T && tar cf - --exclude=testsuite/tests --exclude='*.o' --exclude='*.mod' --exclude='*.a' --exclude=ALF.out .) | tar xf - -C $D
cd $D
. ./configure.sh HELMAGPU noMPI GPU=nvidia NO-INTERACTIVE > cfg.log 2>&1
for v in ALF_FLAGS_MODULES ALF_FLAGS_PROG ALF_FLAGS_ANA ALF_FLAGS_QRREF; do
  eval "export $v=\"\$(echo \"\$$v\" | sed -e 's/ -Mfprelaxed//g' -e 's/ -fast//g')\""
done
echo "$ALF_FLAGS_PROG"
make clean >/dev/null 2>&1; make -j16 lib program > build.log 2>&1 && echo strict ok || tail -3 build.log
