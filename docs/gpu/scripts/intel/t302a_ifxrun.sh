#!/bin/bash
# second half of t302a_ifxid.sh: run hub6x6_cb (shortened: NSweep/NBin as argv 1/2) with base and candidate, compare data.h5/info
W=$HOME/alf-gpu/work/t302a/dev; NS=${1:-10}; NB=${2:-10}
. $HOME/alf-gpu/venv/bin/activate
R=$W/ifxid_runs; rm -rf $R; mkdir -p $R
cd $W/src/testsuite/gpu_reference && python3 prep.py $R/prep hub6x6_cb --nsweep $NS --nbin $NB --alf-dir $W/src > $R/prep.log 2>&1 || { tail $R/prep.log; exit 1; }
. $HOME/alf-gpu/env/intel.sh >/dev/null 2>&1
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
for v in base cand; do
  mkdir -p $R/$v; cp -a $R/prep/hub6x6_cb/. $R/$v/
  ( cd $R/$v; ( time $W/ifx_$v/Prog/ALF.out > run.log 2>&1 ) 2> time.txt; echo "$v exit $?: $(grep real time.txt)" )
done
cmp $R/base/data.h5 $R/cand/data.h5 && echo "data.h5 IDENTICAL ($(sha256sum $R/base/data.h5 | cut -c1-16))" || echo "data.h5 DIFFER"
cmp $R/base/info $R/cand/info > /dev/null && echo "info identical" || diff <(grep -v -i "time\|date\|wall" $R/base/info) <(grep -v -i "time\|date\|wall" $R/cand/info) | head
