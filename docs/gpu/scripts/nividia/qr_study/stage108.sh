#!/bin/bash
# usage: stage108.sh NAME TARGET_S NBIN SEEDS "cases" builds...   -> prepares NAME_prep, runs NAME_run (JOBS env, default 8)
export TMPDIR=$HOME/alf-gpu/work/t108/tmp
. ~/alf-gpu/env/gnu.sh
cd ~/alf-gpu/work/t108
name=$1; target=$2; nbin=$3; seeds=$4; cases=$5; shift 5
specs=$(python qs/mkcases108.py $target $nbin $seeds "$cases")
python qs/prep.py ${name}_prep src $specs > ${name}_prep.log 2>&1
jobs=()
for s in $specs; do for b in "$@"; do jobs+=("${s%%:*}:$b"); done; done
# longest first (large L / beta) for load balance
printf '%s\n' "${jobs[@]}" | sort -t_ -k1.2,1.3nr -k2.2,2nr > ${name}_jobs.txt
qs/run108.sh ${name}_prep ${name}_run ${JOBS:-8} $(cat ${name}_jobs.txt)
