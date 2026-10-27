#!/bin/bash
export TMPDIR=$HOME/alf-gpu/work/t108/tmp
. ~/alf-gpu/env/gnu.sh
cd ~/alf-gpu/work/t108
C="6:4:4:10 6:10:4:10 6:20:4:10 6:4:8:5 6:10:8:5 8:10:4:10 8:20:4:10 8:10:8:5 8:20:8:5 12:10:4:10 12:10:4:5 12:20:4:5 12:10:8:5 6:20:8:5:dense"
rm -rf main_run/*_s2 main_run/*_s3
rm -rf main2_prep; mkdir main2_prep
for d in main_prep/*_s1; do cp -r $d main2_prep/; done
qs/reseed108.sh main2_prep "$C"
for c in $C; do
  IFS=: read L b U w d <<< "$c"; n=L${L}_b${b}_U${U}_w${w}${d:+_dense}
  for s in 2 3; do for v in up v1 v2 v3n32 v3n8; do echo ${n}_s$s:$v; done; done
done | sort -t_ -k1.2,1.3nr > main2_jobs.txt
wc -l main2_jobs.txt
qs/run108.sh main2_prep main_run 8 $(cat main2_jobs.txt)
