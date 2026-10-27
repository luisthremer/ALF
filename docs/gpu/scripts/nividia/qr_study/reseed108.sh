#!/bin/bash
# Ticket 108: pyALF writes the SAME seeds file for every case, so _s1/_s2/_s3 of stage "main" were identical runs.
# This builds real extra seeds: _s2/_s3 get the seed from line 11/51 of the s1 seeds file as their first line.
# usage: reseed108.sh PREPDIR "L:beta:U:Nwrap[:dense] ..."   (parameters copied from the _s1 dir)
prep=$1
for c in $2; do
  IFS=: read L b U w d <<< "$c"
  n=L${L}_b${b}_U${U}_w${w}${d:+_dense}
  for s in 2 3; do
    mkdir -p $prep/${n}_s$s && cp $prep/${n}_s1/parameters $prep/${n}_s$s/
    if [ $s = 2 ]; then k=11; else k=51; fi
    { sed -n "${k}p" $prep/${n}_s1/seeds; sed -n '2,$p' $prep/${n}_s1/seeds; } > $prep/${n}_s$s/seeds
  done
done
