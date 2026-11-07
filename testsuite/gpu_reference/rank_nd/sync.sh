#!/bin/bash
# Ticket 401 (development): on the LAPTOP from the repo root: rsync the working tree to HOST:~/alf-gpu/work/t401b/src/.
#   usage: testsuite/gpu_reference/rank_nd/sync.sh nividia|intel
set -e
H=$1; case "$H" in nividia|intel) ;; *) echo "usage: sync.sh nividia|intel"; exit 1;; esac
ssh $H 'mkdir -p ~/alf-gpu/work/t401b/{src,tmp,logs,bin,runs}'
rsync -a --delete --exclude .git --exclude .alf_gpu_stamp --exclude HDF5 --exclude gpu_bench --exclude Documentation \
      --exclude doxygen --exclude Mathematica --exclude .DS_Store --exclude '*.o' --exclude '*.mod' --exclude '*.a' \
      --exclude '*.so' --exclude '*.out' --exclude .deps.mk --exclude tests --exclude '*.pdf' ./ $H:alf-gpu/work/t401b/src/
echo synced $H
