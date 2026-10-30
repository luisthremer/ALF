#!/bin/bash
# Ticket 212 (development only), on the LAPTOP from the repo root: rsync the working tree to nividia ~/alf-gpu/work/t212/dev/src/.
# "base" argument: also export the base commit 1a363070 to dev/base_src/ (for the cost comparison and the GPU=none comparison).
set -e
ssh nividia 'mkdir -p ~/alf-gpu/work/t212/dev/{src,base_src,tmp,logs,bin,runs}'
rsync -a --delete --exclude .git --exclude .alf_gpu_stamp --exclude HDF5 --exclude gpu_bench --exclude Documentation \
      --exclude doxygen --exclude Mathematica --exclude .DS_Store --exclude '*.o' --exclude '*.mod' --exclude '*.a' \
      --exclude '*.so' --exclude '*.out' --exclude .deps.mk --exclude tests --exclude '*.pdf' ./ nividia:alf-gpu/work/t212/dev/src/
for a in "$@"; do
  if [ "$a" = base ]; then
    ssh nividia 'rm -rf -- ~/alf-gpu/work/t212/dev/base_src && mkdir -p ~/alf-gpu/work/t212/dev/base_src'
    git archive 1a363070 | ssh nividia 'tar -x -C ~/alf-gpu/work/t212/dev/base_src'
  fi
done
echo synced
