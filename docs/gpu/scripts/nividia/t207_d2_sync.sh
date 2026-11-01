#!/bin/bash
# Ticket 207 (rebase on 212, development only), on the LAPTOP from the repo root: rsync the working tree to nividia
# ~/alf-gpu/work/t207/dev2/src/ and the dev2 scripts to dev2/.  "base" argument: also export the base commit 3f8ec8fa
# (gpu-port with 206 + 212) to dev2/base_src/ (GPU=none identity, non-interference and before/after timing).
set -e
ssh nividia 'mkdir -p ~/alf-gpu/work/t207/dev2/{src,base_src,tmp,logs,bin,runs}'
rsync -a --delete --exclude .git --exclude .alf_gpu_stamp --exclude HDF5 --exclude gpu_bench --exclude Documentation \
      --exclude doxygen --exclude Mathematica --exclude .DS_Store --exclude '*.o' --exclude '*.mod' --exclude '*.a' \
      --exclude '*.so' --exclude '*.out' --exclude .deps.mk --exclude tests --exclude '*.pdf' ./ nividia:alf-gpu/work/t207/dev2/src/
rsync -a docs/gpu/scripts/nividia/t207_d2_* docs/gpu/scripts/nividia/t207_firstdiv.py docs/gpu/scripts/nividia/t207_diagstat.py \
      nividia:alf-gpu/work/t207/dev2/
for a in "$@"; do
  if [ "$a" = base ]; then
    ssh nividia 'rm -rf -- ~/alf-gpu/work/t207/dev2/base_src && mkdir -p ~/alf-gpu/work/t207/dev2/base_src'
    git archive 3f8ec8fa | ssh nividia 'tar -x -C ~/alf-gpu/work/t207/dev2/base_src'
  fi
done
echo synced
