#!/bin/bash
# Ticket 203, run on the LAPTOP from the repo root: rsync the working tree to nividia ~/alf-gpu/work/t203/src_new/
# and (once, "head" argument) export git HEAD to src_head/ for the GPU=none identity comparison.
set -e
R=nividia:alf-gpu/work/t203
ssh nividia 'mkdir -p ~/alf-gpu/work/t203/{src_new,src_head,tmp,logs,bin,runs}'
rsync -a --delete --exclude .git --exclude .alf_gpu_stamp --exclude HDF5 --exclude docs --exclude gpu_bench --exclude Documentation \
      --exclude doxygen --exclude Mathematica --exclude .DS_Store \
      --exclude '*.o' --exclude '*.mod' --exclude '*.a' --exclude '*.so' --exclude '*.out' --exclude .deps.mk \
      --exclude tests --exclude '*.pdf' ./ $R/src_new/
if [ "$1" = head ]; then
  ssh nividia 'rm -rf -- ~/alf-gpu/work/t203/src_head && mkdir -p ~/alf-gpu/work/t203/src_head'
  git archive HEAD | ssh nividia 'tar -x -C ~/alf-gpu/work/t203/src_head --exclude=docs --exclude=Documentation --exclude=gpu_bench --exclude=Mathematica'
fi
echo synced
