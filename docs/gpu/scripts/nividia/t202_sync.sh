#!/bin/bash
# Ticket 202, run on the LAPTOP from the repo root: rsync the source tree to nividia ~/alf-gpu/work/t202/src/.
# Excludes .git, docs, HDF5, build products.
set -e
R=nividia:alf-gpu/work/t202
ssh nividia 'mkdir -p ~/alf-gpu/work/t202/src ~/alf-gpu/work/t202/tmp ~/alf-gpu/work/t202/logs'
rsync -a --delete --exclude .git --exclude HDF5 --exclude docs --exclude gpu_bench --exclude Documentation \
      --exclude doxygen --exclude Mathematica --exclude .DS_Store \
      --exclude '*.o' --exclude '*.mod' --exclude '*.a' --exclude '*.so' --exclude '*.out' --exclude .deps.mk \
      --exclude tests --exclude '*.pdf' ./ $R/src/
echo synced
