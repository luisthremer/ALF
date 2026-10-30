#!/bin/bash
# Ticket 207 (dev3: diagnosis of the default-path paired spread, development only), on the LAPTOP from the repo root: rsync the working tree to nividia
# ~/alf-gpu/work/t207/dev3/src/ and the dev3 scripts (plus t207_diagstat.py) to dev3/.
set -e
ssh nividia 'mkdir -p ~/alf-gpu/work/t207/dev3/{src,tmp,logs,bin,runs}'
rsync -a --delete --exclude .git --exclude .alf_gpu_stamp --exclude HDF5 --exclude gpu_bench --exclude Documentation \
      --exclude doxygen --exclude Mathematica --exclude .DS_Store --exclude '*.o' --exclude '*.mod' --exclude '*.a' \
      --exclude '*.so' --exclude '*.out' --exclude .deps.mk --exclude tests --exclude '*.pdf' ./ nividia:alf-gpu/work/t207/dev3/src/
rsync -a docs/gpu/scripts/nividia/t207_d3_* docs/gpu/scripts/nividia/t207_diagstat.py \
      nividia:alf-gpu/work/t207/dev3/
echo synced
