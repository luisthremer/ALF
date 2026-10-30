#!/bin/bash
# Ticket 207 (dev4: rounding-order nulls of the paired test, development only), on the LAPTOP from the repo root: rsync the working
# tree to nividia ~/alf-gpu/work/t207/dev4/src/ and the dev4 scripts (plus t207_diagstat.py) to dev4/.  The source state (HEAD and the
# sha256 of the uncommitted diff) goes to dev4/src/.d4_srcstate and is copied into every build's log.
set -e
ssh nividia 'mkdir -p ~/alf-gpu/work/t207/dev4/{src,tmp,logs,bin,runs}'
{ echo "head $(git rev-parse HEAD)"; echo "diff_sha256 $(git diff HEAD | shasum -a 256 | cut -c1-64)"; git status --short; } > .d4_srcstate
rsync -a --delete --exclude .git --exclude .alf_gpu_stamp --exclude HDF5 --exclude gpu_bench --exclude Documentation \
      --exclude doxygen --exclude Mathematica --exclude .DS_Store --exclude '*.o' --exclude '*.mod' --exclude '*.a' \
      --exclude '*.so' --exclude '*.out' --exclude .deps.mk --exclude tests --exclude '*.pdf' ./ nividia:alf-gpu/work/t207/dev4/src/
rm -f .d4_srcstate
git diff HEAD > /dev/null
rsync -a docs/gpu/scripts/nividia/t207_d4_* docs/gpu/scripts/nividia/t207_diagstat.py nividia:alf-gpu/work/t207/dev4/
echo synced
