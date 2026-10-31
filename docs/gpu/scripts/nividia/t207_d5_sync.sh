#!/bin/bash
# Ticket 207 (dev5: accuracy test against a quad reference, development only), on the LAPTOP from the repo root: rsync the working
# tree to nividia ~/alf-gpu/work/t207/dev5/src/ and the dev5 scripts (plus t207_diagstat.py) to dev5/.  The source state (HEAD and the
# sha256 of the uncommitted diff) goes to dev5/src/.d5_srcstate and is copied into every build's log.
set -e
ssh nividia 'mkdir -p ~/alf-gpu/work/t207/dev5/{src,tmp,logs,bin,runs}'
{ echo "head $(git rev-parse HEAD)"; echo "diff_sha256 $(git diff HEAD | shasum -a 256 | cut -c1-64)"; git status --short; } > .d5_srcstate
rsync -a --delete --exclude .git --exclude .alf_gpu_stamp --exclude HDF5 --exclude gpu_bench --exclude Documentation \
      --exclude doxygen --exclude Mathematica --exclude .DS_Store --exclude '*.o' --exclude '*.mod' --exclude '*.a' \
      --exclude '*.so' --exclude '*.out' --exclude .deps.mk --exclude tests --exclude '*.pdf' ./ nividia:alf-gpu/work/t207/dev5/src/
rm -f .d5_srcstate
git diff HEAD > /dev/null
rsync -a docs/gpu/scripts/nividia/t207_d5_* docs/gpu/scripts/nividia/t207_diagstat.py nividia:alf-gpu/work/t207/dev5/
echo synced
