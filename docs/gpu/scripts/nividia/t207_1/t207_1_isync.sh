#!/bin/bash
# Ticket 207-1 (development only), on the LAPTOP from the checkout root: rsync the working tree and the intel scripts to
# intel:~/alf-gpu/work/t207-1/ .   usage: t207_1_isync.sh SCRIPTDIR
set -e
SD=$1
ssh intel 'mkdir -p ~/alf-gpu/work/t207-1/{src,tmp,logs,bin,runs}'
{ echo "head $(git rev-parse HEAD)"; echo "diff_sha256 $(git diff HEAD | shasum -a 256 | cut -c1-64)"; git status --short; } > .srcstate
rsync -a --delete --exclude .git --exclude .alf_gpu_stamp --exclude HDF5 --exclude gpu_bench --exclude Documentation \
      --exclude doxygen --exclude Mathematica --exclude .DS_Store --exclude '*.o' --exclude '*.mod' --exclude '*.a' \
      --exclude '*.so' --exclude '*.out' --exclude .deps.mk --exclude tests --exclude '*.pdf' ./ intel:alf-gpu/work/t207-1/src/
rm -f .srcstate
rsync -a $SD/t207_1_i*.sh intel:alf-gpu/work/t207-1/
echo synced
