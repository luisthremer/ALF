#!/bin/bash
# Ticket 205, run on the LAPTOP from the repo root: rsync the tree to intel:~/alf-gpu/work/t205/src_new/ and the scripts, and
# copy the prepared run directories (parameters + seeds; tiny) from nividia (argument "runs").
set -e
R=intel:alf-gpu/work/t205
ssh intel 'mkdir -p ~/alf-gpu/work/t205/{src_new,tmp,logs,bin,runs}'
rsync -a --delete --exclude .git --exclude .alf_gpu_stamp --exclude HDF5 --exclude docs --exclude gpu_bench --exclude Documentation \
      --exclude doxygen --exclude Mathematica --exclude .DS_Store \
      --exclude '*.o' --exclude '*.mod' --exclude '*.a' --exclude '*.so' --exclude '*.out' --exclude .deps.mk \
      --exclude tests --exclude '*.pdf' ./ $R/src_new/
rsync -a docs/gpu/scripts/nividia/t205_intel_build.sh docs/gpu/scripts/nividia/t205_intel_run.sh $R/
if [ "$1" = runs ]; then
  ssh nividia 'cd ~/alf-gpu/work/t205/runs && tar cf - base/hub6x6_cb base/hub6x6_dense base/hub6x6_cb_flux' | ssh intel 'cd ~/alf-gpu/work/t205/runs && tar xf -'
fi
echo synced
