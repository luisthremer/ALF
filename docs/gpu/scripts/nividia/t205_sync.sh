#!/bin/bash
# Ticket 205, run on the LAPTOP from the repo root: rsync the working tree to nividia ~/alf-gpu/work/t205/src_new/
# ("head" argument: also export git HEAD to src_head/ for the GPU=none comparison; "h5" argument: copy the HDF5 builds
# of ~/alf-gpu/src/ALF-*/HDF5 to work/t205/h5 once, no rebuild).
set -e
R=nividia:alf-gpu/work/t205
ssh nividia 'mkdir -p ~/alf-gpu/work/t205/{src_new,src_head,tmp,logs,bin,runs,h5}'
rsync -a --delete --exclude .git --exclude .alf_gpu_stamp --exclude HDF5 --exclude docs --exclude gpu_bench --exclude Documentation \
      --exclude doxygen --exclude Mathematica --exclude .DS_Store \
      --exclude '*.o' --exclude '*.mod' --exclude '*.a' --exclude '*.so' --exclude '*.out' --exclude .deps.mk \
      --exclude tests --exclude '*.pdf' ./ $R/src_new/
rsync -a docs/gpu/scripts/nividia/t205_* $R/
for a in "$@"; do
  if [ "$a" = head ]; then
    ssh nividia 'rm -rf -- ~/alf-gpu/work/t205/src_head && mkdir -p ~/alf-gpu/work/t205/src_head'
    git archive HEAD | ssh nividia 'tar -x -C ~/alf-gpu/work/t205/src_head --exclude=docs --exclude=Documentation --exclude=gpu_bench --exclude=Mathematica'
  fi
  if [ "$a" = h5 ]; then
    ssh nividia 'W=$HOME/alf-gpu/work/t205; for d in $HOME/alf-gpu/src/ALF-nvhpc/HDF5/nvfortran_26.9-0 $HOME/alf-gpu/src/ALF-gnu/HDF5/GCC_version_13.3.0; do cp -a "$d" "$W/h5/"; done; du -sh "$W/h5"'
  fi
done
echo synced
