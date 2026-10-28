#!/bin/bash
# Ticket 207, nividia: GPU=none identity of the program output: GNU noMPI HDF5 builds of HEAD (src_head) and of the working
# tree (src_new), same seeds, 3 configs (20 sweeps x 2 bins: hub6x6_cb, hub6x6_dense, hub6x6_cb_symm), data.h5 bit for bit.
#   usage: t207_ident.sh   (needs bin/ALF_gnu_head.out and bin/ALF_gnu_new.out, built with
#          t207_build.sh src_head gnu_head gnu GNU noMPI HDF5   and   t207_build.sh src_new gnu_new gnu GNU noMPI HDF5)
W=$HOME/alf-gpu/work/t207
for c in hub6x6_cb hub6x6_dense hub6x6_cb_symm; do
  TAG=id_head bash $W/t207_run.sh gnu_head $c
  TAG=id_new bash $W/t207_run.sh gnu_new $c
  bash $W/t207_h5cmp.sh $W/runs/id_head/$c/data.h5 $W/runs/id_new/$c/data.h5
  if diff <(grep -v -i "cpu\|time\|date\|git\|version" $W/runs/id_head/$c/info) <(grep -v -i "cpu\|time\|date\|git\|version" $W/runs/id_new/$c/info) > /dev/null; then
    echo "$c: info (without time/git/version lines) identical"
  else
    echo "$c: info differs"
  fi
done
