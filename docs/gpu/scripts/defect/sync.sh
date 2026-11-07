#!/bin/bash
# usage: sync.sh gpu|cpu   -> rsync the working tree to nividia:~/alf-gpu/work/fg1/src_<x>
SRC=/Users/luis/Documents/_WORK/Coding/CLAUDE_ALF_GPU_SUBMATRIX_POINT_DEFECT/ALF_GPU_SUBMATRIX/
ssh nividia "mkdir -p ~/alf-gpu/work/fg1/src_$1 ~/alf-gpu/work/fg1/bin ~/alf-gpu/work/fg1/logs ~/alf-gpu/work/fg1/tmp"
rsync -a --delete --exclude .git --exclude '*.o' --exclude '*.mod' --exclude '*.out' --exclude Documentation --exclude doxygen --exclude gpu_bench --exclude docs --exclude Mathematica --exclude testsuite --exclude Analysis --exclude HDF5 $SRC nividia:~/alf-gpu/work/fg1/src_$1/
