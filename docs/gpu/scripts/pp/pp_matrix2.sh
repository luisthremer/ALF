#!/bin/bash
cd ~/alf-gpu/work/pp; B=~/alf-gpu/work/pp/bin
TAG=optskip ALF_GPU_DEBUG=1 ALF_GPU_OBS_SKIP=1,3 PP_BIN=$B/ALF-gpu-pp-opt.out bash bench_v2.sh 3 10 25 64
TAG=optall ALF_GPU_DEBUG=1 ALF_GPU_OBS_SKIP=1,3 PP_NWRAP=20 PP_BIN=$B/ALF-gpu-pp-opt.out bash bench_v2.sh 3 10 25 64
