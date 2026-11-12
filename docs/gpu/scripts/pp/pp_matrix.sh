#!/bin/bash
cd ~/alf-gpu/work/pp; B=~/alf-gpu/work/pp/bin
pgrep -u luis -f nvidia-cuda-mps-control > /dev/null || { export CUDA_MPS_PIPE_DIRECTORY=$HOME/alf-gpu/work/mps/pipe CUDA_MPS_LOG_DIRECTORY=$HOME/alf-gpu/work/mps/log CUDA_VISIBLE_DEVICES=0; nvidia-cuda-mps-control -d; }
TAG=base PP_BIN=$B/ALF-gpu-pp.out bash bench_v2.sh 3 10 25 64
TAG=opt PP_BIN=$B/ALF-gpu-pp-opt.out bash bench_v2.sh 3 10 25 64
TAG=optskip ALF_GPU_OBS_SKIP=1,3 PP_BIN=$B/ALF-gpu-pp-opt.out bash bench_v2.sh 3 10 25 64
TAG=optnw20 PP_NWRAP=20 PP_BIN=$B/ALF-gpu-pp-opt.out bash bench_v2.sh 3 10 25 64
TAG=optall ALF_GPU_OBS_SKIP=1,3 PP_NWRAP=20 PP_BIN=$B/ALF-gpu-pp-opt.out bash bench_v2.sh 3 10 25 64
TAG=base2 PP_BIN=$B/ALF-gpu-pp.out bash bench_v2.sh 3 10 25 64
