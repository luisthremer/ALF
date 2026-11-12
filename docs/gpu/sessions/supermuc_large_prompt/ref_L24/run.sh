#!/bin/bash
. ~/alf-gpu/env/nvhpc.sh; export TMPDIR=~/alf-gpu/work/prodmpi/tmp
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=0
export CUDA_MPS_PIPE_DIRECTORY=$HOME/alf-gpu/work/mps/pipe CUDA_MPS_LOG_DIRECTORY=$HOME/alf-gpu/work/mps/log
nvidia-cuda-mps-control -d; sleep 2
echo "start $(date)"; mpirun -n 4 --bind-to core ~/alf-gpu/work/prodmpi/bin/ALF-gpu-mpi.out > run.log 2>&1; echo "alf rc=$? end $(date)"
echo quit | nvidia-cuda-mps-control; sleep 3; echo "MPS left: $(pgrep -u $USER -f nvidia-cuda-mps | wc -l)"
