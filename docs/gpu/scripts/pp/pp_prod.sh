#!/bin/bash
# L=24 production: U = 3.0, 5.25, 8.0 with 3 MPI walkers each (PARALLEL_PARAMS), GPU submatrix N_delay 64, CUDA MPS 25 %.
R=~/alf-gpu/work/pp; export TMPDIR=$R/tmp
. ~/alf-gpu/env/nvhpc.sh
export CUDA_VISIBLE_DEVICES=0 CUDA_MPS_PIPE_DIRECTORY=$HOME/alf-gpu/work/mps/pipe CUDA_MPS_LOG_DIRECTORY=$HOME/alf-gpu/work/mps/log
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1 CUDA_MPS_ACTIVE_THREAD_PERCENTAGE=25
export PP_NWRAP=20
d=$R/runs/L24_prod
[ -e $d ] && { echo "$d exists"; exit 1; }
~/alf-gpu/venv/bin/python $R/pp_setup.py $d 3 10 100000 64 23.8 || exit 1
cd $d; date > STARTED
mpirun -n 9 --bind-to core $R/bin/ALF-gpu-pp.out > run.log 2>&1; echo "rc=$?" > DONE; date >> DONE
echo quit | nvidia-cuda-mps-control
