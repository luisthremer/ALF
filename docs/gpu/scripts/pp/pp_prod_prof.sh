#!/bin/bash
# production run (see pp_prod.sh) with nsys profiling of all 9 ranks for a window after start-up; the run continues unprofiled afterwards
R=~/alf-gpu/work/pp; export TMPDIR=$R/tmp
. ~/alf-gpu/env/nvhpc.sh
export CUDA_VISIBLE_DEVICES=0 CUDA_MPS_PIPE_DIRECTORY=$HOME/alf-gpu/work/mps/pipe CUDA_MPS_LOG_DIRECTORY=$HOME/alf-gpu/work/mps/log
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1 CUDA_MPS_ACTIVE_THREAD_PERCENTAGE=25
d=$R/runs/L24_prod
[ -e $d ] && { echo "$d exists"; exit 1; }
~/alf-gpu/venv/bin/python $R/pp_setup.py $d 3 10 100000 64 23.8 || exit 1
cd $d; mkdir -p prof; date > STARTED
mpirun -n 9 --bind-to core nsys profile -o prof/rank%q{OMPI_COMM_WORLD_RANK} -t cuda,nvtx,osrt --cuda-memory-usage=false --sample=none --cpuctxsw=none --delay=${PROF_DELAY:-100} --duration=${PROF_DUR:-150} --kill=none --force-overwrite=true $R/bin/ALF-gpu-pp.out > run.log 2>&1; echo "rc=$?" > DONE; date >> DONE
echo quit | nvidia-cuda-mps-control
