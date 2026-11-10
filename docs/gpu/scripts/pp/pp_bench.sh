#!/bin/bash
# usage: pp_bench.sh NPER NSWEEP [MPS_PCT]   (3 parameter sets x NPER walkers)
R=~/alf-gpu/work/pp; export TMPDIR=$R/tmp
. ~/alf-gpu/env/nvhpc.sh
export CUDA_VISIBLE_DEVICES=0 CUDA_MPS_PIPE_DIRECTORY=$HOME/alf-gpu/work/mps/pipe CUDA_MPS_LOG_DIRECTORY=$HOME/alf-gpu/work/mps/log
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1
[ -n "$3" ] && export CUDA_MPS_ACTIVE_THREAD_PERCENTAGE=$3
d=$R/runs/bench_${TAG:-x}_n$1_p${3:-100}_d${4:-32}; python3 -c "import shutil,os;shutil.rmtree('$d') if os.path.isdir('$d') else None"
~/alf-gpu/venv/bin/python $R/pp_setup.py $d $1 $2 1 ${4:-32} > /dev/null || exit 1
cd $d; s=$(date +%s.%N)
mpirun -n $((3*$1)) --bind-to core ${PP_BIN:-$R/bin/ALF-gpu-pp.out} > run.log 2>&1; rc=$?
e=$(date +%s.%N)
echo "nper=$1 ranks=$((3*$1)) pct=${3:-100} nsweep=$2 rc=$rc wall=$(echo "$e - $s" | bc) s"
tail -2 run.log | cut -c1-150; grep -h -E "GPU path|Precision Green|Acceptance" Temp_0/info | head -3
