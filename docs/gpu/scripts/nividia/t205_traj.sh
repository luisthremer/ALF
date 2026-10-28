#!/bin/bash
# Ticket 205, nividia: trajectory test (validation level 2).  usage: t205_traj.sh CONFIG
#   runs/traj/<CONFIG> (1 sweep per bin, 20 bins, from t205_prep.sh-like prep) with
#     nv_none  : GPU=none nvfortran build                       -> runs/tr_none/<CONFIG>
#     nv_def   : GPU=nvidia, ALF_GPU_NOUPDATE=1 (host path)     -> runs/tr_host/<CONFIG>, trace_cpu
#     nv_def   : GPU=nvidia, device update path                 -> runs/tr_gpu/<CONFIG>,  trace_gpu
#   then compare traces (first differing WRAPGRUP/DO call and site), data.h5 none vs host (bit-identical expected) and
#   none vs gpu (per-bin Ener_scal).
W=$HOME/alf-gpu/work/t205
C=$1
BASE=traj TAG=tr_none bash $W/t205_run.sh nv_none $C
BASE=traj TAG=tr_host bash $W/t205_run.sh nv_def $C ALF_GPU_NOUPDATE=1 ALF_GPU_TRACE=trace_cpu
BASE=traj TAG=tr_gpu  bash $W/t205_run.sh nv_def $C ALF_GPU_TRACE=trace_gpu
. ~/alf-gpu/env/gnu.sh
python3 -I $W/t205_trace_cmp.py $W/runs/tr_host/$C/trace_cpu $W/runs/tr_gpu/$C/trace_gpu 40
bash $W/t205_h5cmp.sh $W/runs/tr_none/$C/data.h5 $W/runs/tr_host/$C/data.h5
python3 -I $W/t205_ener.py $W/runs/tr_none/$C/data.h5 $W/runs/tr_gpu/$C/data.h5 $W/runs/tr_gpu/$C/trace_gpu $W/runs/tr_host/$C/trace_cpu
