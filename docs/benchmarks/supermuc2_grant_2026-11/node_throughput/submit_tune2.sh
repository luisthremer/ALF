#!/bin/bash
# Tuning stage 2: one parameter at a time from the best layout of stage 1 (disjoint core slices), beta = 4, L = 24, 40, 56.
# Variants (tag <layout>_<variant>): nd32/128/256 = N_delay; Intel-opt. only: qr8/32/64 = ALFSH_QRNB, bf3/bf8 = ALFSH_BLKFUSE,
# bpf0 = ALFSH_BLKPF, asy0 = ALFSH_ASYNCBLK; both codes: tl2/4, tr1/4, vb4/16 = ALF_GPU_TL/TR/VB (need ALF_GPU_DEBUG=1). Results in runs_tune2/.
set -e
A=/hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27
export BETA=4 RUNDIR=runs_tune2 EXE_GPU=$A/build_pp_alfgpu_detfix/Prog/ALF.out BINB=$A/bin_batched_65800b82/pp
declare -A BEST=( [gpub_24]=4_2 [gpub_40]=7_2 [gpub_56]=8_1 [gpu_24]=4_3 [gpu_40]=7_2 [gpu_56]=7_2 )
COMMON="nd32:NDELAY=32 nd128:NDELAY=128 nd256:NDELAY=256 tl2:ALF_GPU_DEBUG=1,ALF_GPU_TL=2 tl4:ALF_GPU_DEBUG=1,ALF_GPU_TL=4 tr1:ALF_GPU_DEBUG=1,ALF_GPU_TR=1 tr4:ALF_GPU_DEBUG=1,ALF_GPU_TR=4 vb4:ALF_GPU_DEBUG=1,ALF_GPU_VB=4 vb16:ALF_GPU_DEBUG=1,ALF_GPU_VB=16"
INTEL="qr8:ALFSH_QRNB=8 qr32:ALFSH_QRNB=32 qr64:ALFSH_QRNB=64 bf3:ALFSH_BLKFUSE=3 bf8:ALFSH_BLKFUSE=8 bpf0:ALFSH_BLKPF=0 asy0:ALFSH_ASYNCBLK=0"
for L in 56 40 24; do
  nb=6; [ $L = 56 ] && nb=4
  t=00:50:00; [ $L = 56 ] && t=01:30:00
  for side in gpub gpu; do
    lay=${BEST[${side}_$L]}; items=""; n=0
    vars="$COMMON"; [ $side = gpub ] && vars="$COMMON $INTEL"
    for v in $vars; do items="$items $side:$L:$lay:$nb:${v%%:*}:PINMODE=excl,${v#*:}"; n=$((n+1)); done
    sbatch --export=BETA,RUNDIR,EXE_GPU,BINB -N $n -J nt_t2_${side}_$L -p general -t $t multi2.sh $items
  done
done
