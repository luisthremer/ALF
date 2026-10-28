#!/bin/bash
# Ticket 203, nividia: build GPU=nvidia, then GPU=none in the SAME tree without 'make clean' (Prog/.alf_gpu_stamp check),
# then back to nvidia.  Reports which Prog objects were rebuilt and whether the binary contains gpu_state symbols.
W=$HOME/alf-gpu/work/t203; export TMPDIR=$W/tmp; : "${HOME:?}"
rm -rf -- "$W/tree_switch"; mkdir -p "$W/tree_switch"; rsync -a "$W/src_new/" "$W/tree_switch/"; cd "$W/tree_switch" || exit 1
step() {  # $1 = label, rest = configure args
  l=$1; shift
  bash -c ". ~/alf-gpu/env/nvhpc.sh; . ./configure.sh PGI noMPI $* NO-INTERACTIVE > $W/logs/sw_conf_$l.log 2>&1 && make -j4 program > $W/logs/sw_make_$l.log 2>&1"
  echo "$l: make rc=$? Prog objects compiled: $(grep -c 'pgfortran -c' $W/logs/sw_make_$l.log), gpu_state compiled: $(grep -c 'gpu_state_mod.F90' $W/logs/sw_make_$l.log), gpu_state symbols in ALF.out: $(nm Prog/ALF.out | grep -ci gpu_state), stamp=[$(cat Prog/.alf_gpu_stamp)]"
  tail -2 $W/logs/sw_make_$l.log | head -1
}
step nvidia1 GPU=nvidia
step none    GPU=none
step nvidia2 GPU=nvidia
