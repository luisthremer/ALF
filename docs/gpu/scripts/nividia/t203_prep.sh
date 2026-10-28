#!/bin/bash
# Ticket 203, nividia: create the run directories (parameters + one seeds file shared by every run) in work/t203/runs/base.
W=$HOME/alf-gpu/work/t203; export TMPDIR=$W/tmp; . ~/alf-gpu/env/gnu.sh
cd $W/src_new/testsuite/gpu_reference || exit 1
: "${HOME:?}"
rm -rf -- "$W/runs/base"; mkdir -p "$W/runs"
python3 prep.py $W/runs/base hub6x6_cb hub6x6_dense hub12x12_cb hub6x6_cb_mzoff hub6x6_cb_symm --nsweep ${1:-20} --nbin ${2:-2} --alf-dir $W/src_new
