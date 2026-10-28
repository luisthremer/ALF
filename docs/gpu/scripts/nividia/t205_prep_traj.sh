#!/bin/bash
# Ticket 205, nividia: run directories for the trajectory test (1 sweep per bin, NBIN bins) in work/t205/runs/traj.
#   usage: t205_prep_traj.sh [NBIN]
W=$HOME/alf-gpu/work/t205; export TMPDIR=$W/tmp; . ~/alf-gpu/env/gnu.sh
cd $W/src_new/testsuite/gpu_reference || exit 1
: "${HOME:?}"
rm -rf -- "$W/runs/traj"; mkdir -p "$W/runs"
python3 prep.py $W/runs/traj hub6x6_cb hub12x12_cb --nsweep 1 --nbin ${1:-20} --alf-dir $W/src_new
