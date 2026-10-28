#!/bin/bash
# Ticket 205, run on the LAPTOP: sync + rebuild GPU=nvidia (nv_def) on nividia.  usage: t205_cycle.sh
cd "$(dirname "$0")/../../../.." || exit 1
bash docs/gpu/scripts/nividia/t205_sync.sh >/dev/null && ssh nividia 'cd ~/alf-gpu/work/t205 && bash t205_build.sh src_new nv_def nvhpc PGI noMPI HDF5 GPU=nvidia' 2>&1 | tail -15
