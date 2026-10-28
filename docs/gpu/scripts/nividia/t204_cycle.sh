#!/bin/bash
# Ticket 204, run on the LAPTOP: sync + rebuild GPU=nvidia (nv_def) on nividia.  usage: t204_cycle.sh
cd "$(dirname "$0")/../../../.." || exit 1
bash docs/gpu/scripts/nividia/t204_sync.sh >/dev/null && ssh nividia 'cd ~/alf-gpu/work/t204 && bash t204_build.sh src_new nv_def nvhpc PGI noMPI HDF5 GPU=nvidia' 2>&1 | tail -15
