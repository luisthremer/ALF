#!/bin/bash
# Ticket 207, nividia: run several binaries concurrently on one prepared config (identical seeds), then print the precision lines.
#   usage: t207_batch.sh CONFIG BIN [BIN ...]     (<= 4 binaries: 4 cores)   -> runs/BIN/CONFIG, table on stdout
W=$HOME/alf-gpu/work/t207
C=$1; shift
for b in "$@"; do TAG=$b bash $W/t207_run.sh $b $C & done; wait
bash $W/t207_prectab.sh $C "$@"
