#!/bin/bash
# Ticket 206, nividia: run several binaries concurrently on one prepared config (identical seeds), then print the precision lines.
#   usage: t206_batch.sh CONFIG BIN [BIN ...]     (<= 4 binaries: 4 cores)   -> runs/BIN/CONFIG, table on stdout
W=$HOME/alf-gpu/work/t206
C=$1; shift
for b in "$@"; do TAG=$b bash $W/t206_run.sh $b $C & done; wait
bash $W/t206_prectab.sh $C "$@"
