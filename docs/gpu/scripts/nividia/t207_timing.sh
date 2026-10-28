#!/bin/bash
# Ticket 207, nividia: time of one tau pass (without observables) device vs host (1 core), ALF_GPU_TAUTIME=1, and s/sweep.
#   usage: t207_timing.sh tau BIN CONFIG [CONFIG ...]     -> logs/timing_tau_BIN.txt   (BIN = nv_def | nv_log; production flags, no QR_PREPIVOT)
#          t207_timing.sh sweep BIN CONFIG [CONFIG ...]   -> logs/timing_sweep_BIN.txt (wall time of a short run / number of sweeps; with Ltau = 1)
# GPU runs one after the other (one GPU process at a time on the A30).
W=$HOME/alf-gpu/work/t207
M=$1; B=$2; shift 2
case $M in
  tau)
    out=$W/logs/timing_tau_$B.txt; : > "$out"
    for c in "$@"; do
      nvidia-smi --query-gpu=memory.used --format=csv,noheader >> "$out"
      r=$(TAG=tt_$B timeout 900 bash $W/t207_run.sh $B $c ALF_GPU_TAUTIME=1 ALF_GPU_EXIT=1)
      echo "== $r" >> "$out"
      grep -E -A3 "GPU tau timing" $W/runs/tt_$B/$c/run.log | cut -c1-300 >> "$out"
    done
    cat "$out";;
  sweep)
    out=$W/logs/timing_sweep_$B.txt; : > "$out"
    for c in "$@"; do
      nsw=$(grep -i -E "^Nsweep" $W/runs/base/$c/parameters | head -1 | sed -E 's/.*= *([0-9]+).*/\1/')
      nb=$(grep -i -E "^Nbin " $W/runs/base/$c/parameters | head -1 | sed -E 's/.*= *([0-9]+).*/\1/')
      nvidia-smi --query-gpu=memory.used --format=csv,noheader >> "$out"
      r=$(TAG=ts_$B timeout 900 bash $W/t207_run.sh $B $c)
      echo "$r  sweeps=$((nsw*nb))" >> "$out"
    done
    cat "$out";;
esac
