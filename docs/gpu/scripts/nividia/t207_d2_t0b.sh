#!/bin/bash
# Ticket 207 (rebase on 212, development only), nividia: T0 production runs (dev seeds) for same-seed vs CPU oracle, run-twice
# determinism, "the tau pass leaves the sweep unchanged" (Ltau = 1 vs Ltau = 0, equal-time datasets) and non-interference vs the
# base binary (3f8ec8fa).  usage: t207_d2_t0b.sh RUN gpu|cpu   -> dev2/logs/t0b_<RUN>_<gpu|cpu>.txt, runs in dev2/runs/t0b<RUN>_*
W=$HOME/alf-gpu/work/t207/dev2; R=$1; K=$2; out=$W/logs/t0b_${R}_$K.txt
case "$R" in ""|*/*|*..*) echo bad run label; exit 1;; esac
r() { echo "$(TAG=t0b${R}_$1 bash $W/t207_d2_run.sh $2 $3)" >> $out; }
if [ "$K" = cpu ]; then
  r cpp6 nv_none_pp hub6x6_cb_fd; r clog6 nv_none_log hub6x6_cb_fd
  r cpp12 nv_none_pp hub12x12_cb_fd; r clog12 nv_none_log hub12x12_cb_fd
else
  for b in nv_def nv_log; do
    r ${b}_6a $b hub6x6_cb_fd; r ${b}_6b $b hub6x6_cb_fd; r ${b}_6lt0 $b hub6x6_cb_fd_lt0
    r ${b}_b10a $b hub12x12_cb_b10_s; r ${b}_b10b $b hub12x12_cb_b10_s; r ${b}_b10lt0 $b hub12x12_cb_b10_s_lt0
    r base_${b}_6 base_$b hub6x6_cb_fd
  done
  r nv_def_12 nv_def hub12x12_cb_fd; r nv_log_12 nv_log hub12x12_cb_fd
fi
echo done >> $out
