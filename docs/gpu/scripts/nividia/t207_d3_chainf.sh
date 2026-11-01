#!/bin/bash
# Ticket 207 dev3, nividia: after the variant / amplitude runs: self-tests and timing of fix candidates (TAUSWAP), the
# ALF_GPU_TAUDIAG_SEED demonstration (short runs), then the null-variant + candidate chains on dev1 / dev2 (binary nv_def_pp_f).
#   usage: t207_d3_chainf.sh [short|long|all]   -> dev3/logs/chainf.txt
W=$HOME/alf-gpu/work/t207/dev3; cd $W || exit 1; out=$W/logs/chainf.txt; B=nv_def_pp_f
r() { local tag=$1 set=$2 cfg=$3; shift 3; SET=$set TAG=$tag bash $W/t207_d3_run.sh $B $cfg "$@" >> $out; }
P=${1:-all}
if [ "$P" = short ] || [ "$P" = all ]; then
  for m in 0 8 512 21 1047; do
    r f_st4_m$m dev1 hub12x12_cb ALF_GPU_EXIT=1 ALF_GPU_SELFTEST=4 ALF_GPU_TAUSWAP=$m
    r f_tim_m$m dev1 hub12x12_cb_tm ALF_GPU_EXIT=1 ALF_GPU_DEBUG=1 ALF_GPU_TAUTIME=1 ALF_GPU_TAUSWAP=$m
    r f_tim20_m$m dev1 hub20x20_cb_tm ALF_GPU_EXIT=1 ALF_GPU_DEBUG=1 ALF_GPU_TAUTIME=1 ALF_GPU_TAUSWAP=$m
  done
  D="ALF_GPU_DEBUG=1 ALF_GPU_TAUDIAG=2 ALF_GPU_TAUDIAG_FROM=1"
  r f_seed1a smoke1 hub12x12_cb_p $D ALF_GPU_TAUDIAG_SEED=1
  r f_seed2 smoke1 hub12x12_cb_p $D ALF_GPU_TAUDIAG_SEED=2
  r f_seed1b smoke1 hub12x12_cb_p $D ALF_GPU_TAUDIAG_SEED=1
  r f_seedu1a smoke1 hub12x12_cb_p $D
  r f_seedu1b smoke1 hub12x12_cb_p $D
  r f_seedu2 smoke2 hub12x12_cb_p $D
fi
if [ "$P" = long ] || [ "$P" = all ]; then
  V="ALF_GPU_DEBUG=1 ALF_GPU_TAUDIAG=2 ALF_GPU_TAUDIAG_FROM=331 ALF_GPU_TAUDIAG_SEED=20261029 ALF_GPU_TAUDIAG_VAR=23,1047,2071,3095,8,512"
  for s in dev1 dev2; do r f_v3_$s $s hub12x12_cb_p $V; done
fi
echo "done $P" >> $out
