#!/bin/bash
# Ticket 207 (rebase on 212, development only), nividia: T0 of acceptance 3 (paired same-trajectory test with its two code
# controls, precision means) and 4 (compare2), dev seeds (devN = what T1 uses for accN).
#   usage: t207_d2_t0c.sh RUN gpu|cpu  -> dev2/logs/t0c_<RUN>_<gpu|cpu>.txt, runs in dev2/runs/t0c<RUN>_*
W=$HOME/alf-gpu/work/t207/dev2; R=$1; K=$2; out=$W/logs/t0c_${R}_$K.txt
case "$R" in ""|*/*|*..*) echo bad run label; exit 1;; esac
r() { local id=$1 set=$2; shift 2; echo "$(SET=$set TAG=t0c${R}_$id bash $W/t207_d2_run.sh "$@")" >> $out; }
PD="ALF_GPU_DEBUG=1 ALF_GPU_TAUDIAG=2 ALF_GPU_TAUDIAG_FROM=331"
if [ "$K" = cpu ]; then
  for b in nv_none nv_none_pp nv_none_log; do r pr_${b}_b4 dev1 $b hub12x12_cb_p; r pr_${b}_b10 dev1 $b hub12x12_cb_b10_p; done
  r c2_cpu dev2 nv_none s8_U4_b2_cb_t0; r c2_cpu_null dev3 nv_none s8_U4_b2_cb_t0
elif [ "$K" = gpu ]; then
  r pd_def_b4 dev1 nv_def_pp hub12x12_cb_p $PD
  r pd_log_b4 dev1 nv_log hub12x12_cb_p $PD
  r ctl_pd_tauswap8_log dev1 nv_log hub12x12_cb_p $PD ALF_GPU_TAUSWAP=8
  r ctl_pd_tauns0_def dev1 nv_def_pp hub12x12_cb_p $PD ALF_GPU_TAU_NS=0
  for b in nv_def nv_log; do r pr_${b}_b4 dev1 $b hub12x12_cb_p; r pr_${b}_b10 dev1 $b hub12x12_cb_b10_p; done
  r ctl_pr_ns0_def_b4 dev1 nv_def hub12x12_cb_p ALF_GPU_DEBUG=1 ALF_GPU_NS=0
  r c2_gpu dev1 nv_def s8_U4_b2_cb_t0
fi
echo done >> $out
