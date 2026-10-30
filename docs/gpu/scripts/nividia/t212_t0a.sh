#!/bin/bash
# Ticket 212 (development only), nividia: T0 self-tests, faults, refusals, invalid values (all short runs).
# usage: t212_t0a.sh   -> dev/logs/t0a.txt (one line per run: rc and key lines), runs in dev/runs/t0a/<id>
W=$HOME/alf-gpu/work/t212/dev; out=$W/logs/t0a.txt; : > $out
r() { # id binary config env...
  id=$1; b=$2; c=$3; shift 3
  res=$(TAG=t0a_$id bash $W/t212_run.sh $b $c "$@")
  echo "== $id: $res env: $*" >> $out
  grep -hE "PASSED|FAILED|GPU debug env|GPU fault|GPU build|refused|invalid" $W/runs/t0a_$id/$c/run.log | cut -c1-300 | sort -u >> $out
}
E="ALF_GPU_EXIT=1"
for c in hub6x6_cb; do r st2_$c nv_def $c $E ALF_GPU_SELFTEST=2; done
for c in hub6x6_cb hub6x6_cb_flux; do r upd_$c nv_def $c $E ALF_GPU_UPDTEST=1; done
for c in hub6x6_cb hub12x12_cb hub12x12_cb_b10 hub6x6_cb_flux hub6x6_dense; do r st3pp_$c nv_def_pp $c $E ALF_GPU_SELFTEST=3; r st3log_$c nv_log $c $E ALF_GPU_SELFTEST=3; done
r eps11pp nv_def_pp hub6x6_cb $E ALF_GPU_SELFTEST=3 ALF_GPU_FAULT_STAB=eps:1e-11
r eps11log nv_log hub6x6_cb $E ALF_GPU_SELFTEST=3 ALF_GPU_FAULT_STAB=eps:1e-11
for f in perm phase nan eps:1e-9; do r fpp_$f nv_def_pp hub6x6_cb $E ALF_GPU_SELFTEST=3 ALF_GPU_FAULT_STAB=$f; r flog_$f nv_log hub6x6_cb $E ALF_GPU_SELFTEST=3 ALF_GPU_FAULT_STAB=$f; done
r fhop nv_def hub6x6_cb $E ALF_GPU_SELFTEST=2 ALF_GPU_FAULT=hop
r fupd nv_def hub6x6_cb $E ALF_GPU_UPDTEST=1 ALF_GPU_FAULT=upd
r invstab nv_def_pp hub6x6_cb $E ALF_GPU_SELFTEST=3 ALF_GPU_FAULT_STAB=xyz
r invflt nv_def hub6x6_cb $E ALF_GPU_SELFTEST=2 ALF_GPU_FAULT=xyz
for kv in ALF_GPU_FAULT=upd ALF_GPU_FAULT_STAB=perm ALF_GPU_NS=0 ALF_GPU_STABMIX=1 ALF_GPU_STABCHECK=1 ALF_GPU_TL=8 ALF_GPU_TR=8 ALF_GPU_VB=4 ALF_GPU_NOUPDATE=1 ALF_GPU_TRACE=trace.txt ALF_GPU_TIMING=1 ALF_GPU_STABTIME=1 ALF_GPU_EXIT=1 ALF_GPU_TAUSWAP=8 ALF_GPU_NSS=0; do
  r ref_${kv%%=*} nv_def hub6x6_cb_fd $kv; done
r acc_ns nv_def hub6x6_cb_fd ALF_GPU_NS=0 ALF_GPU_DEBUG=1 ALF_GPU_EXIT=1
r acc_tl nv_def hub6x6_cb_fd ALF_GPU_TL=8 ALF_GPU_DEBUG=1 ALF_GPU_EXIT=1
echo done >> $out
