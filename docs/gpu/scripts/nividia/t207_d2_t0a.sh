#!/bin/bash
# Ticket 207 (rebase on 212, development only), nividia: T0 self-tests (SELFTEST=4, 6 configs x default/LOG), fault controls incl.
# the graded fault, production refusal of every 207 variable, invalid fault values.  -> dev2/logs/t0a_<RUN>.txt, runs in
# dev2/runs/t0a<RUN>_<id>/ .   usage: t207_d2_t0a.sh RUN [selftest|faults|refusal ...]   (RUN: a fresh label, e.g. 1)
W=$HOME/alf-gpu/work/t207/dev2; R=$1; shift; out=$W/logs/t0a_$R.txt
case "$R" in ""|*/*|*..*) echo bad run label; exit 1;; esac
r() { # id binary config env...
  id=$1; b=$2; c=$3; shift 3
  res=$(TAG=t0a${R}_$id bash $W/t207_d2_run.sh $b $c "$@")
  echo "== $id: $res" >> $out
  grep -hE "PASSED|FAILED|GPU debug env|GPU fault|GPU build|refused|invalid|tau self-test|GPU tau pass|branch" $W/runs/t0a${R}_$id/$c/run.log | cut -c1-260 | uniq >> $out
}
E="ALF_GPU_EXIT=1"
for part in "$@"; do case $part in
selftest)
  for c in hub6x6_cb hub6x6_dense hub12x12_cb hub12x12_cb_b10 hub6x6_cb_flux hub6x6_cb_symm; do
    r st4pp_$c nv_def_pp $c $E ALF_GPU_SELFTEST=4; r st4log_$c nv_log $c $E ALF_GPU_SELFTEST=4; done ;;
faults)
  for f in perm prop prec blocks eps:1e-9 eps:1e-11; do
    r fpp_${f/:/_} nv_def_pp hub6x6_cb $E ALF_GPU_SELFTEST=4 ALF_GPU_FAULT_TAU=$f
    r flog_${f/:/_} nv_log hub6x6_cb $E ALF_GPU_SELFTEST=4 ALF_GPU_FAULT_TAU=$f; done
  for f in xyz eps:abc eps:-1 eps:; do r inv_${f/:/_} nv_def_pp hub6x6_cb $E ALF_GPU_SELFTEST=4 ALF_GPU_FAULT_TAU=$f; done ;;
refusal)
  for kv in ALF_GPU_TAU_NS=0 ALF_GPU_TAUSWAP=8 ALF_GPU_TAUDIAG=2 ALF_GPU_TAUDIAG_FROM=331 ALF_GPU_TAUTIME=1 ALF_GPU_TAUTEST_V=1 \
            ALF_GPU_FAULT_TAU=perm ALF_GPU_FAULT_TAU=eps:1e-9 ALF_GPU_FAULT_TAU=xyz ALF_GPU_TAUDIAG_SEED=1; do
    r ref_$(echo ${kv} | tr "=:" "__") nv_def hub6x6_cb_fd $kv; done ;;
esac; done
echo done >> $out
