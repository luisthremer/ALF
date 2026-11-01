#!/bin/bash
# Ticket 207 dev4 (rounding-order nulls of the paired test, development only), nividia: the run sequence (one GPU process at a time,
# each under flock ~/alf-gpu/gpu.lock via t207_d4_run.sh).   usage: t207_d4_chain.sh smoke|ident|main   -> dev4/logs/chain_<phase>.txt
#   smoke: short LOG run with all chains (FROM=1);
#   ident: production-like runs (no debug variable) of HEAD (base_*) and the working tree, data.h5 compared bit for bit;
#   main:  hub12x12_cb_p on dev1..dev4, default (nv_def_pp) and LOG (nv_log), the gate's TAUDIAG=2 / FROM=331 plus V (0, 8), four
#          amplitude-1 chains and the null chains.
W=$HOME/alf-gpu/work/t207/dev4; cd $W || exit 1
P=$1; out=$W/logs/chain_$P.txt
r() { local tag=$1 set=$2 bin=$3 cfg=$4; shift 4; SET=$set TAG=$tag bash $W/t207_d4_run.sh $bin $cfg "$@" >> $out 2>&1; }
D="ALF_GPU_DEBUG=1 ALF_GPU_TAUDIAG=2 ALF_GPU_TAUDIAG_VAR=0,8 ALF_GPU_TAUDIAG_AMP=1,1,1,1"
ND="ALF_GPU_TAUDIAG_NUL=0,16,100004,200004,100001,200001,2,1608,4808,100005,200005"
NL="ALF_GPU_TAUDIAG_NUL=0,16,100004,200004,100001,200001,2,100005,200005"
case $P in
smoke)
  r s1_log smoke1 nv_log hub12x12_cb_p $D ALF_GPU_TAUDIAG_FROM=1 $NL ;;
ident)
  for b in nv_def_pp base_nv_def_pp nv_log base_nv_log; do r id_$b dev1 $b hub12x12_cb_b10_s; done
  for p in def_pp log; do
    echo "h5cmp $p: $(bash $W/t207_d4_h5cmp.sh runs/id_nv_$p/hub12x12_cb_b10_s/data.h5 runs/id_base_nv_$p/hub12x12_cb_b10_s/data.h5) rc=$?" >> $out
  done ;;
main)
  for s in dev1 dev2 dev3 dev4; do
    r m_def_$s $s nv_def_pp hub12x12_cb_p $D ALF_GPU_TAUDIAG_FROM=331 $ND
    r m_log_$s $s nv_log hub12x12_cb_p $D ALF_GPU_TAUDIAG_FROM=331 $NL
  done ;;
*) echo "bad phase"; exit 1 ;;
esac
echo "done $P $(date -Is)" >> $out
