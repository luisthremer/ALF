#!/bin/bash
# Ticket 207 dev5 (accuracy test against a quad reference, development only), nividia: run sequence.  GPU runs one at a time
# (t207_d5_run.sh: flock ~/alf-gpu/gpu.lock, nvidia-smi empty), each followed by its CPU-only offline processing
# (t207_d5_proc.sh, 4 threads; the GPU is idle then).   usage: t207_d5_chain.sh PHASE [SETS]  -> dev5/logs/chain_<PHASE>.txt
#   smoke: hub12x12_cb_p, last 2 passes (FROM=399), every point, default and LOG
#   ident: production-like runs (no debug variable) of base (c1cc528e) and the working tree, data.h5 compared bit for bit
#   main:  hub12x12_cb_p on SETS (default dev1..dev4), default (nv_def_pp) and LOG (nv_log), FROM=331, every 2nd point,
#          device variant 8 (= the planted defect: device ZUNMQR without the Newton-Schulz step) and the nvfortran host CGR2_2
#   ctl:   production control runs on dev1 (default ALF_GPU_TAU_NS=0, LOG ALF_GPU_TAUSWAP=8), FROM=331, every 2nd point
#   din:   TAUDIAG=2 runs on dev1 with the host chain's state in the dumps (206 input part), every 4th point
W=$HOME/alf-gpu/work/t207/dev5; cd $W || exit 1
P=$1; shift; SETS=${*:-dev1 dev2 dev3 dev4}; out=$W/logs/chain_$P.txt
BD=${BD:-nv_def_pp2}; BL=${BL:-nv_log2}     # GPU binaries (default / LOG); smoke used nv_def_pp / nv_log (dump version 1)
NUL=0,4,100004,1,100001,2,1608,4808,5,100005; NULL=0,4,100004,1,100001,2,5,100005   # LOG: no compact-WY codes (flag 8)
r() { local tag=$1 set=$2 bin=$3 cfg=$4; shift 4; SET=$set TAG=$tag bash $W/t207_d5_run.sh $bin $cfg "$@" >> $out 2>&1; }
pr() { local nl=$NUL; [ "$2" = log ] && nl=$NULL; bash $W/t207_d5_proc.sh $W/runs/$1/hub12x12_cb_p $2 $3 $nl $4 $5 >> $out 2>&1; }
Q="ALF_GPU_DEBUG=1 ALF_GPU_TAUDIAG_QHOST=1"
case $P in
smoke)
  r s_def dev1 $BD hub12x12_cb_p $Q ALF_GPU_TAUDIAG_QREF=1 ALF_GPU_TAUDIAG_FROM=399 ALF_GPU_TAUDIAG_QVAR=8
  pr s_def def 1 0 0
  r s_log dev1 $BL hub12x12_cb_p $Q ALF_GPU_TAUDIAG_QREF=1 ALF_GPU_TAUDIAG_FROM=399 ALF_GPU_TAUDIAG_QVAR=8
  pr s_log log 1 0 0 ;;
ident)
  for b in $BD base_nv_def_pp $BL base_nv_log; do r id_$b dev1 $b hub12x12_cb_b10_s; done
  echo "h5cmp def: $(bash $W/t207_d5_h5cmp.sh runs/id_$BD/hub12x12_cb_b10_s/data.h5 runs/id_base_nv_def_pp/hub12x12_cb_b10_s/data.h5) rc=$?" >> $out
  echo "h5cmp log: $(bash $W/t207_d5_h5cmp.sh runs/id_$BL/hub12x12_cb_b10_s/data.h5 runs/id_base_nv_log/hub12x12_cb_b10_s/data.h5) rc=$?" >> $out ;;
main)
  for s in $SETS; do
    r m_def_$s $s $BD hub12x12_cb_p $Q ALF_GPU_TAUDIAG_QREF=2 ALF_GPU_TAUDIAG_FROM=331 ALF_GPU_TAUDIAG_QVAR=8
    pr m_def_$s def 20 0 40
    r m_log_$s $s $BL hub12x12_cb_p $Q ALF_GPU_TAUDIAG_QREF=2 ALF_GPU_TAUDIAG_FROM=331 ALF_GPU_TAUDIAG_QVAR=8
    pr m_log_$s log 20 0 40
  done ;;
ctl)
  r c_def_dev1 dev1 $BD hub12x12_cb_p $Q ALF_GPU_TAUDIAG_QREF=2 ALF_GPU_TAUDIAG_FROM=331 ALF_GPU_TAU_NS=0
  pr c_def_dev1 def 0 0 40
  r c_log_dev1 dev1 $BL hub12x12_cb_p $Q ALF_GPU_TAUDIAG_QREF=2 ALF_GPU_TAUDIAG_FROM=331 ALF_GPU_TAUSWAP=8
  pr c_log_dev1 log 0 0 40 ;;
din)
  r d_def_dev1 dev1 $BD hub12x12_cb_p $Q ALF_GPU_TAUDIAG=2 ALF_GPU_TAUDIAG_QREF=4 ALF_GPU_TAUDIAG_FROM=331
  pr d_def_dev1 def 0 1 10
  r d_log_dev1 dev1 $BL hub12x12_cb_p $Q ALF_GPU_TAUDIAG=2 ALF_GPU_TAUDIAG_QREF=4 ALF_GPU_TAUDIAG_FROM=331
  pr d_log_dev1 log 0 1 10 ;;
*) echo "bad phase"; exit 1 ;;
esac
echo "done $P $(date -Is)" >> $out
