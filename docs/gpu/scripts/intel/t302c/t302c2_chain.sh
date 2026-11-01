#!/bin/bash
# Ticket 302c-2 (development only, dev seeds), on intel: QREF run sequence on the B580, the method of 207-3 (t207_3_chain.sh) on nividia.
# Every GPU run (t302c2_run.sh: load check, flock ~/alf-gpu/gpu.lock) is followed by its CPU-only processing with the unchanged 207-2
# pipeline (t207_d5_proc.sh of commit 5d876291 with D5W = t302c-2: spec, manifest, binary128 program, strict check, then thinning to every
# KEEP-th dump) before the next run starts.
#   Decision 35: stride 3.  beta = 4 (hub12x12_cb_p, 400 passes, NSTM 4, N_FL 2, C = 8): FROM 296 -> passes 296..400, 840 calls, 280 dumps,
#   210 with B.  beta = 10 (hub12x12_cb_b10_p, 160 passes, NSTM 10, C = 20): FROM 126 -> passes 126..160, 700 calls, 234 dumps, 210 with B.
#   Keep rule: P = 8 resp. 20, gcd(39, P) = 1 -> KEEP 39.  QS 20 (quad self-check every 20th dump).  HOSTIN 0.
#   Null codes (offline, gfortran + OpenBLAS oracle objects): default 0,4,100004,1,100001,2,1608,4808,5,100005; LOG without the compact-WY
#   codes (flag 8 is refused under STABLOG): 0,4,100004,1,100001,2,5,100005.
#   usage: t302c2_chain.sh PHASE [ARGS]  -> t302c-2/logs/chain_<PHASE>.txt
#     smoke            hub12x12_cb_p dev1, FROM 399 (passes 399..400: 6 dumps), main env, default and LOG, KEEP 0
#     main  B SETS...  main runs m{B}_{def,log}_<set> (B = 4 | 10): QHOST=1 (the in-run ifx/MKL host CGR2_2), QVAR=8 (stand-in control)
#     ctl   B SET      production controls c{B}_def_<set> (ALF_GPU_TAU_NS=0) and c{B}_log_<set> (ALF_GPU_TAUSWAP=8); no QVAR (207-3 form)
W=$HOME/alf-gpu/work/t302c-2; cd $W || exit 1; export TMPDIR=$W/tmp
P=$1; shift; out=$W/logs/chain_$P.txt
BD=${BD:-in_def_pp}; BL=${BL:-in_log}
NUL=0,4,100004,1,100001,2,1608,4808,5,100005; NULL=0,4,100004,1,100001,2,5,100005
STRIDE=3; KEEP=39; QS=20
cfg() { case $1 in 4) echo hub12x12_cb_p;; 10) echo hub12x12_cb_b10_p;; *) echo bad;; esac; }
frm() { case $1 in 4) echo 296;; 10) echo 126;; *) echo bad;; esac; }
r() { local tag=$1 set=$2 bin=$3 cfg=$4; shift 4; SET=$set TAG=$tag bash $W/t302c2_run.sh $bin $cfg "$@" >> $out 2>&1; }
pr() { local nl=$NUL; [ "$2" = log ] && nl=$NULL
       D5W=$W D5RUNS=$W/runs bash $W/t207_d5_proc.sh $W/runs/$1/$3 $2 $QS $nl 0 ${4:-$KEEP} >> $out 2>&1
       echo "proc exit $? for $1; df / $(df -h / | tail -1 | awk '{print $4}') free; $(du -sh $W/runs/$1 | cut -f1) in $1" >> $out; }
echo "start $P $* $(date -Is)" >> $out
case $P in
smoke)
  Q="ALF_GPU_DEBUG=1 ALF_GPU_TAUDIAG_QHOST=1 ALF_GPU_TAUDIAG_QREF=$STRIDE ALF_GPU_TAUDIAG_FROM=399 ALF_GPU_TAUDIAG_QVAR=8"
  r s_def dev1 $BD hub12x12_cb_p $Q; pr s_def def hub12x12_cb_p 0
  r s_log dev1 $BL hub12x12_cb_p $Q; pr s_log log hub12x12_cb_p 0 ;;
main)
  b=$1; shift; c=$(cfg $b); [ "$c" = bad ] && { echo "bad beta" >> $out; exit 1; }
  Q="ALF_GPU_DEBUG=1 ALF_GPU_TAUDIAG_QHOST=1 ALF_GPU_TAUDIAG_QREF=$STRIDE ALF_GPU_TAUDIAG_FROM=$(frm $b)"
  for s in "$@"; do
    r m${b}_def_$s $s $BD $c $Q ALF_GPU_TAUDIAG_QVAR=8; pr m${b}_def_$s def $c
    r m${b}_log_$s $s $BL $c $Q ALF_GPU_TAUDIAG_QVAR=8; pr m${b}_log_$s log $c
  done ;;
ctl)
  b=$1; s=$2; c=$(cfg $b); [ "$c" = bad ] && { echo "bad beta" >> $out; exit 1; }
  Q="ALF_GPU_DEBUG=1 ALF_GPU_TAUDIAG_QHOST=1 ALF_GPU_TAUDIAG_QREF=$STRIDE ALF_GPU_TAUDIAG_FROM=$(frm $b)"
  r c${b}_def_$s $s $BD $c $Q ALF_GPU_TAU_NS=0; pr c${b}_def_$s def $c
  r c${b}_log_$s $s $BL $c $Q ALF_GPU_TAUSWAP=8; pr c${b}_log_$s log $c ;;
*) echo "bad phase" >> $out; exit 1 ;;
esac
echo "done $P $(date -Is)" >> $out
