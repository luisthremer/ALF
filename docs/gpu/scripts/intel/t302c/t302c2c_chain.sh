#!/bin/bash
# Ticket 302c-2 round 3 (development only, dev seeds), on intel: QREF runs on the B580 at HEAD d3e82a36 (final 207-3 merged), the method of
# round 2 (t302c2b_chain.sh, unchanged settings): every GPU run (t302c2_run.sh: load check, flock ~/alf-gpu/gpu.lock) is followed by its
# CPU-only processing with the approved 207-2 pipeline (t207_d5_proc.sh of 4118a1de) before the next run starts.
#   Decision 35: stride 3.  beta = 4 (hub12x12_cb_p): FROM 296, 280 dumps, 210 with B, P = 8.  beta = 10 (hub12x12_cb_b10_p): FROM 126,
#   234 dumps, 210 with B, P = 20.  QS 20, HOSTIN 0.  Main runs follow the keep rule (KEEP 39 at beta 4, KEEP 11 at beta 10); fault runs keep
#   KEEP_F = 117 (gcd 1 with P = 8 and 20) as in round 2 (live-evidence budget; the statistics use the full qref.txt).
#   Null codes: default 0,4,100004,1,100001,2,1608,4808,5,100005; LOG without the compact-WY codes 1608, 4808.
#   Before every run: >= 5 GB free on / and the work dir below 1.4 GB (live-evidence budget 1.5 GB), else the chain stops.
#   usage: [PATHS="def log"] t302c2c_chain.sh PHASE ARGS  -> t302c-2c/logs/chain_<PHASE>.txt
#     main  B SETS...               m{B}_<path>_<set> for every path in PATHS: QHOST=1, QVAR=8 (207-3 main env)
#     flt   B PATH SET X...         f{B}_<path>_<x>_<set>: ALF_GPU_FAULT_TAU=eps:X (no QVAR)
#     fltat B PATH SET NST:NFE X... f{B}_<path>_<x>_at<nst>x<nfe>_<set>: plus ALF_GPU_FAULT_TAU_AT=NST:NFE (NST = 0: flavour-only fault)
W=$HOME/alf-gpu/work/t302c-2c; cd $W || exit 1; export TMPDIR=$W/tmp T302C2_W=$W
P=$1; shift; out=$W/logs/chain_$P.txt
BD=${BD:-in_def_pp}; BL=${BL:-in_log}; PATHS=${PATHS:-def log}
NUL=0,4,100004,1,100001,2,1608,4808,5,100005; NULL=0,4,100004,1,100001,2,5,100005
STRIDE=3; QS=20; KEEP_F=${KEEP_F:-117}
cfg() { case $1 in 4) echo hub12x12_cb_p;; 10) echo hub12x12_cb_b10_p;; *) echo bad;; esac; }
frm() { case $1 in 4) echo 296;; 10) echo 126;; *) echo bad;; esac; }
kmain() { case $1 in 4) echo 39;; 10) echo 11;; esac; }
space() { local f u; f=$(df -k --output=avail / | tail -1); u=$(du -sk $W | cut -f1)
          [ "$f" -ge 5242880 ] && [ "$u" -le 1468006 ] || { echo "STOP: disk $f kB free on /, $u kB in t302c-2c (need >= 5 GB free, <= 1.4 GB used)" >> $out; exit 3; }; }
r() { local tag=$1 set=$2 bin=$3 cfg=$4; shift 4; space; SET=$set TAG=$tag bash $W/t302c2_run.sh $bin $cfg "$@" >> $out 2>&1; }
pr() { local nl=$NUL; [ "$2" = log ] && nl=$NULL
       D5W=$W D5RUNS=$W/runs bash $W/t207_d5_proc.sh $W/runs/$1/$3 $2 $QS $nl 0 $4 >> $out 2>&1
       echo "proc exit $? for $1 (KEEP $4); df / $(df -h / | tail -1 | awk '{print $4}') free; $(du -sh $W/runs/$1 | cut -f1) in $1; $(du -sh $W | cut -f1) in t302c-2c" >> $out; }
xtag() { echo "$1" | tr '.-' 'pm'; }
bin() { [ "$1" = log ] && echo $BL || echo $BD; }
echo "start $P $* $(date -Is)" >> $out
case $P in
main)
  b=$1; shift; c=$(cfg $b); [ "$c" = bad ] && { echo "bad beta" >> $out; exit 1; }
  Q="ALF_GPU_DEBUG=1 ALF_GPU_TAUDIAG_QHOST=1 ALF_GPU_TAUDIAG_QREF=$STRIDE ALF_GPU_TAUDIAG_FROM=$(frm $b)"
  for s in "$@"; do
    for p in $PATHS; do [ -d runs/$s/$c ] || { echo "missing runs/$s/$c (prep)" >> $out; exit 1; }
      r m${b}_${p}_$s $s $(bin $p) $c $Q ALF_GPU_TAUDIAG_QVAR=8; pr m${b}_${p}_$s $p $c $(kmain $b); done
  done ;;
flt)
  b=$1; p=$2; s=$3; shift 3; c=$(cfg $b); [ "$c" = bad ] && { echo "bad beta" >> $out; exit 1; }
  Q="ALF_GPU_DEBUG=1 ALF_GPU_TAUDIAG_QHOST=1 ALF_GPU_TAUDIAG_QREF=$STRIDE ALF_GPU_TAUDIAG_FROM=$(frm $b)"
  for x in "$@"; do t=f${b}_${p}_$(xtag $x)_$s; r $t $s $(bin $p) $c $Q ALF_GPU_FAULT_TAU=eps:$x; pr $t $p $c $KEEP_F; done ;;
fltat)
  b=$1; p=$2; s=$3; pos=$4; shift 4; c=$(cfg $b); [ "$c" = bad ] && { echo "bad beta" >> $out; exit 1; }
  Q="ALF_GPU_DEBUG=1 ALF_GPU_TAUDIAG_QHOST=1 ALF_GPU_TAUDIAG_QREF=$STRIDE ALF_GPU_TAUDIAG_FROM=$(frm $b)"
  for x in "$@"; do t=f${b}_${p}_$(xtag $x)_at$(echo $pos | tr ':' 'x')_$s
    r $t $s $(bin $p) $c $Q ALF_GPU_FAULT_TAU=eps:$x ALF_GPU_FAULT_TAU_AT=$pos; pr $t $p $c $KEEP_F; done ;;
*) echo "bad phase" >> $out; exit 1 ;;
esac
echo "done $P $(date -Is)" >> $out
