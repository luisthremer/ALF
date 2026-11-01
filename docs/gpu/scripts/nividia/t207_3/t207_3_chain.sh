#!/bin/bash
# Ticket 207-3 (development only, dev seeds), nividia: QREF run sequence at beta = 4 (hub12x12_cb_p, 400 passes).  GPU runs one at a time
# (t207_3_run.sh: flock ~/alf-gpu/gpu.lock and no other compute process), each followed by its CPU-only processing with the unchanged
# pipeline (t207_d5_proc.sh with D5W = t207-3: spec, manifest, quad program, strict check, then thinning to every KEEP-th dump).
#   Decision 35: stride 3 (cycles all 8 (nst, flavour) positions, C = NSTM*N_FL = 8), FROM = 296: passes 296..400 = 105 passes, 840 calls,
#   280 dumps, 210 with B (nst < 4).  Keep rule: P = 8, gcd(39, 8) = 1 -> KEEP 39.  QS 20 (quad self-check every 20th dump).
#   usage: t207_3_chain.sh PHASE [ARGS]  -> t207-3/logs/chain_<PHASE>.txt
#   main SETS...      main runs m_{def,log}_<set>: QVAR=8 (stand-in control) and QHOST=1
#   ctl  SET          production controls c_def_<set> (ALF_GPU_TAU_NS=0) and c_log_<set> (ALF_GPU_TAUSWAP=8); one env form: no QVAR
#   flt  PATH SET X.. graded output fault f_<path>_<tag(X)>_<set> (ALF_GPU_FAULT_TAU=eps:X, no QVAR)
#   fltat PATH SET NST:NFE X..  localized graded fault (ALF_GPU_FAULT_TAU=eps:X ALF_GPU_FAULT_TAU_AT=NST:NFE, needs BD=nv_def_pp2 / BL=nv_log2)
#   attr SET          default path, TAUDIAG=2 (needed for host-twin masks) + QVAR=8,16,32,1,128,2,4,256 (task 3 attribution)
W=${T207_3W:-$HOME/alf-gpu/work/t207-3};   # T207_3W: work dir of a follow-up agent (t207-3b)
cd $W || exit 1
P=$1; shift; out=$W/logs/chain_$P.txt
BD=${BD:-nv_def_pp}; BL=${BL:-nv_log}
NUL=0,4,100004,1,100001,2,1608,4808,5,100005; NULL=0,4,100004,1,100001,2,5,100005   # LOG: no compact-WY codes (flag 8)
FROM=296; STRIDE=3; KEEP=39; QS=20
# processing: the redesigned pipeline of gpu-port-207-2 305488ec in $W/g3 (records qref_sha256 and prog_src, needed by its gate); the
# runs m_{def,log}_dev1 were processed with the 207-2 fix-round pipeline (a0912a29) in $W (PROC_DIR=$W to repeat that)
PROC=${PROC_DIR:-$W/g3}/t207_d5_proc.sh
# disk guard (follow-up t207-3b): a run writes ~1.6 GB of dumps before its processing thins them; it starts only if / has at least
# MINFREE_MB free (default 6700 = 5 GB floor + 1.6 GB), else it is skipped (logged) and so is its processing
r() { local tag=$1 set=$2 bin=$3 cfg=$4; shift 4
      local fm=$(df -B1M --output=avail / | tail -1 | tr -d ' ')
      if [ "$fm" -lt "${MINFREE_MB:-6700}" ]; then echo "SKIPPED $tag: / has $fm MB free < ${MINFREE_MB:-6700} MB" >> $out; return 1; fi
      echo "df / $fm MB free before $tag" >> $out
      SET=$set TAG=$tag bash $W/t207_3_run.sh $bin $cfg "$@" >> $out 2>&1; }
pr() { [ -f $W/runs/$1/hub12x12_cb_p/provenance.json ] || { echo "no run $1, not processed" >> $out; return 1; }
       local nl=$NUL; [ "$2" = log ] && nl=$NULL; D5W=$W D5RUNS=$W/runs bash $PROC $W/runs/$1/hub12x12_cb_p $2 $QS $nl 0 $KEEP >> $out 2>&1
       echo "df / $(df -h / | tail -1 | awk '{print $4}') free after $1" >> $out; }
Q="ALF_GPU_DEBUG=1 ALF_GPU_TAUDIAG_QHOST=1 ALF_GPU_TAUDIAG_QREF=$STRIDE ALF_GPU_TAUDIAG_FROM=$FROM"
xtag() { echo "$1" | tr '.-' 'pm'; }
echo "start $P $* $(date -Is)" >> $out
case $P in
main)
  for s in "$@"; do
    r m_def_$s $s $BD hub12x12_cb_p $Q ALF_GPU_TAUDIAG_QVAR=8; pr m_def_$s def
    r m_log_$s $s $BL hub12x12_cb_p $Q ALF_GPU_TAUDIAG_QVAR=8; pr m_log_$s log
  done ;;
main_def)   # default path only (end-to-end check of the approved gate 4118a1de with PROC_DIR=$W/g5)
  for s in "$@"; do r m_def_$s $s $BD hub12x12_cb_p $Q ALF_GPU_TAUDIAG_QVAR=8; pr m_def_$s def; done ;;
ctl)
  s=$1
  r c_def_$s $s $BD hub12x12_cb_p $Q ALF_GPU_TAU_NS=0; pr c_def_$s def
  r c_log_$s $s $BL hub12x12_cb_p $Q ALF_GPU_TAUSWAP=8; pr c_log_$s log ;;
flt)
  p=$1; s=$2; shift 2; b=$BD; [ "$p" = log ] && b=$BL
  for x in "$@"; do t=f_${p}_$(xtag $x)_$s; r $t $s $b hub12x12_cb_p $Q ALF_GPU_FAULT_TAU=eps:$x; pr $t $p; done ;;
fltat)
  p=$1; s=$2; pos=$3; shift 3; b=$BD; [ "$p" = log ] && b=$BL
  for x in "$@"; do t=f_${p}_$(xtag $x)_at$(echo $pos | tr ':' 'x')_$s; r $t $s $b hub12x12_cb_p $Q ALF_GPU_FAULT_TAU=eps:$x ALF_GPU_FAULT_TAU_AT=$pos; pr $t $p; done ;;
attr)
  s=$1
  r a_def_$s $s $BD hub12x12_cb_p $Q ALF_GPU_TAUDIAG=2 ALF_GPU_TAUDIAG_QVAR=8,16,32,1,128,2,4,256; pr a_def_$s def ;;
*) echo "bad phase"; exit 1 ;;
esac
echo "done $P $(date -Is)" >> $out
