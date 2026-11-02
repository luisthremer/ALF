#!/bin/bash
# Ticket 207-4 (development only, dev seeds), nividia: QREF run sequence at beta = 10 (hub12x12_cb_b10_p: NSweep 20 x NBin 8 = 160 passes,
# NSTM 10, N_FL 2), the method of t207_3_chain.sh (beta = 4).  GPU runs one at a time (chain holds ~/alf-gpu/gpu.lock through run + processing;
# t207_4_run.sh: no other compute process, disk check inside the lock), each followed by its CPU-only processing with the approved pipeline g5/t207_d5_proc.sh (gpu-port-207-2
# 4118a1de, unchanged; D5W = t207-4: spec, manifest, binary128 program, strict check, then thinning to every KEEP-th dump) before the next.
#   Decision 35: stride 3; C = NSTM*N_FL = 20, P = C/gcd(C, 3) = 20 positions.  FROM 126: passes 126..160 = 35 passes, 700 calls,
#   234 dumps, 210 with B (as 302c-2).  Keep rule proposed by 207-3: all P positions are kept if gcd(m, P) = 1 and m <= (ndumps-1)/(P-1)
#   = 233/19 = 12.3, m = cover*keep -> KEEP 11 (gcd(11, 20) = 1): dumps 1, 12, ..., 232 = 22 kept, positions 11 j mod 20, j = 0..21, all 20.
#   QS 20 (quad self-check every 20th dump), HOSTIN 0.
#   usage: t207_4_chain.sh PHASE [ARGS]  -> t207-4/logs/chain_<PHASE>.txt
#   main SETS...      main runs m_{def,log}_<set>: QVAR=8 (stand-in control) and QHOST=1
#   ctl  SET          production controls c_def_<set> (ALF_GPU_TAU_NS=0) and c_log_<set> (ALF_GPU_TAUSWAP=8); no QVAR
#   flt  PATH SET X.. graded sweep output fault f_<path>_<tag(X)>_<set> (ALF_GPU_FAULT_TAU=eps:X, no QVAR)
#   fltat PATH SET NST:NFE X..  localized graded fault (ALF_GPU_FAULT_TAU=eps:X ALF_GPU_FAULT_TAU_AT=NST:NFE)
W=$HOME/alf-gpu/work/t207-4; cd $W || exit 1; export TMPDIR=$W/tmp
P=$1; shift; out=$W/logs/chain_$P.txt
BD=${BD:-nv_def_pp}; BL=${BL:-nv_log}; CFG=hub12x12_cb_b10_p
NUL=0,4,100004,1,100001,2,1608,4808,5,100005; NULL=0,4,100004,1,100001,2,5,100005   # LOG: no compact-WY codes (flag 8)
FROM=126; STRIDE=3; KEEP=11; QS=20
PROC=$W/g5/t207_d5_proc.sh
NEED=${DISK_NEED_KB:-6500000}      # 5 GB floor (machine rule) + one unthinned beta = 10 run: measured 234 x 5976632 B = 1.40 GB of dumps + ~6 MB
                                   # outputs, rounded up to 1.5 GB (was 6900000 until 2026-11-02 06:50: it held the lock waiting at 6.63 GB free)
# Every run and its processing hold gpu.lock together (pattern of 207-3c b2c71200, HOLDLOCK): fd 9 of this shell, t207_4_run.sh gets
# T207_4_LOCKED=1, so no other agent's GPU run (and its dumps) starts while this run's dumps are on disk; released after the thinning.
r() { local tag=$1 set=$2 bin=$3; shift 3
      echo "df / $(df -h / | tail -1 | awk '{print $4}') free before $tag (lock + disk check in t207_4_run.sh)" >> $out
      SET=$set TAG=$tag DISK_NEED_KB=$NEED bash $W/t207_4_run.sh $bin $CFG "$@" 9>&- >> $out 2>&1; }
rp() { local tag=$1 path=$2; shift 2
       exec 9> "$HOME/alf-gpu/gpu.lock"; flock 9; echo "gpu.lock taken for $tag + processing $(date -Is)" >> $out
       export T207_4_LOCKED=1; r $tag "$@"; pr $tag $path; unset T207_4_LOCKED
       flock -u 9; exec 9>&-; echo "gpu.lock released after $tag $(date -Is)" >> $out; }
pr() { local nl=$NUL; [ "$2" = log ] && nl=$NULL
       [ -f $W/runs/$1/$CFG/provenance.json ] || { echo "no run $1: not processed" >> $out; return 1; }
       grep -q '"rc": 0,' $W/runs/$1/$CFG/provenance.json || { echo "run $1 rc != 0, not processed" >> $out; return 1; }
       D5W=$W D5RUNS=$W/runs bash $PROC $W/runs/$1/$CFG $2 $QS $nl 0 $KEEP >> $out 2>&1
       echo "proc exit $? for $1; df / $(df -h / | tail -1 | awk '{print $4}') free; $(du -sh $W/runs/$1 | cut -f1) in $1" >> $out; }
Q="ALF_GPU_DEBUG=1 ALF_GPU_TAUDIAG_QHOST=1 ALF_GPU_TAUDIAG_QREF=$STRIDE ALF_GPU_TAUDIAG_FROM=$FROM"
xtag() { echo "$1" | tr '.-' 'pm'; }
echo "start $P $* $(date -Is)" >> $out
case $P in
main)
  for s in "$@"; do
    rp m_def_$s def $s $BD $Q ALF_GPU_TAUDIAG_QVAR=8
    rp m_log_$s log $s $BL $Q ALF_GPU_TAUDIAG_QVAR=8
  done ;;
ctl)
  s=$1
  rp c_def_$s def $s $BD $Q ALF_GPU_TAU_NS=0
  rp c_log_$s log $s $BL $Q ALF_GPU_TAUSWAP=8 ;;
flt)
  p=$1; s=$2; shift 2; b=$BD; [ "$p" = log ] && b=$BL
  for x in "$@"; do t=f_${p}_$(xtag $x)_$s; rp $t $p $s $b $Q ALF_GPU_FAULT_TAU=eps:$x; done ;;
fltat)
  p=$1; s=$2; pos=$3; shift 3; b=$BD; [ "$p" = log ] && b=$BL
  for x in "$@"; do t=f_${p}_$(xtag $x)_at$(echo $pos | tr ':' 'x')_$s; rp $t $p $s $b $Q ALF_GPU_FAULT_TAU=eps:$x ALF_GPU_FAULT_TAU_AT=$pos; done ;;
*) echo "bad phase"; exit 1 ;;
esac
echo "done $P $(date -Is)" >> $out
