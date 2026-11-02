#!/bin/bash
# Ticket 207-4 (development only; after t207_3_gate.sh), nividia: run the unmodified approved gate g5/t207_d5_gate.py (gpu-port-207-2
# 4118a1de) on the beta = 10 stride-3 cover-1 run sets with a spec written by t207_4_spec.py (bounds are calibration inputs, not a
# pre-registration).  The gate re-runs the offline program on every kept dump (KEEP 11: 22 per run) and must reproduce qref.txt exactly.
#   usage: t207_4_gate.sh TAG PATH "MAINSETS" "CONTROLS" [extra t207_4_spec.py args (--gate ..., --control ...)]
#     CONTROLS: space-separated  ctl:SET  (TAU_NS=0 / TAUSWAP=8)  and  flt:X:SET  (FAULT_TAU=eps:X); '-' for none.  The localized
#     fault (ALF_GPU_FAULT_TAU_AT) is not in the gate's FAULTS table (it expects eps:<x>@<nst>:<nfe>), so it cannot enter the gate.
#   -> gate/{spec,report,gate}_<TAG>_<PATH>.*  (refuses to overwrite)
W=$HOME/alf-gpu/work/t207-4; cd $W || exit 1; export TMPDIR=$W/tmp
TAG=$1; P=$2; MAINS=$3; CTLS=$4; shift 4
PY=~/alf-gpu/venv/bin/python3; GD=g5; CFG=hub12x12_cb_b10_p
q=qref_pp; [ "$P" = log ] && q=qref_log
PS=$(cut -d' ' -f1 bin/$q.out.sha256); SS=$(sha256sum $GD/t207_d5_qref.F90 | cut -d' ' -f1)
R=(); D=(); C=(); n=0; nc=0; declare -A NF
for s in $MAINS; do R+=(--rundir main=runs/m_${P}_$s/$CFG); D+=(runs/m_${P}_$s/$CFG); n=$((n+1)); done
[ "$CTLS" = - ] && CTLS=""
for c in $CTLS; do
  IFS=: read -r kind a b <<< "$c"
  case $kind in
  ctl) R+=(--rundir ctl=runs/c_${P}_$a/$CFG); D+=(runs/c_${P}_$a/$CFG); nc=$((nc+1)) ;;
  flt) t=$(echo "$a" | tr '.-' 'pm'); role=flt_$t
       R+=(--rundir $role=runs/f_${P}_${t}_$b/$CFG); D+=(runs/f_${P}_${t}_$b/$CFG); NF[$a]=$(( ${NF[$a]:-0} + 1 )) ;;
  *) echo "bad control $c"; exit 1 ;;
  esac
done
[ $nc -gt 0 ] && C+=(--ctl $nc)
for x in "${!NF[@]}"; do C+=(--flt $x:${NF[$x]}); done
S=gate/spec_${TAG}_$P.json; O=gate/report_${TAG}_$P.json; T=gate/gate_${TAG}_$P.txt
{ [ -e $S ] || [ -e $O ]; } && { echo "refusing: $S or $O exists"; exit 1; }
$PY -I t207_4_spec.py $S $P $n "${C[@]}" --prog-sha $PS --src-sha $SS "${R[@]}" "$@" || exit 1
s0=$(date +%s)
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 $PY -I $GD/t207_d5_gate.py --spec $S --spec-sha $(sha256sum $S | cut -d' ' -f1) --prog bin/$q.out --report $O "${D[@]}" > $T 2>&1; rc=$?
echo "gate $TAG $P rc=$rc wall=$(( $(date +%s) - s0 ))s"; grep -E "^  gate |^  control |^GATE|reprocess|authentic" $T | cut -c1-320
