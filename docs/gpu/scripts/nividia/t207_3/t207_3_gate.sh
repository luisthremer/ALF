#!/bin/bash
# Ticket 207-3 (development only), nividia: run the unmodified redesigned gate (GD=g3: copy of gpu-port-207-2 305488ec, the version that
# processed the runs; GD=g4: ffab796f, GD=g5: 4118a1de (approved); g4/g5 need runs processed by them) on genuine
# stride-3 cover-1 run sets with a spec written by t207_3_spec.py (bounds are calibration candidates, not a pre-registration).
#   usage: t207_3_gate.sh TAG PATH "MAINSETS" CTLSET|- [extra t207_3_spec.py args...]   -> gate/{spec,report,gate}_<TAG>_<PATH>.*
W=${T207_3W:-$HOME/alf-gpu/work/t207-3};   # T207_3W: work dir of a follow-up agent (t207-3b)
cd $W || exit 1
TAG=$1; P=$2; MAINS=$3; CTL=$4; shift 4
PY=~/alf-gpu/venv/bin/python3
GD=${GD:-g3}
q=qref_pp; [ "$P" = log ] && q=qref_log
PS=$(cut -d' ' -f1 bin/$q.out.sha256); SS=$(sha256sum $GD/t207_d5_qref.F90 | cut -d' ' -f1)
R=(); D=(); n=0
for s in $MAINS; do R+=(--rundir main=runs/m_${P}_$s/hub12x12_cb_p); D+=(runs/m_${P}_$s/hub12x12_cb_p); n=$((n+1)); done
C=()
if [ "$CTL" != - ]; then R+=(--rundir ctl=runs/c_${P}_$CTL/hub12x12_cb_p); D+=(runs/c_${P}_$CTL/hub12x12_cb_p); C=(--ctl 1 --control ctl:run:ctl:acc_pTT); fi
S=gate/spec_${TAG}_$P.json; O=gate/report_${TAG}_$P.json; T=gate/gate_${TAG}_$P.txt
[ -e $S ] || [ -e $O ] && { echo "refusing: $S or $O exists"; exit 1; }
if [ $GD != g3 ]; then
  $PY -I t207_3_spec.py $S $P $n "${C[@]}" --control var8:variant:var8:acc_pTT --prog-sha $PS --src-sha $SS "${R[@]}" "$@" || exit 1
  $PY -I $GD/t207_d5_gate.py --spec $S --spec-sha $(sha256sum $S | cut -d' ' -f1) --prog bin/$q.out --report $O "${D[@]}" > $T 2>&1; rc=$?
else
  $PY -I t207_3_spec.py $S $P $n "${C[@]}" --control var8:variant:var8:acc_pTT --prog-sha $PS --src-sha $SS "$@" || exit 1
  $PY -I g3/t207_d5_gate.py --spec $S --report $O "${D[@]}" > $T 2>&1; rc=$?
fi
echo "gate $TAG $P rc=$rc"; grep -E "^  gate |^  control |^GATE" $T | cut -c1-300
