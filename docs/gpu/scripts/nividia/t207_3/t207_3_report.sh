#!/bin/bash
# Ticket 207-3 (development only), nividia: generated calibration tables from the stamped, processed runs (t207_3_cal.py with the strict
# reader of QP/ (default g3 = gpu-port-207-2 305488ec; g5 = 4118a1de, approved), per path: all main runs (dev1..dev4), the production
# control (dev5), the graded faults (every f_<path>_* run).  Runs are taken from $W/runs and, if PREV is set, also from PREV/runs (read
# only; a follow-up agent's work dir W reads the previous agent's runs in place).  -> logs/report_cal_<path>[_<SUF>].{txt,json} and the
# proposed-gate summary logs/summary_gates[_<SUF>].{txt,json} (t207_3_summary.py).  Runs that do not exist yet are skipped (listed).
#   QP3 (optional): a second reader directory for the runs of its spec version (g5: runs processed by the approved pipeline 4118a1de,
#   e.g. QP=g3 QP3=g5 reads the 305488ec-processed runs with g3 and the 4118a1de-processed runs with g5; t207_3_cal.py T207_3_QPDIR_V3)
#   usage: [PREV=dir] [QP=g3|g5] [QP3=g5] [SUF=suffix] t207_3_report.sh
W=${T207_3W:-$HOME/alf-gpu/work/t207-3};   # T207_3W: work dir of a follow-up agent (t207-3b)
cd $W || exit 1
PY=~/alf-gpu/venv/bin/python3
QP=${QP:-g3}; SUF=${SUF:+_$SUF}; V3=(); [ -n "$QP3" ] && V3=(T207_3_QPDIR_V3=$W/$QP3)
ROOTS=(runs); [ -n "$PREV" ] && ROOTS+=("$PREV/runs")
fr() { local r; for r in "${ROOTS[@]}"; do [ -f $r/$1/hub12x12_cb_p/qref_proc.json ] && { echo $r/$1/hub12x12_cb_p; return; }; done; }
for p in def log; do
  args=()
  for s in dev1 dev2 dev3 dev4; do d=$(fr m_${p}_$s); [ -n "$d" ] && args+=("main:m_${p}_$s=$d") || echo "missing m_${p}_$s"; done
  d=$(fr c_${p}_dev5); [ -n "$d" ] && args+=("ctl:c_${p}_dev5=$d") || echo "missing c_${p}_dev5"
  for r in "${ROOTS[@]}"; do
    for d in $r/f_${p}_*/hub12x12_cb_p; do [ -f $d/qref_proc.json ] && args+=("flt:$(basename $(dirname $d))=$d"); done
  done
  out=logs/report_cal_$p$SUF
  rm -f $out.json
  env "${V3[@]}" T207_3_STRATA=1 T207_3_QPDIR=$W/$QP $PY -I t207_3_cal.py $out.json $p "${args[@]}" > $out.txt 2>&1; echo "report $p rc=$? ($out.txt)"
done
$PY -I t207_3_summary.py logs/summary_gates$SUF.json logs/report_cal_def$SUF.json logs/report_cal_log$SUF.json > logs/summary_gates$SUF.txt 2>&1
echo "summary rc=$? (logs/summary_gates$SUF.txt)"
