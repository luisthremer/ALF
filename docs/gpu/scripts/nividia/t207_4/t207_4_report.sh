#!/bin/bash
# Ticket 207-4 (development only; after t207_3_report.sh), nividia: generated calibration tables at beta = 10 from the stamped, processed
# runs: t207_4_cal.py (= t207_3_cal.py) with the strict reader of g5/ (the approved gpu-port-207-2 4118a1de), per path: the main runs
# (dev11..dev13), the production control (dev14), every graded fault run f_<path>_* (dev15..dev18), then t207_4_sep.py (separation of
# the candidate gates by decision 37's rule) on both paths.  -> logs/report_cal_<path>.{txt,json}, logs/report_sep.txt.
# Runs that do not exist (or are not processed) are skipped and listed.
W=$HOME/alf-gpu/work/t207-4; cd $W || exit 1; export TMPDIR=$W/tmp
PY=~/alf-gpu/venv/bin/python3; CFG=hub12x12_cb_b10_p
J=()
for p in def log; do
  args=()
  for s in dev11 dev12 dev13; do d=runs/m_${p}_$s/$CFG; [ -f $d/qref_proc.json ] && args+=("main:m_${p}_$s=$d") || echo "missing $d"; done
  d=runs/c_${p}_dev14/$CFG; [ -f $d/qref_proc.json ] && args+=("ctl:c_${p}_dev14=$d") || echo "missing $d"
  for d in runs/f_${p}_*/$CFG; do [ -f $d/qref_proc.json ] && args+=("flt:$(basename $(dirname $d))=$d"); done
  out=logs/report_cal_$p
  [ -e $out.json ] && mv $out.json $out.prev.json && mv $out.txt $out.prev.txt      # the previous report is kept once (regenerated)
  T207_3_STRATA=1 T207_3_QPDIR=$W/g5 $PY -I t207_4_cal.py $out.json $p "${args[@]}" > $out.txt 2>&1; echo "report $p rc=$? ($out.txt)"
  [ -f $out.json ] && J+=($out.json)
done
$PY -I t207_4_sep.py "${J[@]}" > logs/report_sep.txt 2>&1; echo "separation rc=$? (logs/report_sep.txt)"
