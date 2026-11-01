#!/bin/bash
# Ticket 207-3 (development only), nividia: generated calibration tables from the stamped, processed runs (t207_3_cal.py with the strict
# reader of g3/ = gpu-port-207-2 305488ec), per path: all main runs (dev1..dev4), the production control (dev5), the graded faults
# (dev6 near, dev7 far, dev8 localized).  -> logs/report_cal_<path>.{txt,json}.  Runs that do not exist yet are skipped (listed).
W=$HOME/alf-gpu/work/t207-3; cd $W || exit 1
PY=~/alf-gpu/venv/bin/python3
for p in def log; do
  args=()
  for s in dev1 dev2 dev3 dev4; do d=runs/m_${p}_$s/hub12x12_cb_p; [ -f $d/qref_proc.json ] && args+=("main:m_${p}_$s=$d") || echo "missing $d"; done
  d=runs/c_${p}_dev5/hub12x12_cb_p; [ -f $d/qref_proc.json ] && args+=("ctl:c_${p}_dev5=$d") || echo "missing $d"
  for d in runs/f_${p}_*/hub12x12_cb_p; do [ -f $d/qref_proc.json ] && args+=("flt:$(basename $(dirname $d))=$d"); done
  out=logs/report_cal_$p
  rm -f $out.json
  T207_3_STRATA=1 T207_3_QPDIR=$W/g3 $PY -I t207_3_cal.py $out.json $p "${args[@]}" > $out.txt 2>&1; echo "report $p rc=$? ($out.txt)"
done
