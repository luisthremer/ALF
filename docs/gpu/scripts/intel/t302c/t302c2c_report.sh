#!/bin/bash
# Ticket 302c-2 round 3 (development only), on intel: all tables, generated from processed runs, per beta and path -> t302c-2c/tab/
#   gates_b<beta>_<path>.{json,log}   t302c2c_gates.py over
#       main runs: round 2's m<beta>_<path>_dev5, _dev6 (t302c-2b/runs, read only; they count for d3e82a36 only after logs/bitid.txt shows
#                  bit identity) and round 3's m<beta>_<path>_dev* (t302c-2c/runs; the bit-identity reruns of dev5 are left out: they are
#                  the same runs)
#       controls:  round 2's f<beta>_<path>_* (08559f0b binaries) and round 3's f<beta>_<path>_* (d3e82a36), processed ones only
#   cal_b<beta>_<path>.{json,txt}     t207_3_cal.py (b2c71200; per run, per pair, flavour subsets), same runs, cross-check
#   summary.md                        t302c2c_summary.py (tables (a)-(d))
cd ~/alf-gpu/work/t302c-2c || exit 1; export TMPDIR=$PWD/tmp
O=$HOME/alf-gpu/work/t302c-2b
P=~/alf-gpu/venv/bin/python3; mkdir -p tab
for b in 4 10; do
  c=hub12x12_cb_p; [ $b = 10 ] && c=hub12x12_cb_b10_p
  for p in def log; do
    g=()
    for d in $O/runs/m${b}_${p}_dev5/$c $O/runs/m${b}_${p}_dev6/$c; do [ -f $d/qref_proc.json ] && g+=("main:r2_$(basename $(dirname $d))=$d"); done
    for d in runs/m${b}_${p}_dev*/$c; do case $d in */m${b}_${p}_dev5/*) continue;; esac
      [ -f $d/qref_proc.json ] && g+=("main:$(basename $(dirname $d))=$PWD/$d"); done
    for d in $O/runs/f${b}_${p}_*/$c; do [ -f $d/qref_proc.json ] && g+=("flt:r2_$(basename $(dirname $d))=$d"); done
    for d in runs/f${b}_${p}_*/$c; do [ -f $d/qref_proc.json ] && g+=("flt:$(basename $(dirname $d))=$PWD/$d"); done
    q=bin/qref_pp.out; [ $p = log ] && q=bin/qref_log.out
    $P -I t302c2c_gates.py tab/gates_b${b}_$p.json b${b}_$p $p $(cut -d' ' -f1 $q.sha256) "${g[@]}" > tab/gates_b${b}_$p.log 2>&1
    echo "gates b${b}_$p rc=$? (${#g[@]} runs)"
    rm -f tab/cal_b${b}_$p.json
    T207_3_STRATA=1 T207_3_QPDIR=$PWD $P -I t207_3_cal.py tab/cal_b${b}_$p.json $p "${g[@]}" > tab/cal_b${b}_$p.txt 2>&1
    echo "cal b${b}_$p rc=$?"
  done
done
$P -I t302c2c_summary.py tab > tab/summary.md; echo "summary rc=$?"
