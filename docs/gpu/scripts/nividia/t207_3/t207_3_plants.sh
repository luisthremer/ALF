#!/bin/bash
# planted inputs for t207_3_summary.py (207-3f); run on nividia.  Reviewer inputs read only from t207-3rev3/in.
F=$HOME/alf-gpu/work/t207-3f; I=$HOME/alf-gpu/work/t207-3rev3/in; P=$F/tmp/pl; PY=$HOME/alf-gpu/venv/bin/python3
mkdir -p $P; cd $F
$PY -I - $I $P <<'E'
import json, sys, copy
I, P = sys.argv[1:3]
d = json.load(open(I + '/report_cal_def_final.json')); l = json.load(open(I + '/report_cal_log_final.json'))
def w(n, o): json.dump(o, open('%s/%s.json' % (P, n), 'w'))
# sd 0: the device pTT mean identical in all four main runs of def
x = copy.deepcopy(d)
for r in [r['label'] for r in x['runs'] if r['role'] == 'main']:
    x['tables']['pTT']['per_run'][r]['dev']['mean'] = 0.05
w('sd0_def', x)
# duplicated main run: two main runs only, the second a copy of the first (same seeds)
x = copy.deepcopy(d)
m = [r for r in x['runs'] if r['role'] == 'main']
dup = copy.deepcopy(m[0]); dup['label'] = m[0]['label'] + 'b'
x['runs'] = [r for r in x['runs'] if r['role'] != 'main'] + [m[0], dup]
w('dup_def', x)
# one main run only
x = copy.deepcopy(d); m = [r for r in x['runs'] if r['role'] == 'main']
x['runs'] = [r for r in x['runs'] if r['role'] != 'main'] + [m[0]]
w('one_def', x)
# LOG: a pair-level device excess is not in the reviewer files; device above the bound in the MAIN POOL (def pTT_pooled)
x = copy.deepcopy(d); x['tables']['pTT']['pooled']['dev']['mean'] = 0.16; w('devpool_def', x)
E
run() { # name def log [env]
  n=$1; shift
  env "$@" > $P/$n.txt 2> $P/$n.err; echo "rc=$?" > $P/$n.rc
}
echo "== baseline (reviewer's unmodified finals)"
$PY -I t207_3_summary.py $P/o_base.json $I/report_cal_def_final.json $I/report_cal_log_final.json > $P/o_base.txt 2>&1; echo "rc=$?"
echo "== P1 device above the bound in one def main run (p1_def: m_def_dev1 pTT 0.16)"
$PY -I t207_3_summary.py $P/o_p1.json $I/p1_def.json $I/report_cal_log_final.json > $P/o_p1.txt 2>&1; echo "rc=$?"
grep -E "^  def pTT_pooled " $P/o_p1.txt | grep -o "correct side [A-Z ]*\|-> [A-Z ]*$" | tr '\n' ' '; echo
echo "== P1b device above the bound in the main pool (devpool_def)"
$PY -I t207_3_summary.py $P/o_p1b.json $P/devpool_def.json $I/report_cal_log_final.json > $P/o_p1b.txt 2>&1; echo "rc=$?"
grep -E "^  def pTT_pooled " $P/o_p1b.txt | grep -o "correct side [A-Z ]*\|-> [A-Z ]*$" | tr '\n' ' '; echo
echo "== P1c device flavour-2 above the flavour-max bound in one LOG run (p3_log)"
$PY -I t207_3_summary.py $P/o_p3.json $I/report_cal_def_final.json $I/p3_log.json > $P/o_p3.txt 2>&1; echo "rc=$?"
grep -E "^  log pTT_flavour_max " $P/o_p3.txt | grep -o "correct side [A-Z ]*\|-> [A-Z ]*$" | tr '\n' ' '; echo
echo "== P2 one LOG control seed set removed (p2_log)"
$PY -I t207_3_summary.py $P/o_p2.json $I/report_cal_def_final.json $I/p2_log.json > $P/o_p2.txt 2>&1; echo "rc=$?"
grep -E "control kind|claim B" $P/o_p2.txt | grep "gate set log" | cut -c1-260
echo "== P3 sd 0 (device pTT identical in the 4 main runs)"
$PY -I t207_3_summary.py $P/o_sd0.json $P/sd0_def.json $I/report_cal_log_final.json > $P/o_sd0.txt 2>&1; echo "rc=$?"; head -c 300 $P/o_sd0.txt
echo "== P4 duplicated main run"
$PY -I t207_3_summary.py $P/o_dup.json $P/dup_def.json $I/report_cal_log_final.json > $P/o_dup.txt 2>&1; echo "rc=$?"; head -c 300 $P/o_dup.txt
echo "== P5 one main run"
$PY -I t207_3_summary.py $P/o_one.json $P/one_def.json $I/report_cal_log_final.json > $P/o_one.txt 2>&1; echo "rc=$?"; head -c 300 $P/o_one.txt
echo "== P6 T207_3_TEST_BOUND unknown path / gate / malformed"
for tb in "xyz:pTT_pooled=0.1" "def:nogate=0.1" "def:pTT_by_stratum=0.1" "def:pTT_pooled" "def:pTT_flavour_max=0.15"; do
  T207_3_TEST_BOUND=$tb $PY -I t207_3_summary.py $P/o_tb.json $I/report_cal_def_final.json $I/report_cal_log_final.json > $P/o_tb.txt 2>&1; echo "  [$tb] rc=$? $(head -c 160 $P/o_tb.txt | head -1)"; done
