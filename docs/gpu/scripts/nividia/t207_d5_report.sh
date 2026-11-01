#!/bin/bash
# Ticket 207 dev5 (accuracy test against a quad reference, development only), nividia: tables of the calibration, generated from the
# qref.txt files of the runs (t207_d5_agg.py): per path and measure (pTT, pmax, emax, eF) over the seed sets dev1..dev4 (role main:
# device 'dev', stand-in control 'var8'), plus the production control runs (role ctl) when present; the 206 input part (din runs);
# the production-identity check; GPU minutes of all dev5 GPU runs.  Output: dev5/logs/report_*.txt / *.json
# Ticket 207-2: every qref.txt is read through the strict parser t207_d5_qparse.py (a run whose qref.txt fails it makes the report fail,
# the error is in the output), the reference ensemble is the explicit list 'def' / 'log' (t207_d5_agg.py --ref), v1 dumps (smoke runs
# s_*) are accepted only in the quad self-check / timing summary.  Exit status 0 only if nothing failed.  [D5W=workdir]
W=${D5W:-$HOME/alf-gpu/work/t207/dev5}; cd $W || exit 1
P=~/alf-gpu/venv/bin/python3
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); export QPARSE=$HERE/t207_d5_qparse.py
rc=0
for p in def log; do
  args=()
  for s in dev1 dev2 dev3 dev4; do
    f=runs/m_${p}_$s/hub12x12_cb_p/qref.txt; [ -s $f ] && args+=("main:$s=$f")
  done
  f=runs/c_${p}_dev1/hub12x12_cb_p/qref.txt; [ -s $f ] && args+=("ctl:ctl1=$f")
  for m in pTT pmax emax eF; do
    $P -I $HERE/t207_d5_agg.py logs/report_agg_${p}_$m.json $m --ref $p "${args[@]}" > logs/report_agg_${p}_$m.txt 2>&1 || { echo "report FAILED: agg $p $m (see logs/report_agg_${p}_$m.txt)"; rc=1; }
  done
  f=runs/d_${p}_dev1/hub12x12_cb_p/qref.txt
  if [ -s $f ]; then
    $P -I $HERE/t207_d5_stat.py logs/report_din_${p}.json din=$f > logs/report_din_${p}.txt 2>&1 || { echo "report FAILED: stat din $p"; rc=1; }
  fi
done
{
  echo "== production-like identity (no debug variable), working tree vs base c1cc528e"
  grep h5cmp logs/chain_ident.txt 2>/dev/null
  grep -h "GPU debug env" runs/id_*/hub12x12_cb_b10_s/run.log 2>/dev/null | sort | uniq -c
  echo "== quad self-check over all processed runs (qself: row-permuted QR input; qbr: other branch), max over blocks"
  $P -I - <<'PY'
import glob, importlib.util, os, sys
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('qp', os.environ['QPARSE']); QP = importlib.util.module_from_spec(spec); spec.loader.exec_module(QP)
runs = {}
for f in sorted(glob.glob('runs/*/hub12x12_cb_p/qref.txt')):
    tag = f.split('/')[1]
    try:
        runs[tag] = QP.load(f, allow_v1=tag.startswith('s_'))
    except QP.QError as e:
        print('  %s: QPARSE ERROR: %s' % (tag, e)); sys.exit(1)
for k in ('qself', 'qbr'):
    v = [max(p.m[k][:4]) for r in runs.values() for p in r.points if k in p.m]
    print('  %s: n %d, max %.3e' % (k, len(v), max(v) if v else float('nan')))
t = sorted(p.hd['tq'] for r in runs.values() for p in r.points)
print('  quad CGR2_2 wall per call (4 threads, N = 144): n %d, median %.2f s, max %.2f s' % (len(t), t[len(t)//2], t[-1]))
PY
  echo "== 206 input part (TAUDIAG=2 runs d_<path>_dev1, every 4th point): geometric means; din = quad G(device input) - quad G(host"
  echo "   chain input); cpu / cpuh = host CGR2_2 on the device / host input against its own quad reference; dev = device CGR2_2"
  $P -I - <<'PY'
import glob, importlib.util, math, os, statistics, sys
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('qp', os.environ['QPARSE']); QP = importlib.util.module_from_spec(spec); spec.loader.exec_module(QP)
for f in sorted(glob.glob('runs/d_*_dev1/hub12x12_cb_p/qref.txt')):
    tag = f.split('/')[1]
    try:
        run = QP.load(f)
    except QP.QError as e:
        print('  %s: QPARSE ERROR: %s' % (tag, e)); sys.exit(1)
    print('  ' + tag)
    for m in ('dev', 'cpu', 'cpuh', 'din'):
        out = []
        for i, lab in ((4, 'eF'), (5, 'pTT')):
            v = [p.m[m][i] for p in run.points if m in p.m and (i < 5 or p.hasb)]
            if not v or min(v) <= 0:
                print('  %s: QPARSE ERROR: member %s %s: no point or a non-positive number' % (tag, m, lab)); sys.exit(1)
            out.append('%s n %3d gm %.3e max %.3e' % (lab, len(v), 10 ** statistics.fmean(math.log10(x) for x in v), max(v)))
        print('    %-5s %s' % (m, ' | '.join(out)))
    r = [math.log10(p.m['cpu'][5] / p.m['cpuh'][5]) for p in run.points if p.hasb]
    print('    mean log10(pTT cpu on device input / pTT cpuh on host input) %+.3f (n %d)' % (statistics.fmean(r), len(r)))
PY
  echo "== offline processing (CPU, 4 threads)"
  for f in runs/*/hub12x12_cb_p/qref_proc.json; do $P -I -c "import json,sys; d=json.load(open(sys.argv[1])); print('  %-40s rc %d check_rc %s wall %5d s dumps %d lines %d' % (sys.argv[1].split('/')[1], d['rc'], d.get('check_rc', 'n/a (older pipeline)'), d['wall_s'], d['dumps'], d['lines']))" $f; done
  echo "== GPU minutes (wall of the locked section, all dev5 GPU runs)"
  $P -I - <<'PY'
import glob, json
t = 0
for f in sorted(glob.glob('runs/*/*/provenance.json')):
    p = json.load(open(f))
    if p['build'].startswith(('gnu', 'nv_none')):
        continue
    t += p['wall_s']; print('  %-40s %7.1f s rc %d' % (f.split('/')[1] + '/' + f.split('/')[2], p['wall_s'], p['rc']))
print('  total %.1f s = %.1f GPU-min' % (t, t / 60))
PY
} > logs/report_checks.txt 2>&1
if grep -q "QPARSE ERROR\|Traceback" logs/report_checks.txt; then echo "report FAILED: see logs/report_checks.txt"; rc=1; fi
if [ $rc -eq 0 ]; then echo "report done"; else echo "report finished WITH ERRORS"; fi
exit $rc
