#!/bin/bash
# Ticket 207 dev5 (accuracy test against a quad reference, development only), nividia: tables of the calibration, generated from the
# qref.txt files of the runs (t207_d5_agg.py): per path and measure (pTT, pmax, emax, eF) over the seed sets dev1..dev4 (role main:
# device 'dev', stand-in control 'var8'), plus the production control runs (role ctl) when present; the 206 input part (din runs);
# the production-identity check; GPU minutes of all dev5 GPU runs.  Output: dev5/logs/report_*.txt / *.json
W=$HOME/alf-gpu/work/t207/dev5; cd $W || exit 1
P=~/alf-gpu/venv/bin/python3
for p in def log; do
  args=()
  for s in dev1 dev2 dev3 dev4; do
    f=runs/m_${p}_$s/hub12x12_cb_p/qref.txt; [ -s $f ] && args+=("main:$s=$f")
  done
  f=runs/c_${p}_dev1/hub12x12_cb_p/qref.txt; [ -s $f ] && args+=("ctl:ctl1=$f")
  for m in pTT pmax emax eF; do
    $P -I t207_d5_agg.py logs/report_agg_${p}_$m.json $m "${args[@]}" > logs/report_agg_${p}_$m.txt 2>&1
  done
  f=runs/d_${p}_dev1/hub12x12_cb_p/qref.txt
  [ -s $f ] && $P -I t207_d5_stat.py logs/report_din_${p}.json din=$f > logs/report_din_${p}.txt 2>&1
done
{
  echo "== production-like identity (no debug variable), working tree vs base c1cc528e"
  grep h5cmp logs/chain_ident.txt 2>/dev/null
  grep -h "GPU debug env" runs/id_*/hub12x12_cb_b10_s/run.log 2>/dev/null | sort | uniq -c
  echo "== quad self-check over all processed runs (qself: row-permuted QR input; qbr: other branch), max over blocks"
  $P -I - <<'PY'
import glob, re
NUM = re.compile(r'[-+]?\d\.\d+E[-+]\d+')
for k in ('qself', 'qbr'):
    v = []
    for f in sorted(glob.glob('runs/*/hub12x12_cb_p/qref.txt')):
        for line in open(f):
            for s in line.split(' | ')[1:]:
                if s.startswith(k + ' '):
                    v.append(max(float(x) for x in NUM.findall(s[len(k):])[:4]))
    print('  %s: n %d, max %.3e' % (k, len(v), max(v) if v else float('nan')))
t = []
for f in sorted(glob.glob('runs/*/hub12x12_cb_p/qref.txt')):
    for line in open(f):
        m = re.search(r' tq +([0-9.]+)', line)
        if m: t.append(float(m.group(1)))
t.sort()
print('  quad CGR2_2 wall per call (4 threads, N = 144): n %d, median %.2f s, max %.2f s' % (len(t), t[len(t)//2], t[-1]))
PY
  echo "== 206 input part (TAUDIAG=2 runs d_<path>_dev1, every 4th point): geometric means; din = quad G(device input) - quad G(host"
  echo "   chain input); cpu / cpuh = host CGR2_2 on the device / host input against its own quad reference; dev = device CGR2_2"
  $P -I - <<'PY'
import glob, math, re, statistics
NUM = re.compile(r'[-+]?\d\.\d+E[-+]\d+')
for f in sorted(glob.glob('runs/d_*_dev1/hub12x12_cb_p/qref.txt')):
    rows = []
    for line in open(f):
        d = {}
        for s in line.split(' | ')[1:]:
            t = s.split(); d[t[0]] = [float(x) for x in NUM.findall(s[len(t[0]):])]
        rows.append(d)
    print('  ' + f.split('/')[1])
    for m in ('dev', 'cpu', 'cpuh', 'din'):
        out = []
        for i, lab in ((4, 'eF'), (5, 'pTT')):
            v = [r[m][i] for r in rows if m in r and len(r[m]) > i and r[m][i] > 0]
            out.append('%s n %3d gm %.3e max %.3e' % (lab, len(v), 10 ** statistics.fmean(math.log10(x) for x in v), max(v)))
        print('    %-5s %s' % (m, ' | '.join(out)))
    r = [math.log10(x['cpu'][5] / x['cpuh'][5]) for x in rows if len(x['cpu']) > 5 and x['cpu'][5] > 0 and x['cpuh'][5] > 0]
    print('    mean log10(pTT cpu on device input / pTT cpuh on host input) %+.3f (n %d)' % (statistics.fmean(r), len(r)))
PY
  echo "== offline processing (CPU, 4 threads)"
  for f in runs/*/hub12x12_cb_p/qref_proc.json; do $P -I -c "import json,sys; d=json.load(open(sys.argv[1])); print('  %-40s rc %d wall %5d s dumps %d lines %d' % (sys.argv[1].split('/')[1], d['rc'], d['wall_s'], d['dumps'], d['lines']))" $f; done
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
echo report done
