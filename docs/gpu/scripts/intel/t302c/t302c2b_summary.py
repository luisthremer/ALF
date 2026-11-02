#!/usr/bin/env python3
"""Ticket 302c-2 round 2 (development only): compact markdown tables generated from tab/gates_<tag>.json (t302c2b_gates.py), for the report.
Per tag and candidate (measure / statistic / pool): device, correct-variant maximum (member), var8, every fault control, and per control
the decision-37 range max(device, correct max) <= b < control, or '-' (no bound separates).  No number is typed by hand.
usage: t302c2b_summary.py TABDIR [CANDIDATE ...]   (default candidates below)"""
import json, os, sys

D = sys.argv[1]
CANDS = sys.argv[2:] or ['pTT/S_ens/pooled', 'pTT/S_ens/each_run', 'pTT/S_ens/by_stratum', 'pTT/S_ens_q95/pooled', 'pTT/S_ens_max/pooled',
                         'emax/S_ens/pooled', 'emax/S_ens/each_run', 'emax/S_ens/by_stratum', 'emax/S_ens_q95/pooled', 'emax/S_ens_max/pooled',
                         'eF/S_ens/pooled', 'eF/S_ens/by_stratum', 'eF/S_ens_q95/pooled', 'eF/S_ens_max/pooled']
ORDER = ['2em15', '1em14', '1em12', '1em10', '2em12_at2x2']
NAME = {'2em15': 'eps 2e-15', '1em14': 'eps 1e-14', '1em12': 'eps 1e-12', '1em10': 'eps 1e-10', '2em12_at2x2': 'loc 2e-12@(2,2)'}
for tag in ('b4_def', 'b4_log', 'b10_def', 'b10_log'):
    f = os.path.join(D, 'gates_%s.json' % tag)
    if not os.path.exists(f):
        continue
    o = json.load(open(f))
    runs = [r['label'] for r in o['runs']]
    print('\n### %s  (runs: %s; ALF binaries %s; binary128 program %s)\n' % (tag, ', '.join(runs), ' '.join(b[:16] for b in o['alf_binaries']), o['prog_sha256'][:16]))
    ctl_labels = {}
    for r in o['runs']:
        if r['role'] == 'flt':
            for k in ORDER:
                if '_%s_' % k in r['label']:
                    ctl_labels[k] = r['label']
    hdr = '| candidate | device | correct max (member) | var8 | ' + ' | '.join(NAME[k] for k in ORDER if k in ctl_labels) + ' |'
    print(hdr)
    print('|' + '---|' * (hdr.count('|') - 1))
    for c in CANDS:
        x = o['candidates'][c]
        cells = []
        for k in ORDER:
            if k not in ctl_labels:
                continue
            ct = x['controls'][ctl_labels[k]]
            cells.append('%+.3f (%s)' % (ct['value'], ('%+.3f..%+.3f' % tuple(ct['range'])) if ct['separates'] else 'none'))
        print('| %s | %+.3f | %+.3f (%s) | %+.3f | %s |' % (c, x['device']['value'], x['correct_max']['value'], x['correct_max']['member'], x['var8']['value'], ' | '.join(cells)))
    print('\nper main run S_ens(pTT), top correct members: ' + '; '.join(
        '%s: %s' % (l, ', '.join('%s %+.3f' % kv for kv in sorted(((k, v) for k, v in d.items() if k in o['refs']), key=lambda kv: -kv[1])[:4]))
        for l, d in o['members_pTT_per_run'].items()))
