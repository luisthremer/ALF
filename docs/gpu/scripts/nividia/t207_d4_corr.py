#!/usr/bin/env python3
# Ticket 207 dev4 (development only): correctness of the null chains over all runs, from the stat.json files of t207_d4_stat.py:
# per chain the median over runs of the per-run median, and the max over runs of the per-run max, of Nd (max |GTT_chain - GTT_host| /
# max |GTT_host| of the new GTT at the same point and input), and the bitwise-equal counts; the device (Ad) and the 1-ulp chains alongside.
#   usage: t207_d4_corr.py stat.json [...]
import json, statistics as st, sys
agg = {}
for fn in sys.argv[1:]:
    s = json.load(open(fn))
    path = 'log' if '/m_log_' in fn else 'def'
    for k, v in s['correctness'].items():
        a = agg.setdefault((path, k), dict(med=[], mx=[], zero=0, n=0))
        a['med'].append(v['median']); a['mx'].append(v['max']); a['zero'] += v['zero']; a['n'] += v['n']
print('%-5s %-14s %12s %12s %16s' % ('path', 'chain', 'median(med)', 'max(max)', 'bitwise equal'))
for (p, k), a in agg.items():
    print('%-5s %-14s %12.2e %12.2e %9d/%-6d' % (p, k, st.median(a['med']), max(a['mx']), a['zero'], a['n']))
