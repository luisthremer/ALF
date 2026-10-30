#!/usr/bin/env python3
# Ticket 207 dev4 (rounding-order nulls of the paired test, development only): tables over seed sets, generated from the stat.json
# files of t207_d4_stat.py (no hand-made numbers).
#   usage: t207_d4_agg.py OUT.json LABEL=stat.json [LABEL=stat.json ...]     (one path per call, e.g. dev1=runs/m_def_dev1/.../stat.json)
# Per row (a / b): per set the sd ratio to the gate reference sd(P2/P1), the mean and the three conditions; over the sets min / mean /
# max of the ratio; pooled over all sets: sd of the pooled log ratios divided by the pooled reference sd.  Also the ratio of each row
# to an alternative reference: the median sd of the host-input rounding-order nulls / host (codes without flag 16, code 0 excluded).
import json, math, statistics as st, sys
out, args = sys.argv[1], sys.argv[2:]
sets = []
for a in args:
    lab, fn = a.split('=', 1)
    sets.append((lab, json.load(open(fn))))
keys = []
for _, s in sets:
    for r in s['rows']:
        k = (r['label'], r['a'], r['b'])
        if k not in keys:
            keys.append(k)
def get(s, k):
    for r in s['rows']:
        if (r['label'], r['a'], r['b']) == k:
            return r
    return None
def altref(s):
    v = [r['sd'] for r in s['rows'] if r['label'].startswith('null ') and r['b'] == 'host'
         and int(r['label'].split()[1]) not in (0, 16) and not (int(r['label'].split()[1]) % 100) & 16]
    return st.median(v) if v else float('nan')
res = dict(sets=[l for l, _ in sets], files={l: s['file'] for l, s in sets}, refs={}, rows=[])
for l, s in sets:
    res['refs'][l] = dict(s['refs'], alt_null_median_sd=altref(s), seed=s['seed'], used=s['used'])
for k in keys:
    per = {}
    for l, s in sets:
        r = get(s, k)
        if r is None:
            continue
        per[l] = dict(ratio=r['ratio'], mean=r['mean'], sd=r['sd'], ratio_alt=r['sd'] / altref(s), cond=[r['cond_mean'], r['cond_sd'], r['cond_counts']],
                      counts=r['counts'])
    rs = [v['ratio'] for v in per.values()]
    xa, xr = [], []                                         # pooled over the sets: log ratios of the row and of the reference
    for l, s in sets:
        if get(s, k) is None:
            continue
        xa += [math.log10(p[k[1]] / p[k[2]]) for p in s['series']]
        xr += [math.log10(p['P2'] / p['P1']) for p in s['series']]
    pooled = dict(n=len(xa), mean=st.mean(xa), sd=st.pstdev(xa), ratio=st.pstdev(xa) / st.pstdev(xr))
    res['rows'].append(dict(label=k[0], a=k[1], b=k[2], per_set=per, ratio_min=min(rs), ratio_mean=st.mean(rs), ratio_max=max(rs),
                            pooled=pooled, sd_fail_sets=sum(not v['cond'][1] for v in per.values()),
                            mean_fail_sets=sum(not v['cond'][0] for v in per.values()),
                            counts_fail_sets=sum(not v['cond'][2] for v in per.values())))
json.dump(res, open(out, 'w'), indent=1)
labs = [l for l, _ in sets]
print('references per set: ' + '; '.join('%s P2/P1 %.4f%s, alt %.4f' % (l, res['refs'][l]['P2/P1'],
      ''.join(', %s %.4f' % (k, v) for k, v in res['refs'][l].items() if '/' in k and k != 'P2/P1'), res['refs'][l]['alt_null_median_sd']) for l in labs))
print('%-62s ' % 'row (sd ratio to P2/P1; mean; conditions m/s/c)' + ' '.join('%-24s' % l for l in labs) + '  min   mean   max  | pooled | alt-ref ratio per set')
for r in res['rows']:
    cells = []
    for l in labs:
        v = r['per_set'].get(l)
        cells.append('%-24s' % ('-' if v is None else '%5.3f %+.3f %s' % (v['ratio'], v['mean'], ''.join('.' if c else 'F' for c in v['cond']))))
    alt = ' '.join('%5.2f' % r['per_set'][l]['ratio_alt'] for l in labs if l in r['per_set'])
    print('%-62s ' % (r['label'] + ': ' + r['a'] + '/' + r['b'])[:62] + ' '.join(cells) +
          ' %5.3f %5.3f %5.3f | pooled %5.3f (mean %+.3f, n %d) | %s' % (r['ratio_min'], r['ratio_mean'], r['ratio_max'],
          r['pooled']['ratio'], r['pooled']['mean'], r['pooled']['n'], alt))
