#!/usr/bin/env python3
"""Ticket 207 dev5 (development only): gate simulation for the proposed statistic S_ens(pTT) pooled over two seed sets.
For every pair of the four dev seed sets and each path: pooled S_ens of the device, of the stand-in control (var8) and of every
correct variant, with a bootstrap over points (sets resampled separately); share of bootstrap means above the bound.
usage: t207_d5_pairs.py BOUND def|log  (reads runs/m_<path>_devN/hub12x12_cb_p/qref.txt)"""
import itertools, math, random, sys
sys.path.insert(0, '.')
import importlib.util
spec = importlib.util.spec_from_file_location('agg', 't207_d5_agg.py'); agg = importlib.util.module_from_spec(spec); spec.loader.exec_module(agg)
B = float(sys.argv[1]); path = sys.argv[2]
per = {}
for s in ('dev1', 'dev2', 'dev3', 'dev4'):
    res, corr = agg.per_run(agg.parse('runs/m_%s_%s/hub12x12_cb_p/qref.txt' % (path, s)), 'pTT')
    per[s] = (res, corr)
rng = random.Random(7)
print('path %s, bound %.2f: S_ens(pTT) pooled over two seed sets' % (path, B))
print('  %-10s %8s %8s %8s %10s %10s %8s' % ('pair', 'dev', 'var8', 'nullmax', 'P(dev>B)', 'P(var8<=B)', 'P(nmax>B)'))
for a, b in itertools.combinations(sorted(per), 2):
    ra, ca = per[a]; rb, cb = per[b]
    nulls = [m for m in ca if m != 'cpu' and m in cb]
    def pooled(m):
        return ra[m]['re'] + rb[m]['re']
    def bmean(m, k=2000):
        va, vb = ra[m]['re'], rb[m]['re']
        out = []
        for _ in range(k):
            x = [va[rng.randrange(len(va))] for _ in va] + [vb[rng.randrange(len(vb))] for _ in vb]
            out.append(sum(x) / len(x))
        return out
    md = sum(pooled('dev')) / len(pooled('dev')); mv = sum(pooled('var8')) / len(pooled('var8'))
    nm = max(sum(pooled(m)) / len(pooled(m)) for m in nulls)
    pd = sum(1 for x in bmean('dev') if x > B) / 2000
    pv = sum(1 for x in bmean('var8') if x <= B) / 2000
    pn = max(sum(1 for x in bmean(m, 1000) if x > B) / 1000 for m in nulls)
    print('  %-10s %+8.3f %+8.3f %+8.3f %10.4f %10.4f %8.4f' % (a + '+' + b, md, mv, nm, pd, pv, pn))
