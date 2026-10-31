#!/usr/bin/env python3
"""Ticket 207 dev5 (development only): calibration table of the accuracy test against the quad reference.

Per run (seed set) and pooled over the runs of one path, for every method m and a measure (default pTT: max|B dGTT Bi|, the error
of the CGR2_2 output as the next precision check receives it; emax: pointwise max-relative error):
  S_cpu(m) = mean over points of log10(e_m / e_cpu)            (reference: the host CGR2_2 of the oracle algorithm, same inputs)
  S_ens(m) = mean over points of log10(e_m / g)                 (g: geometric mean of e over the correct host variants at the point,
                                                                  cpu, nvh and the null codes except n0; leave-one-out for a member)
with a bootstrap 95 % interval over points, and the share of bootstrap means above each candidate bound.

usage: t207_d5_agg.py OUT.json MEASURE ROLE:LABEL=path/qref.txt ...
  ROLE: main (device = 'dev', stand-in control = 'var8') | ctl (production control run: its 'dev' is the defect)
"""
import json, math, random, re, statistics, sys

NUM = re.compile(r'[-+]?\d\.\d+E[-+]\d+')
CORRECT_EXCL = {'dev', 'var8', 'n0', 'u1', 'u2', 'qself', 'qbr', 'din', 'cpuh', 'cpuhx'}


def parse(fn):
    pts = []
    for line in open(fn):
        if not line.startswith('QREF'):
            continue
        segs = line.rstrip('\n').split(' | ')
        m = {}
        for s in segs[1:]:
            t = s.split()
            m[t[0]] = [float(x) for x in NUM.findall(s[len(t[0]):])]
        pts.append(m)
    return pts


def meas(e, which):
    if which == 'emax':
        return max(e[0:4])
    if which == 'eF':
        return e[4]
    if len(e) < 8 or e[5] < 0:
        return None
    return e[5] if which == 'pTT' else max(e[5:8])


def q(v, p):
    v = sorted(v)
    k = (len(v) - 1) * p
    i = int(math.floor(k)); j = min(i + 1, len(v) - 1)
    return v[i] + (v[j] - v[i]) * (k - i)


def boot_means(v, n=4000, seed=20261030):
    rng = random.Random(seed)
    L = len(v)
    return [sum(v[rng.randrange(L)] for _ in range(L)) / L for _ in range(n)]


def logs(pts, which):
    """per point: dict method -> log10 e (only points where every method has the measure)"""
    out = []
    for d in pts:
        r = {}
        ok = True
        for k, e in d.items():
            if k in ('qself', 'qbr', 'din', 'cpuh', 'cpuhx'):
                continue
            x = meas(e, which)
            if x is None:
                ok = False; break
            r[k] = math.log10(max(x, 1e-300))
        if ok and 'cpu' in r:
            out.append(r)
    return out


def per_run(pts, which):
    L = logs(pts, which)
    methods = sorted(L[0]) if L else []
    corr = [m for m in methods if m not in CORRECT_EXCL]
    res = {}
    for m in methods:
        rc = [p[m] - p['cpu'] for p in L]
        re_ = []
        for p in L:
            members = [k for k in corr if k != m]
            g = sum(p[k] for k in members) / len(members)
            re_.append(p[m] - g)
        res[m] = dict(n=len(rc), rc=rc, re=re_)
    return res, corr


def summ(v):
    b = boot_means(v)
    return dict(n=len(v), mean=statistics.fmean(v), lo=q(b, 0.025), hi=q(b, 0.975), sd=statistics.pstdev(v), boot=b)


def main():
    outf, which = sys.argv[1], sys.argv[2]
    runs = []
    for a in sys.argv[3:]:
        role, rest = a.split(':', 1)
        lab, fn = rest.split('=', 1)
        runs.append((role, lab, parse(fn)))
    table = {}
    pooled = {}
    corr_all = set()
    for role, lab, pts in runs:
        res, corr = per_run(pts, which)
        corr_all |= set(corr)
        table[lab] = {}
        for m, d in res.items():
            name = m
            if role == 'ctl' and m == 'dev':
                name = 'ctl_dev'
            elif role == 'ctl':
                continue                            # the other methods of a control run are not part of the calibration
            sc, se = summ(d['rc']), summ(d['re'])
            table[lab][name] = dict(S_cpu=sc, S_ens=se)
            pooled.setdefault(name, dict(rc=[], re=[]))
            pooled[name]['rc'] += d['rc']; pooled[name]['re'] += d['re']
    pool = {m: dict(S_cpu=summ(d['rc']), S_ens=summ(d['re'])) for m, d in pooled.items()}
    # print
    labs = [lab for _, lab, _ in runs]
    order = ['dev', 'var8', 'ctl_dev', 'cpu', 'nvh'] + sorted(m for m in pool if m.startswith('n') and m not in ('nvh', 'n0')) + ['u1', 'u2']
    for ref in ('S_cpu', 'S_ens'):
        print('== measure %s, statistic %s (mean log10 ratio; [95%% bootstrap interval over points])' % (which, ref))
        print('  %-9s' % 'method' + ''.join('%24s' % l for l in labs) + '%24s' % 'pooled')
        for m in order:
            if m not in pool:
                continue
            row = '  %-9s' % m
            for l in labs:
                s = table[l].get(m, {}).get(ref)
                row += '%24s' % ('%+.3f [%+.3f,%+.3f]' % (s['mean'], s['lo'], s['hi']) if s else '-')
            s = pool[m][ref]
            row += '%24s' % ('%+.3f [%+.3f,%+.3f]' % (s['mean'], s['lo'], s['hi']))
            print(row)
        # separation and false-alarm view
        nulls = [m for m in pool if m in corr_all and m != 'cpu']
        per_set_null_max = {l: max(table[l][m][ref]['mean'] for m in nulls if m in table[l]) for l in labs if any(m in table[l] for m in nulls)}
        print('  max over correct variants per set: ' + ', '.join('%s %+.3f' % (l, v) for l, v in per_set_null_max.items()))
        for cand in (0.10, 0.15, 0.20, 0.25):
            fa = []
            for l in labs:
                for m in nulls:
                    if m in table[l]:
                        b = table[l][m][ref]['boot']
                        fa.append(sum(1 for x in b if x > cand) / len(b))
            dv = [table[l]['dev'][ref]['mean'] for l in labs if 'dev' in table[l]]
            cv = [table[l][c][ref]['mean'] for l in labs for c in ('var8', 'ctl_dev') if c in table[l]]
            print('  bound %.2f: correct variants: mean bootstrap exceedance %.4f, max %.4f, (set, variant) pairs %d, of them with mean above the bound: %d;'
                  ' device per set max %+.3f (k_min %.2f); controls per set min %+.3f' % (
                      cand, statistics.fmean(fa) if fa else float('nan'), max(fa) if fa else float('nan'), len(fa),
                      sum(1 for l in labs for m in nulls if m in table[l] and table[l][m][ref]['mean'] > cand),
                      max(dv) if dv else float('nan'), 10 ** (cand - max(dv)) if dv else float('nan'), min(cv) if cv else float('nan')))
    for l in table:
        for m in table[l]:
            for ref in ('S_cpu', 'S_ens'):
                table[l][m][ref].pop('boot')
    for m in pool:
        for ref in ('S_cpu', 'S_ens'):
            pool[m][ref].pop('boot')
    json.dump(dict(measure=which, runs=table, pooled=pool), open(outf, 'w'), indent=1)


if __name__ == '__main__':
    main()
