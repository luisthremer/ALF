#!/usr/bin/env python3
"""Ticket 207 dev5 (development only): statistics of the accuracy test against the quad reference.

Reads qref.txt files of t207_d5_proc.sh (one 'QREF' line per dump: per method the errors e00 eT0 eTT e0T eF against the quad
reference of the same inputs) and prints, per run and pooled per path, for every method m the per-point log ratios
r = log10(e_m / e_cpu) (cpu = host CGR2_2 of the oracle algorithm on the same inputs) under two error measures
(emax = max of the four block max-relative errors, eF = Frobenius), with bootstrap intervals over points.

usage: t207_d5_stat.py OUT.json LABEL=path/qref.txt [LABEL=...]
"""
import json, math, random, re, statistics, sys

NUM = re.compile(r'[-+]?\d\.\d+E[-+]\d+')


def parse(fn):
    pts = []
    for line in open(fn):
        if not line.startswith('QREF'):
            continue
        segs = line.rstrip('\n').split(' | ')
        h = segs[0].split()
        hd = {h[i]: h[i + 1] for i in range(1, len(h) - 1, 2)}
        m = {}
        for s in segs[1:]:
            t = s.split()
            m[t[0]] = [float(x) for x in NUM.findall(s[len(t[0]):])]
        pts.append((hd, m))
    return pts


def q(v, p):
    v = sorted(v)
    if not v:
        return float('nan')
    k = (len(v) - 1) * p
    i = int(math.floor(k)); j = min(i + 1, len(v) - 1)
    return v[i] + (v[j] - v[i]) * (k - i)


def boot(v, f, n=2000, seed=12345):
    rng = random.Random(seed)
    b = [f([v[rng.randrange(len(v))] for _ in v]) for _ in range(n)]
    return q(b, 0.025), q(b, 0.975)


MEAS = ('emax', 'eF', 'eTT', 'pTT', 'pmax')


def meas(e, which):
    # emax: max of the four block max-relative errors; eF: Frobenius; eTT: GTT block; pTT: max|B dTT Bi| (error carried to the next
    # precision check, absolute); pmax: max of the three propagated numbers (TT, T0, 0T)
    if which == 'emax':
        return max(e[0:4])
    if which == 'eF':
        return e[4]
    if which == 'eTT':
        return e[2]
    if len(e) < 8 or e[5] < 0:
        return None
    return e[5] if which == 'pTT' else max(e[5:8])


def stats(pts, m, which, ref='cpu'):
    r = []
    for hd, d in pts:
        if m in d and ref in d:
            a, b = meas(d[m], which), meas(d[ref], which)
            if a is None or b is None:
                continue
            if a > 0 and b > 0:
                r.append(math.log10(a / b))
            elif a == 0 and b == 0:
                r.append(0.0)
    if not r:
        return None
    em = [meas(d[m], which) for hd, d in pts if m in d and meas(d[m], which) is not None]
    out = dict(n=len(r), mean=statistics.fmean(r), median=q(r, 0.5), q90=q(r, 0.9), q10=q(r, 0.1),
               sd=statistics.pstdev(r), frac_gt2=sum(1 for x in r if x > math.log10(2)) / len(r),
               frac_gt3=sum(1 for x in r if x > math.log10(3)) / len(r),
               gm_err=10 ** statistics.fmean([math.log10(max(x, 1e-300)) for x in em]), med_err=q(em, 0.5), max_err=max(em))
    out['mean_ci'] = boot(r, statistics.fmean)
    out['median_ci'] = boot(r, lambda v: q(v, 0.5))
    return out


def main():
    outf = sys.argv[1]
    runs = {}
    for a in sys.argv[2:]:
        lab, fn = a.split('=', 1)
        runs[lab] = parse(fn)
    res = {}
    for lab, pts in runs.items():
        methods = sorted({k for hd, d in pts for k in d})
        res[lab] = {}
        for which in MEAS:
            res[lab][which] = {}
            for m in methods:
                if m in ('qself', 'qbr', 'din', 'cpuh', 'cpuhx'):
                    em = [meas(d[m], which) for hd, d in pts if m in d and meas(d[m], which) is not None]
                    if em:
                        res[lab][which][m] = dict(n=len(em), med=q(em, 0.5), max=max(em))
                    continue
                s = stats(pts, m, which)
                if s:
                    res[lab][which][m] = s
    json.dump(res, open(outf, 'w'), indent=1)
    for lab in res:
        for which in MEAS:
            print('== %s  measure %s' % (lab, which))
            print('  %-8s %4s %8s %8s %8s %8s %6s %6s %6s %10s %10s %10s  %s' % ('method', 'n', 'mean_r', 'med_r', 'q90_r', 'sd_r', 'f>2x', 'f>3x', '', 'gm_err', 'med_err', 'max_err', 'mean_r 95% CI'))
            for m, s in res[lab][which].items():
                if 'mean' not in s:
                    print('  %-8s %4d %s med %.3e max %.3e' % (m, s['n'], ' ' * 46, s['med'], s['max']))
                    continue
                print('  %-8s %4d %+8.3f %+8.3f %+8.3f %8.3f %6.3f %6.3f %6s %10.3e %10.3e %10.3e  [%+.3f, %+.3f]' % (
                    m, s['n'], s['mean'], s['median'], s['q90'], s['sd'], s['frac_gt2'], s['frac_gt3'], '', s['gm_err'],
                    s['med_err'], s['max_err'], s['mean_ci'][0], s['mean_ci'][1]))


if __name__ == '__main__':
    main()
