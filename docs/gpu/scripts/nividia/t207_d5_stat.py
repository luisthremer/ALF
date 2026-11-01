#!/usr/bin/env python3
"""Ticket 207 dev5 (development only): statistics of the accuracy test against the quad reference.

Reads qref.txt files of t207_d5_proc.sh through the strict parser t207_d5_qparse.py (one 'QREF' line per dump: per method the
errors e00 eT0 eTT e0T eF and three propagated numbers against the quad reference of the same inputs; a deviation from the expected
content is an error, nothing is skipped) and prints, per run and pooled per path, for every method m the per-point log ratios
r = log10(e_m / e_cpu) (cpu = host CGR2_2 of the oracle algorithm on the same inputs) under two error measures
(emax = max of the four block max-relative errors, eF = Frobenius), with bootstrap intervals over points.

usage: t207_d5_stat.py OUT.json [--allow-v1] LABEL=path/qref.txt [LABEL=...]
"""
import importlib.util, json, math, os, random, re, statistics, sys
sys.dont_write_bytecode = True      # the loaders below import sibling scripts by path: no __pycache__ in the scripts directory

_spec = importlib.util.spec_from_file_location('t207_d5_qparse', os.path.join(os.path.dirname(os.path.abspath(__file__)), 't207_d5_qparse.py'))
QP = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(QP)
QError = QP.QError


def parse(fn, allow_v1=False):
    """validated run -> (list of (header dict, members dict), number of points with B derived from the run parameters)"""
    run = QP.load(fn, allow_v1)
    return [(p.hd, p.m) for p in run.points], run.n_with_b


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
    if len(e) != 8 or e[5] < 0:
        raise QError('measure %s needs the propagated numbers (point without B)' % which)
    return e[5] if which == 'pTT' else max(e[5:8])


def stats(pts, m, which, nb, ref='cpu'):
    """per-point log ratios of method m to ref.  pTT / pmax: only the points with B (hd['hasb'] == 1; their number must equal the number
    derived from the run parameters, nb); a zero or non-finite error number is an error (nothing is skipped)."""
    sel = [(hd, d) for hd, d in pts if (hd['hasb'] == 1 or which not in ('pTT', 'pmax'))]
    if which in ('pTT', 'pmax') and len(sel) != nb:
        raise QError('%d points with B, %d derived from the run parameters' % (len(sel), nb))
    if not sel:
        return None
    r = []
    em = []
    for hd, d in sel:
        if m not in d or ref not in d:
            raise QError('dump %d: member %s or %s missing' % (hd['file'], m, ref))
        a, b = meas(d[m], which), meas(d[ref], which)
        if not (math.isfinite(a) and math.isfinite(b) and a > 0 and b > 0):
            raise QError('dump %d: error number of %s or %s is zero or not finite (%r, %r)' % (hd['file'], m, ref, a, b))
        r.append(math.log10(a / b))
        em.append(a)
    out = dict(n=len(r), mean=statistics.fmean(r), median=q(r, 0.5), q90=q(r, 0.9), q10=q(r, 0.1),
               sd=statistics.pstdev(r), frac_gt2=sum(1 for x in r if x > math.log10(2)) / len(r),
               frac_gt3=sum(1 for x in r if x > math.log10(3)) / len(r),
               gm_err=10 ** statistics.fmean([math.log10(x) for x in em]), med_err=q(em, 0.5), max_err=max(em))
    out['mean_ci'] = boot(r, statistics.fmean)
    out['median_ci'] = boot(r, lambda v: q(v, 0.5))
    return out


def main():
    args = sys.argv[1:]
    allow_v1 = '--allow-v1' in args
    args = [a for a in args if a != '--allow-v1']
    if len(args) < 2:
        raise QError('usage: t207_d5_stat.py OUT.json [--allow-v1] LABEL=path/qref.txt [LABEL=...]')
    outf = args[0]
    runs = {}
    for a in args[1:]:
        if '=' not in a:
            raise QError('bad run argument %r (LABEL=path)' % a)
        lab, fn = a.split('=', 1)
        if lab in runs:
            raise QError('duplicate label %s' % lab)
        runs[lab] = parse(fn, allow_v1)
    res = {}
    for lab, (pts, nb) in runs.items():
        methods = sorted({k for hd, d in pts for k in d})
        res[lab] = {}
        for which in MEAS:
            res[lab][which] = {}
            for m in methods:
                if m in QP.DIAG_NAMES:
                    # members that are no error statistic of a G (quad self-checks: 5 numbers only; input part: 8): their median / max
                    # over the points they exist at; the propagated measures exist only for the 8-number members, at points with B
                    pr = which in ('pTT', 'pmax')
                    if pr and m in ('qself', 'qbr'):
                        continue
                    em = [meas(d[m], which) for hd, d in pts if m in d and (hd['hasb'] == 1 or not pr)]
                    if em:
                        res[lab][which][m] = dict(n=len(em), med=q(em, 0.5), max=max(em))
                    continue
                s = stats(pts, m, which, nb)
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
    try:
        main()
    except QError as e:
        sys.stderr.write('t207_d5_stat ERROR: %s\n' % e)
        sys.exit(1)
