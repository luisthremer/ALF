#!/usr/bin/env python3
"""Ticket 207-2: verdict of the accuracy test against the quad reference from the aggregate of t207_d5_agg.py.

usage: t207_d5_gate.py --bound B [--require-controls] AGG.json        (agg run with measure pTT)
Verdict on the pooled statistic S_ens(pTT) of the runs in AGG.json (one path):
  INPUT_ERROR (exit 2)  the file is not a complete aggregate: wrong measure, a run whose points used differ from the number derived from
                        its parameters, runs with different reference lists / member lists, missing or non-finite numbers, pooled counts
                        that are not the sum of the runs, a main run without 'dev'; with --require-controls: no control member
  FAIL (exit 1)         the device S_ens (member 'dev' of the main runs) is > B, or a control present in the file (the stand-in
                        'var8' of the main runs, 'ctl_dev' of production control runs) has S_ens <= B
  PASS (exit 0)         otherwise
The bound is an argument and is not stored in the aggregate.  Every number of the decision is printed.
"""
import json, math, sys


class InputError(Exception):
    pass


def need(d, k, what):
    if not isinstance(d, dict) or k not in d:
        raise InputError('%s: key %r missing' % (what, k))
    return d[k]


def fin(x, what):
    if not isinstance(x, (int, float)) or isinstance(x, bool) or not math.isfinite(x):
        raise InputError('%s: %r is not a finite number' % (what, x))
    return x


def decide(agg, bound, require_controls):
    out = []
    if need(agg, 'measure', 'aggregate') != 'pTT':
        raise InputError('measure is %r, the gate uses pTT' % agg['measure'])
    refs = need(agg, 'refs', 'aggregate')
    inputs = need(agg, 'inputs', 'aggregate')
    runs = need(agg, 'runs', 'aggregate')
    pooled = need(agg, 'pooled', 'aggregate')
    if not isinstance(refs, list) or 'cpu' not in refs or len(refs) < 2 or len(set(refs)) != len(refs):
        raise InputError('bad reference list %r' % (refs,))
    if not inputs or set(inputs) != set(runs):
        raise InputError('runs and inputs differ')
    main = [l for l, i in inputs.items() if need(i, 'role', 'input ' + l) == 'main']
    if not main:
        raise InputError('no main run')
    members0 = None
    variants = set()
    for lab, i in inputs.items():
        nu, ne, npt = need(i, 'n_used', lab), need(i, 'n_expected', lab), need(i, 'n_points', lab)
        if not (isinstance(nu, int) and nu == ne and nu > 0 and npt >= nu):
            raise InputError('run %s: %r points used, %r derived from its parameters (of %r points)' % (lab, nu, ne, npt))
        mem = need(i, 'members', lab)
        for r in refs:
            if r not in mem:
                raise InputError('run %s: reference member %s missing' % (lab, r))
        if i['role'] == 'main':
            if 'dev' not in mem:
                raise InputError('main run %s has no member dev' % lab)
            if members0 is None:
                members0 = mem
            elif mem != members0:
                raise InputError('main runs %s and %s have different member lists' % (main[0], lab))
        variants.add(need(i, 'variant', lab))
        out.append('run %-8s %-4s points %d used %d expected %d' % (lab, i['role'], npt, nu, ne))
    if len(variants) != 1:
        raise InputError('runs of different paths: %s' % sorted(variants))

    def pooled_s(name):
        e = need(need(pooled, name, 'pooled'), 'S_ens', 'pooled ' + name)
        mean, n = fin(need(e, 'mean', name), name + ' mean'), need(e, 'n', name)
        nexp = sum(runs[l][name]['S_ens']['n'] for l in runs if name in runs[l])
        if n != nexp or n < 1:
            raise InputError('pooled %s: n %r, sum over the runs %r' % (name, n, nexp))
        for l in runs:
            if name in runs[l]:
                fin(runs[l][name]['S_ens']['mean'], '%s %s mean' % (l, name))
        return mean, n

    dev, ndev = pooled_s('dev')
    if ndev != sum(inputs[l]['n_used'] for l in main):
        raise InputError('pooled dev n %d differs from the points of the main runs' % ndev)
    ctl = [n for n in ('var8', 'ctl_dev') if n in pooled]
    if require_controls and not ctl:
        raise InputError('--require-controls: no control member (var8, ctl_dev) in the aggregate')
    verdict = 'PASS'
    reasons = []
    out.append('device S_ens(pTT) pooled over %d points: %+.4f  (bound %+.4f)' % (ndev, dev, bound))
    if not dev <= bound:
        verdict = 'FAIL'; reasons.append('device S_ens %+.4f > bound %+.4f' % (dev, bound))
    for c in ctl:
        v, n = pooled_s(c)
        ok = v > bound
        out.append('control %-7s S_ens(pTT) pooled over %d points: %+.4f  (must be > %+.4f): %s' % (c, n, v, bound, 'killed' if ok else 'NOT KILLED'))
        if not ok:
            verdict = 'FAIL'; reasons.append('control %s S_ens %+.4f <= bound %+.4f' % (c, v, bound))
    return verdict, reasons, out


def main():
    a = sys.argv[1:]
    req = '--require-controls' in a
    a = [x for x in a if x != '--require-controls']
    if len(a) != 3 or a[0] != '--bound':
        sys.stderr.write(__doc__)
        print('GATE INPUT_ERROR: usage')
        sys.exit(2)
    try:
        bound = float(a[1])
        if not math.isfinite(bound):
            raise ValueError
        with open(a[2]) as f:
            agg = json.load(f)
        verdict, reasons, lines = decide(agg, bound, req)
    except (InputError, ValueError, OSError) as e:
        print('GATE INPUT_ERROR: %s' % e)
        sys.exit(2)
    for l in lines:
        print('  ' + l)
    print('GATE %s%s' % (verdict, (': ' + '; '.join(reasons)) if reasons else ''))
    sys.exit(0 if verdict == 'PASS' else 1)


if __name__ == '__main__':
    main()
