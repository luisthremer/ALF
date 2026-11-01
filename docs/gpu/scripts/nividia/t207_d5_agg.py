#!/usr/bin/env python3
"""Ticket 207 dev5 (development only): calibration table of the accuracy test against the quad reference.

Per run (seed set) and pooled over the runs of one path, for every method m and a measure (default pTT: max|B dGTT Bi|, the error
of the CGR2_2 output as the next precision check receives it; emax: pointwise max-relative error):
  S_cpu(m) = mean over points of log10(e_m / e_cpu)            (reference: the host CGR2_2 of the oracle algorithm, same inputs)
  S_ens(m) = mean over points of log10(e_m / g)                 (g: geometric mean of e over the members of the --ref list at the
                                                                  point; leave-one-out for a member)
with a bootstrap 95 % interval over points, and the share of bootstrap means above each candidate bound.

usage: t207_d5_agg.py OUT.json MEASURE --ref def|log|LIST [--allow-v1] ROLE:LABEL=path/qref.txt ...
  ROLE: main (device = 'dev', stand-in control = 'var8') | ctl (production control run: its 'dev' is the defect)
  --ref: the reference ensemble, an explicit list (ticket 207-2): 'def' = cpu,nvh,n4,n100004,n1,n100001,n2,n1608,n4808,n5,n100005,
         'log' = cpu,nvh,n4,n100004,n1,n100001,n2,n5,n100005, or a comma list of names cpu, nvh, n<code> (code != 0).  An unknown
         name, a duplicate, or a listed name missing from any run is an error (no method is ever treated as correct by default).
Input is read through the strict parser t207_d5_qparse.py (every point, every member, every number checked); any deviation ends the
run with an error and no output file.  Measures pTT / pmax use the points with B (their number must equal the number derived from
the run parameters); emax / eF use all points.  An error number that is zero or not finite is an error (no clamping).
"""
import importlib.util, json, math, os, random, re, statistics, sys
sys.dont_write_bytecode = True      # the loaders below import sibling scripts by path: no __pycache__ in the scripts directory

_spec = importlib.util.spec_from_file_location('t207_d5_qparse', os.path.join(os.path.dirname(os.path.abspath(__file__)), 't207_d5_qparse.py'))
QP = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(QP)
QError = QP.QError
REFS = {'def': QP.REF_DEF, 'log': QP.REF_LOG}


def ref_list(arg):
    """the explicit reference ensemble: validated names, in the given order"""
    names = REFS.get(arg, arg).split(',')
    seen = set()
    for n in names:
        if not (n in ('cpu', 'nvh') or re.match(r'^n[1-9][0-9]*$', n)):
            raise QError('--ref: %r is not a reference name (allowed: cpu, nvh, n<code> with code != 0)' % n)
        if n in seen:
            raise QError('--ref: duplicate name %s' % n)
        seen.add(n)
    if 'cpu' not in seen or len(names) < 2:
        raise QError('--ref: the list must contain cpu and at least one more member')
    return names


def parse(fn, allow_v1=False):
    """the validated run of one qref.txt: list of (Point), n points with B"""
    run = QP.load(fn, allow_v1)
    return run


def meas(e, which):
    if which == 'emax':
        return max(e[0:4])
    if which == 'eF':
        return e[4]
    if which == 'eTT':
        return e[2]
    if len(e) != 8 or e[5] < 0:
        raise QError('measure %s needs the propagated numbers (point without B)' % which)
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


def logs(run, which):
    """per point: dict method -> log10 e, for the points the measure exists for (all, or those with B); returns (list, derived count)"""
    out = []
    for pt in run.points:
        if which in ('pTT', 'pmax') and not pt.hasb:
            continue
        r = {}
        for k, e in pt.m.items():
            if k in QP.DIAG_NAMES:
                continue
            x = meas(e, which)
            if not (math.isfinite(x) and x > 0):
                raise QError('dump %d member %s: error number %r is not positive and finite (log undefined; no clamping)' % (pt.hd['file'], k, x))
            r[k] = math.log10(x)
        out.append(r)
    n_exp = run.n_with_b if which in ('pTT', 'pmax') else len(run.points)
    if len(out) != n_exp or not out:
        raise QError('%d points used for measure %s, %d derived from the run parameters' % (len(out), which, n_exp))
    return out, n_exp


def per_run(run, which, refs):
    L, n_exp = logs(run, which)
    methods = list(L[0])
    for p in L:
        if list(p) != methods:
            raise QError('the member list differs between points')
    for r in refs:
        if r not in methods:
            raise QError('reference member %s is missing in the run (members: %s)' % (r, ' '.join(methods)))
    res = {}
    for m in methods:
        members = [k for k in refs if k != m]
        rc = [p[m] - p['cpu'] for p in L]
        re_ = [p[m] - sum(p[k] for k in members) / len(members) for p in L]
        res[m] = dict(n=len(rc), rc=rc, re=re_)
    return res, n_exp


def summ(v):
    b = boot_means(v)
    return dict(n=len(v), mean=statistics.fmean(v), lo=q(b, 0.025), hi=q(b, 0.975), sd=statistics.pstdev(v), boot=b)


def main():
    args = sys.argv[1:]
    allow_v1 = False
    refarg = None
    rest = []
    i = 0
    while i < len(args):
        if args[i] == '--allow-v1':
            allow_v1 = True
        elif args[i] == '--ref':
            if i + 1 >= len(args):
                raise QError('--ref needs a value')
            refarg = args[i + 1]; i += 1
        else:
            rest.append(args[i])
        i += 1
    if refarg is None or len(rest) < 3 or rest[1] not in ('pTT', 'pmax', 'emax', 'eF', 'eTT'):
        raise QError('usage: t207_d5_agg.py OUT.json pTT|pmax|emax|eF|eTT --ref def|log|LIST [--allow-v1] ROLE:LABEL=path/qref.txt ...')
    outf, which = rest[0], rest[1]
    refs = ref_list(refarg)
    runs = []
    seen = set()
    for a in rest[2:]:
        m = re.match(r'^(main|ctl):([A-Za-z0-9_.-]+)=(.+)$', a)
        if not m:
            raise QError('bad run argument %r (ROLE:LABEL=path)' % a)
        if m.group(2) in seen:
            raise QError('duplicate label %s' % m.group(2))
        seen.add(m.group(2))
        runs.append((m.group(1), m.group(2), m.group(3), parse(m.group(3), allow_v1)))
    table = {}
    pooled = {}
    inputs = {}
    for role, lab, fn, run in runs:
        res, n_exp = per_run(run, which, refs)
        inputs[lab] = dict(role=role, path=fn, qref_sha256=QP.sha256_file(fn), n_points=len(run.points), n_used=res['cpu']['n'], n_expected=n_exp,
                           variant=run.spec['run']['variant'], args=run.spec['args'], env=run.spec['run']['env'], members=sorted(res))
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
    corr_all = set(refs)
    pool = {m: dict(S_cpu=summ(d['rc']), S_ens=summ(d['re'])) for m, d in pooled.items()}
    # print
    labs = [lab for _, lab, _, _ in runs]
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
    json.dump(dict(measure=which, refs=refs, inputs=inputs, runs=table, pooled=pool), open(outf, 'w'), indent=1)


if __name__ == '__main__':
    try:
        main()
    except QError as e:
        sys.stderr.write('t207_d5_agg ERROR: %s\n' % e)
        sys.exit(1)
