#!/usr/bin/env python3
"""Ticket 207-4: copy of t207_3_cal.py (207-3, 428576dd), used at beta = 10 with T207_3_QPDIR = g5 (approved reader 4118a1de); one
addition: OUT.json also holds the per-point values (res['points'], see main()); tables and text output unchanged.
Ticket 207-3 (development only, dev seeds): calibration tables of the two accuracy gates of decision 33 on stride-3 QREF runs.

Reads every run through the strict reader t207_d5_qparse.py (load: every point, member and number checked; unchanged 207-2 code) and adds
what the 207-2 aggregate does not cover (to be merged into the pipeline later): the control kind 'flt' (graded output fault
ALF_GPU_FAULT_TAU=eps:<x>, ticket 207-3 task 1), the pointwise measures emax / eF next to pTT, flavour subsets, and pair pooling.

usage: t207_3_cal.py OUT.json PATH(def|log) ROLE:LABEL=RUNDIR ...
  ROLE attr  env ALF_GPU_{DEBUG=1, TAUDIAG_QHOST=1, TAUDIAG=2, TAUDIAG_QVAR=<masks>, TAUDIAG_QREF, TAUDIAG_FROM} (task 3, host-twin masks)
       main  env exactly ALF_GPU_{DEBUG=1, TAUDIAG_QHOST=1, TAUDIAG_QVAR=8, TAUDIAG_QREF, TAUDIAG_FROM}
       ctl   the same without QVAR, plus TAU_NS=0 (def) or TAUSWAP=8 (log)          (one env form per path: without QVAR)
       flt   the same without QVAR, plus FAULT_TAU=eps:<x> [and FAULT_TAU_AT=<nst>:<nfe>]
  Every run: qref_proc.json rc 0, check_rc 0, cover 1, path as given, null codes contain the frozen reference list; qref.txt loads
  strictly.  Reference ensemble: the frozen list of the path (t207_d5_qparse.REF_DEF / REF_LOG).
Statistic, per measure M in (pTT, emax, eF) and point p:  r_m(p) = log10 M_m(p) - mean_{k in refs, k != m} log10 M_k(p)  (leave-one-out
for a member); S_ens(m) = mean over points (pTT: points with B; emax, eF: all points).  Flavour subsets nfe = 1 / 2 of the same points.
Printed (and in OUT.json): per run and pooled over the main runs, S_ens of dev, var8, every reference member (correct variants),
the controls' dev; pair pooling (every pair of main runs, as a gate on two seed sets would pool); bootstrap 95 % intervals (4000 draws).
"""
import importlib.util, itertools, json, math, os, random, re, statistics, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
QPDIR = os.environ.get('T207_3_QPDIR', HERE)      # directory of the strict reader to use (e.g. the 305488ec copy in g3/)
_spec = importlib.util.spec_from_file_location('t207_d5_qparse', os.path.join(QPDIR, 't207_d5_qparse.py'))
QP = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(QP)
QError = QP.QError
MEASURES = ('pTT', 'emax', 'eF')
CTL = {'def': ('ALF_GPU_TAU_NS', '0'), 'log': ('ALF_GPU_TAUSWAP', '8')}


def role_of(env, path):
    base = {'ALF_GPU_DEBUG': '1', 'ALF_GPU_TAUDIAG_QHOST': '1'}
    for k in ('ALF_GPU_TAUDIAG_QREF', 'ALF_GPU_TAUDIAG_FROM'):
        if not re.match(r'^[0-9]+$', env.get(k, 'x')):
            raise QError('%s missing or not a number' % k)
        base[k] = env[k]
    if env == dict(base, ALF_GPU_TAUDIAG_QVAR='8'):
        return 'main', None
    if env == dict(base, **{CTL[path][0]: CTL[path][1]}):
        return 'ctl', None
    qv = env.get('ALF_GPU_TAUDIAG_QVAR', '')
    if re.match(r'^[0-9]+(,[0-9]+)*$', qv) and env == dict(base, ALF_GPU_TAUDIAG='2', ALF_GPU_TAUDIAG_QVAR=qv):
        return 'attr', None
    f = env.get('ALF_GPU_FAULT_TAU', '')
    m = re.match(r'^eps:([0-9.eE+-]+)$', f)
    if m and env == dict(base, ALF_GPU_FAULT_TAU=f):
        return 'flt', float(m.group(1))
    at = env.get('ALF_GPU_FAULT_TAU_AT', '')
    if m and re.match(r'^[0-9]+:[0-9]+$', at) and env == dict(base, ALF_GPU_FAULT_TAU=f, ALF_GPU_FAULT_TAU_AT=at):
        return 'flt', '%s@%s' % (m.group(1), at)
    raise QError('env %s is none of main / ctl / flt for path %s' % (' '.join('%s=%s' % kv for kv in sorted(env.items())), path))


def meas(e, which):
    if which == 'emax':
        return max(e[0:4])
    if which == 'eF':
        return e[4]
    return e[5]


def boot(v, n=4000, seed=20261101):
    rng = random.Random(seed)
    L = len(v)
    b = sorted(sum(v[rng.randrange(L)] for _ in range(L)) / L for _ in range(n))
    return b


def summ(v):
    b = boot(v)
    return dict(n=len(v), mean=statistics.fmean(v), lo=b[int(0.025 * len(b))], hi=b[int(0.975 * len(b)) - 1], sd=statistics.pstdev(v))


def load(rundir, path):
    q = os.path.join(rundir, 'qref.txt')
    run = QP.load(q)
    with open(os.path.join(rundir, 'qref_proc.json')) as f:
        pr = json.load(f)
    a = run.spec['args']
    if pr.get('rc') != 0 or pr.get('check_rc') != 0 or a['cover'] != 1 or a['path'] != path or run.spec['run']['variant'] != path:
        raise QError('%s: processing not passing / cover != 1 / wrong path' % rundir)
    refs = (QP.REF_DEF if path == 'def' else QP.REF_LOG).split(',')
    refs = ['nvh' if r == 'H' else r for r in refs]      # the approved reader (4118a1de) writes the host member as 'H'
    for r in refs:
        if r not in ('cpu', 'nvh') and int(r[1:]) not in a['nul']:
            raise QError('%s: reference %s not processed' % (rundir, r))
    role, x = role_of(run.spec['run']['env'], path)
    return run, refs, role, x, pr.get('prog')


def ratios(run, refs, which, nfe=None):
    """per method: list of r_m(p)"""
    out = {}
    for pt in run.points:
        if which == 'pTT' and not pt.hasb:
            continue
        if nfe is not None and pt.hd['nfe'] != nfe:
            continue
        lg = {}
        for k, e in pt.m.items():
            if k in QP.DIAG_NAMES:
                continue
            v = meas(e, which)
            if not (math.isfinite(v) and v > 0):
                raise QError('dump %d %s: %s = %r' % (pt.hd['file'], k, which, v))
            lg[k] = math.log10(v)
        for m in lg:
            mem = [k for k in refs if k != m]
            out.setdefault(m, []).append(lg[m] - sum(lg[k] for k in mem) / len(mem))
    return out


def ratios_pos(run, refs, which):
    """per method: list of (position key (nst, nfe, br), r_m(p))"""
    out = {}
    for pt in run.points:
        if which == 'pTT' and not pt.hasb:
            continue
        key = (pt.hd['nst'], pt.hd['nfe'], pt.hd['br'])
        lg = {k: math.log10(meas(e, which)) for k, e in pt.m.items() if k not in QP.DIAG_NAMES}
        for m in lg:
            mem = [k for k in refs if k != m]
            out.setdefault(m, []).append((key, lg[m] - sum(lg[k] for k in mem) / len(mem)))
    return out


def q(v, p):
    v = sorted(v)
    k = (len(v) - 1) * p
    i = int(math.floor(k)); j = min(i + 1, len(v) - 1)
    return v[i] + (v[j] - v[i]) * (k - i)


def strata(runs, refs, which, minpts=8):
    """per-stratum means (strata (nst, nfe) and (nst, nfe, br), pooled over the runs of one role) and tail statistics (q90, q95 of
    the per-point values) for dev, var8, every reference member (leave-one-out) and each control run's dev"""
    res = {}
    mains = [r for r in runs if r['role'] == 'main']
    groups = [('main', mains)] + ([('run:' + r['label'], [r]) for r in mains] if len(mains) > 1 else []) + \
             [(r['label'], [r]) for r in runs if r['role'] in ('ctl', 'flt')]
    for gname, rs in groups:
        if not rs:
            continue
        vals = {}
        for r in rs:
            for m, lst in ratios_pos(r['run'], refs, which).items():
                vals.setdefault(m, []).extend(lst)
        methods = [m for m in vals if m in refs or m in ('dev', 'var8')] if gname == 'main' or gname.startswith('run:') else ['dev']
        g = {}
        for m in methods:
            lst = vals[m]
            d = dict(q90=q([v for _, v in lst], 0.90), q95=q([v for _, v in lst], 0.95), max=max(v for _, v in lst), n=len(lst))
            per_run_q = {}
            for lev, kf in (('nf', lambda k: (k[0], k[1])), ('nfb', lambda k: k)):
                st = {}
                for k, v in lst:
                    st.setdefault(kf(k), []).append(v)
                d[lev] = {'%s' % ('.'.join(str(x) for x in k)): (statistics.fmean(v), len(v)) for k, v in sorted(st.items())}
                ok = [mv for mv, n in d[lev].values() if n >= minpts]
                d[lev + '_max'] = max(ok) if ok else None
            g[m] = d
        res[gname] = g
    return res


def main():
    a = sys.argv[1:]
    if len(a) < 3 or a[1] not in ('def', 'log'):
        sys.exit('usage: t207_3_cal.py OUT.json def|log ROLE:LABEL=RUNDIR ...')
    outf, path = a[0], a[1]
    runs = []
    for arg in a[2:]:
        m = re.match(r'^(main|ctl|flt|attr):([A-Za-z0-9_.-]+)=(.+)$', arg)
        if not m:
            raise QError('bad argument %r' % arg)
        run, refs, role, x, prog = load(m.group(3), path)
        if role != m.group(1):
            raise QError('%s: given role %s, env says %s' % (m.group(2), m.group(1), role))
        runs.append(dict(label=m.group(2), dir=m.group(3), role=role, x=x, run=run, refs=refs, prog=prog,
                         seeds=run.spec['run']['seeds_sha256'], npts=len(run.points), nb=run.n_with_b))
    progs = set(r['prog'] for r in runs)
    if len(progs) != 1:
        raise QError('runs processed by different programs: %s' % progs)
    refs = runs[0]['refs']
    res = dict(path=path, refs=refs, prog=progs.pop(), runs=[{k: v for k, v in r.items() if k != 'run'} for r in runs], tables={})
    for which in MEASURES:
        for nfe in (None, 1, 2):
            key = '%s%s' % (which, '' if nfe is None else '_nfe%d' % nfe)
            per = {}
            for r in runs:
                rr = ratios(r['run'], refs, which, nfe)
                if r['role'] in ('main', 'attr'):
                    per[r['label']] = {m: summ(v) for m, v in rr.items()}
                    per[r['label']]['_raw'] = rr
                else:
                    name = 'ctl_dev' if r['role'] == 'ctl' else 'flt_dev(eps:%s)' % r['x']
                    per[r['label']] = {name: summ(rr['dev']), '_raw': {name: rr['dev']}}
            mains = [r['label'] for r in runs if r['role'] == 'main']
            attrs = [r['label'] for r in runs if r['role'] == 'attr']
            pooled = {}
            for m in per[mains[0]]['_raw'] if mains else []:
                pooled[m] = summ(sum((per[l]['_raw'][m] for l in mains), []))
            pairs = {}
            for l1, l2 in itertools.combinations(mains, 2):
                pairs['%s+%s' % (l1, l2)] = {m: statistics.fmean(per[l1]['_raw'][m] + per[l2]['_raw'][m]) for m in per[l1]['_raw']}
            for l in per:
                per[l].pop('_raw')
            for l in attrs:
                per[l] = {m: per[l][m] for m in per[l]}
            res['tables'][key] = dict(per_run=per, pooled=pooled, pairs=pairs)
    res['strata'] = {w: strata(runs, refs, w) for w in MEASURES}
    # 207-4 addition (the only change against t207_3_cal.py): the per-point values r_m(p) with their position (nst, nfe, br), per measure,
    # run and member (main runs: dev, var8 and every reference member; control runs: dev), so that t207_4_sep.py can evaluate every
    # pooling (one run, a pair of runs, all runs) and statistic (mean, q95, max, by_stratum, per flavour) from the same numbers
    res['points'] = {w: {r['label']: {m: [[k[0], k[1], k[2], v] for k, v in lst] for m, lst in ratios_pos(r['run'], refs, w).items()
                                      if (m == 'dev' or (r['role'] == 'main' and (m in refs or m == 'var8')))} for r in runs} for w in MEASURES}
    json.dump(res, open(outf, 'w'), indent=1)
    # ---- text report (generated from res)
    print('207-3 calibration, path %s, program %s, refs %s' % (path, res['prog'][:16], ','.join(refs)))
    for r in res['runs']:
        print('  run %-22s role %-4s x %-8s seeds %s points %d (with B %d)  %s' % (r['label'], r['role'], r['x'], r['seeds'][:16], r['npts'], r['nb'], r['dir']))
    for key, t in res['tables'].items():
        mains = [r['label'] for r in res['runs'] if r['role'] == 'main']
        print('== S_ens(%s)  mean [95%% bootstrap]' % key)
        print('  %-18s' % 'method' + ''.join('%24s' % l for l in mains) + '%24s' % 'pooled(main)')
        corr = [m for m in refs]
        for m in ['dev', 'var8'] + corr:
            row = '  %-18s' % (m + (' (corr)' if m in corr else ''))
            for l in mains:
                s = t['per_run'][l].get(m)
                row += '%24s' % ('%+.3f [%+.3f,%+.3f]' % (s['mean'], s['lo'], s['hi']) if s else '-')
            s = t['pooled'].get(m)
            row += '%24s' % ('%+.3f [%+.3f,%+.3f]' % (s['mean'], s['lo'], s['hi']) if s else '-')
            print(row)
        if mains:
            cm = {l: max(t['per_run'][l][m]['mean'] for m in corr) for l in mains}
            print('  max correct per run: ' + ', '.join('%s %+.3f (%s)' % (l, v, max(corr, key=lambda m: t['per_run'][l][m]['mean'])) for l, v in cm.items())
                  + '; pooled max %+.3f' % max(t['pooled'][m]['mean'] for m in corr))
            if t['pairs']:
                pc = max(v[m] for v in t['pairs'].values() for m in corr)
                pd = [v['dev'] for v in t['pairs'].values()]
                p8 = [v['var8'] for v in t['pairs'].values()]
                print('  pairs of main runs (%d): correct max %+.3f, dev %+.3f..%+.3f, var8 %+.3f..%+.3f' % (len(t['pairs']), pc, min(pd), max(pd), min(p8), max(p8)))
        for r in res['runs']:
            if r['role'] == 'attr':
                t1 = t['per_run'][r['label']]
                print('  attribution run %s: ' % r['label'] + ', '.join('%s %+.3f [%+.3f,%+.3f]' % (m, t1[m]['mean'], t1[m]['lo'], t1[m]['hi'])
                                                                     for m in ['dev'] + sorted((m for m in t1 if m.startswith('var')), key=lambda z: int(z[3:])) + ['nvh', 'cpu']))
            elif r['role'] != 'main':
                for name, s in t['per_run'][r['label']].items():
                    print('  control %-22s %-18s %+.3f [%+.3f,%+.3f]  n %d' % (r['label'], name, s['mean'], s['lo'], s['hi'], s['n']))


def print_strata(res):
    refs = res['refs']
    for w, st in res['strata'].items():
        print('== strata and tail, measure %s (per-stratum means over the pooled points of the role; strata with < 8 points not in max)' % w)
        for gname, g in st.items():
            if gname.startswith('run:'):
                continue
            for m, d in g.items():
                tag = gname if gname != 'main' else ('main ' + m + (' (corr)' if m in refs else ''))
                print('  %-34s q90 %+.3f q95 %+.3f max %+.3f | max over (nst,nfe) %+.3f | max over (nst,nfe,br) %+.3f | (nst,nfe): %s' % (
                    tag, d['q90'], d['q95'], d['max'], d['nf_max'], d['nfb_max'] if d['nfb_max'] is not None else float('nan'),
                    ' '.join('%s:%+.2f/%d' % (k, v[0], v[1]) for k, v in d['nf'].items())))
        for gname in [g for g in st if g.startswith('run:')]:
            corr = [m for m in st[gname] if m in refs]
            print('  each run %-20s correct max: point %+.3f q95 %+.3f stratum(nst,nfe) %+.3f stratum(nst,nfe,br) %+.3f | dev: q95 %+.3f max %+.3f '
                  'stratum %+.3f / %+.3f | var8: q95 %+.3f stratum %+.3f' % (gname[4:], max(st[gname][m]['max'] for m in corr),
                  max(st[gname][m]['q95'] for m in corr), max(st[gname][m]['nf_max'] for m in corr), max(st[gname][m]['nfb_max'] for m in corr),
                  st[gname]['dev']['q95'], st[gname]['dev']['max'], st[gname]['dev']['nf_max'], st[gname]['dev']['nfb_max'],
                  st[gname]['var8']['q95'], st[gname]['var8']['nf_max']))
        if 'main' in st:
            corr = [m for m in st['main'] if m in refs]
            print('  correct variants: max point value %+.3f, max q90 %+.3f, max q95 %+.3f, max stratum (nst,nfe) %+.3f, max stratum (nst,nfe,br) %+.3f' % (
                max(st['main'][m]['max'] for m in corr), max(st['main'][m]['q90'] for m in corr), max(st['main'][m]['q95'] for m in corr),
                max(st['main'][m]['nf_max'] for m in corr), max(st['main'][m]['nfb_max'] for m in corr)))


if __name__ == '__main__':
    try:
        main()
        if os.environ.get('T207_3_STRATA') == '1':
            print_strata(json.load(open(sys.argv[1])))
    except QError as e:
        sys.stderr.write('t207_3_cal ERROR: %s\n' % e)
        sys.exit(1)
