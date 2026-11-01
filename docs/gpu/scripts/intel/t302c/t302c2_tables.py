#!/usr/bin/env python3
"""Ticket 302c-2 (development only, dev seeds): calibration tables of the device tau pass on the B580, generated from the 207-2 aggregate.

1. For every measure (pTT, emax, eF) and reference list, t207_d5_agg.py (commit 5d876291, unchanged, next to this script) is run on all
   given runs; its JSON and text output are kept (OUTDIR/agg_<tag>_<measure>_<ref>.{json,txt}).  These are the per-set and pooled tables.
2. From the same code (agg.per_run / agg.summ on points read by the strict reader t207_d5_qparse.load) three views agg does not print:
   - flavour subsets: S_ens over the points with nfe = 1 resp. nfe = 2 (per set and pooled);
   - leave-one-set-out: S_ens pooled over all main runs but one (for each one);
   - the separation summary per measure: device max over sets, highest correct member per set, controls (var8, ctl_dev) min.
   The statistic is agg's: per point r_m = log10 e_m - mean_{k in refs, k != m} log10 e_k; S_ens = mean over points (pTT: points with B).
usage: t302c2_tables.py [--prog-sha SHA256] OUTDIR TAG REFNAME=LIST [REFNAME=LIST ...] -- ROLE:LABEL=RUNDIR/qref.txt ...
   REFNAME=LIST: a name for the list and an explicit comma list (agg --ref syntax; 'def' / 'log' are agg's frozen NVIDIA lists).
   The first list is the primary one (the separation summary and the flavour / leave-one-out views use every list).
"""
import importlib.util, json, os, statistics, subprocess, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    sp = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + '.py'))
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


AG = load('t207_d5_agg')
QP = AG.QP
MEAS = ('pTT', 'emax', 'eF')
PY = sys.executable


def fmt(s):
    return '%+.3f [%+.3f,%+.3f] n%d' % (s['mean'], s['lo'], s['hi'], s['n'])


def sub_run(run, nfe):
    pts = [p for p in run.points if nfe is None or p.hd['nfe'] == nfe]
    return QP.Run(run.spec, pts, run.dumps)


def main():
    a = sys.argv[1:]
    psha = []
    if a[:1] == ['--prog-sha']:
        psha = ['--prog-sha', a[1]]          # forwarded to agg: every run's recorded binary128-program sha256 must equal it
        a = a[2:]
    if '--' not in a or len(a) < 5:
        sys.exit(__doc__)
    i = a.index('--')
    outdir, tag, refspecs, runargs = a[0], a[1], a[2:i], a[i + 1:]
    refs = []
    for r in refspecs:
        k, v = r.split('=', 1)
        refs.append((k, AG.ref_list(v), v))
    os.makedirs(outdir, exist_ok=True)
    runs = []
    for ra in runargs:
        role, rest = ra.split(':', 1)
        lab, fn = rest.split('=', 1)
        run = QP.load(fn)
        if QP.env_role(run.spec['run']['env'], run.spec['run']['variant']) != role:
            raise QP.QError('run %s: role %s does not match its env' % (lab, role))
        runs.append((role, lab, fn, run))
    result = dict(tag=tag, refs={k: l for k, l, _ in refs}, agg={}, flavour={}, loso={}, separation={})
    txt = []
    # 1. agg itself
    for meas in MEAS:
        for rname, rl, rv in refs:
            base = os.path.join(outdir, 'agg_%s_%s_%s' % (tag, meas, rname))
            p = subprocess.run([PY, '-I', os.path.join(HERE, 't207_d5_agg.py'), base + '.json', meas, '--ref', rv] + psha + runargs,
                               capture_output=True, text=True)
            open(base + '.txt', 'w').write(p.stdout + p.stderr)
            if p.returncode != 0:
                raise QP.QError('agg failed for %s %s: %s' % (meas, rname, p.stderr.strip()))
            result['agg']['%s/%s' % (meas, rname)] = os.path.basename(base) + '.json'
    # 2. views from agg's functions
    mains = [(lab, run) for role, lab, _, run in runs if role == 'main']
    ctls = [(lab, run) for role, lab, _, run in runs if role == 'ctl']
    for meas in MEAS:
        for rname, rl, _ in refs:
            key = '%s/%s' % (meas, rname)
            # flavour subsets
            fl = {}
            for nfe in (1, 2):
                per, pool = {}, {}
                for lab, run in mains + ctls:
                    res, _ = AG.per_run(sub_run(run, nfe), meas, rl)
                    isctl = (lab, run) in ctls
                    for m, d in res.items():
                        if isctl and m != 'dev':
                            continue
                        name = 'ctl_dev' if isctl else m
                        per.setdefault(lab, {})[name] = AG.summ(d['re'])
                        pool.setdefault(name, []).extend(d['re'])
                fl[nfe] = dict(per={l: {m: {k: v for k, v in s.items() if k != 'boot'} for m, s in d.items()} for l, d in per.items()},
                               pooled={m: {k: v for k, v in AG.summ(v).items() if k != 'boot'} for m, v in pool.items()})
            result['flavour'][key] = fl
            # leave-one-set-out over the main runs
            full = {lab: AG.per_run(run, meas, rl)[0] for lab, run in mains}
            lo = {}
            for out in ([l for l, _ in mains] if len(mains) > 1 else []):
                pool = {}
                for lab, res in full.items():
                    if lab == out:
                        continue
                    for m, d in res.items():
                        pool.setdefault(m, []).extend(d['re'])
                lo[out] = {m: {k: v for k, v in AG.summ(v).items() if k != 'boot'} for m, v in pool.items()}
            result['loso'][key] = lo
            # separation summary (per set means)
            corr = [m for m in rl if m != 'cpu'] + ['cpu']
            sep = {}
            for lab, res in full.items():
                sm = {m: statistics.fmean(d['re']) for m, d in res.items()}
                cm = max(corr, key=lambda m: sm[m])
                sep[lab] = dict(dev=sm['dev'], var8=sm.get('var8'), max_correct=sm[cm], max_correct_member=cm,
                                u_max=max(sm.get('u1', -9), sm.get('u2', -9)))
            for lab, run in ctls:
                res, _ = AG.per_run(run, meas, rl)
                sep[lab] = dict(ctl_dev=statistics.fmean(res['dev']['re']))
            result['separation'][key] = sep
    json.dump(result, open(os.path.join(outdir, 'tables_%s.json' % tag), 'w'), indent=1)
    # text
    labs = [l for l, _ in mains] + [l for l, _ in ctls]
    for meas in MEAS:
        for rname, rl, _ in refs:
            key = '%s/%s' % (meas, rname)
            txt.append('== %s  measure %s  ref %s = %s' % (tag, meas, rname, ','.join(rl)))
            txt.append('  separation (S_ens per set): ' + '; '.join(
                ('%s dev %+.3f var8 %+.3f max-correct %+.3f (%s)' % (l, s['dev'], s['var8'] if s['var8'] is not None else float('nan'),
                                                                     s['max_correct'], s['max_correct_member']))
                if 'dev' in s else '%s ctl_dev %+.3f' % (l, s['ctl_dev']) for l, s in result['separation'][key].items()))
            for nfe in (1, 2):
                f = result['flavour'][key][nfe]
                txt.append('  flavour %d, S_ens per set and pooled:' % nfe)
                names = ['dev', 'var8', 'ctl_dev', 'cpu', 'nvh'] + sorted(m for m in f['pooled'] if m.startswith('n') and m not in ('nvh', 'n0')) + ['u1', 'u2']
                for m in names:
                    if m not in f['pooled']:
                        continue
                    row = '    %-8s' % m + ''.join('%30s' % (fmt(f['per'][l][m]) if m in f['per'].get(l, {}) else '-') for l in labs)
                    txt.append(row + '   pooled ' + fmt(f['pooled'][m]))
            txt.append('  leave-one-set-out pooled S_ens (dev / var8 / highest correct member):')
            for out, d in result['loso'][key].items():
                cm = max((m for m in rl), key=lambda m: d[m]['mean'])
                txt.append('    without %-8s dev %s  var8 %s  %s %s' % (out, fmt(d['dev']), fmt(d['var8']) if 'var8' in d else '-', cm, fmt(d[cm])))
    open(os.path.join(outdir, 'tables_%s.txt' % tag), 'w').write('\n'.join(txt) + '\n')
    print('\n'.join(txt))


if __name__ == '__main__':
    try:
        main()
    except QP.QError as e:
        sys.stderr.write('t302c2_tables ERROR: %s\n' % e)
        sys.exit(1)
