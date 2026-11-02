#!/usr/bin/env python3
"""Ticket 207-3 (development only, dev seeds): summary table of the PROPOSED accuracy gates, generated from the calibration JSON files
of t207_3_cal.py (one per path).  Nothing here is a pre-registration: the bounds below are the proposal the manager brings to the owner
(fixed only in 207-6).

usage: t207_3_summary.py OUT.json CAL_def.json CAL_log.json

Gates (statistic of the per-point values r(p) of t207_3_cal.py; definitions as in t207_d5_gate.py of gpu-port-207-2 4118a1de):
  decisive (proposed):  S_ens(pTT) pooled <= 0.15;  S_ens(emax) pooled <= 0.15;  S_ens(emax) by_stratum <= 0.30;  S_ens_q95(emax) pooled <= 0.60
  candidate (t207-3c):  S_ens(pTT) flavour max <= 0.15 = max(S_ens(pTT) of the flavour-1 points, of the flavour-2 points)  (decision 35)
  report only:          S_ens(pTT) by_stratum, S_ens_q95(pTT) pooled, S_ens_max(pTT) pooled, S_ens_max(emax) pooled, S_ens(emax) flavour max
  pooled     : over all points of the runs of the group
  by_stratum : largest per-stratum mean, stratum = (nst, nfe, br), pooled over the runs of the group; minimum stratum size rule: the one
               of the cal JSON (min_stratum_pts; t207-3c: 1 = every stratum, the rule of the approved gate); the smallest stratum seen
               is reported
  q95        : 95th percentile of the per-point values (linear interpolation)
Groups: main = the main runs pooled (4 per path in the 207-3 data); T2 pooling (t207-3c): ASSUMED that T2 pools 1 or 2 main runs per path
(one or two acceptance seed sets), so every statistic is also given for each main run alone ('run:') and every pair of main runs
('pair:'), and the correct-variant maximum is the largest value over all of them (and over the frozen reference members, leave-one-out);
a control (var8 of the main runs aside) is a single run.  var8 (QVAR=8 stand-in, a member of the main runs) is given at the smallest value
over the same runs and pairs (a control must fail in every pooling).
Margin: bound - value (device, correct variants), value - bound (controls; > 0 = killed).  Separation (decision 37) is stated for the
S_ens(pTT) gates only, over every pTT control (var8, the production control, every graded fault): every correct variant passes in every
T2 pooling and every control fails.  The near-threshold rule of the 207-3 review (margin >= 0.10 above the bound) is flagged per control.
Probes: graded-fault runs measured as response probes, known to lie below the detection threshold and not candidate controls (the 207-3
closing round: LOG eps 2e-15, +0.032), are shown in every table but not counted in the separation; they are listed in PROBES below with
their reason and printed with the separation lines.
"""
import json, statistics, sys

GATES = [  # name, measure, statistic, pool, bound, kind
    ('pTT_pooled', 'pTT', 'S_ens', 'pooled', 0.15, 'decisive'),
    ('pTT_flavour_max', 'pTT', 'S_ens_flmax', 'pooled', 0.15, 'candidate'),
    ('emax_pooled', 'emax', 'S_ens', 'pooled', 0.15, 'decisive'),
    ('emax_by_stratum', 'emax', 'S_ens', 'by_stratum', 0.30, 'decisive'),
    ('emax_q95', 'emax', 'S_ens_q95', 'pooled', 0.60, 'decisive'),
    ('pTT_by_stratum', 'pTT', 'S_ens', 'by_stratum', None, 'report'),
    ('pTT_q95', 'pTT', 'S_ens_q95', 'pooled', None, 'report'),
    ('pTT_max', 'pTT', 'S_ens_max', 'pooled', None, 'report'),
    ('emax_max', 'emax', 'S_ens_max', 'pooled', None, 'report'),
    ('emax_flavour_max', 'emax', 'S_ens_flmax', 'pooled', None, 'report'),
]
NEAR_MARGIN = 0.10
PROBES = {'f_log_2em15_dev9': 'LOG eps 2e-15: response probe below the pTT detection threshold (207-3 closing round, +0.032); not a candidate control'}


def sens(cal, group, member, measure, nfe=None):
    """S_ens (mean) of one member in one group from the tables of the cal JSON"""
    t = cal['tables'][measure + ('' if nfe is None else '_nfe%d' % nfe)]
    if group == 'main':
        return t['pooled'][member]['mean']
    if group.startswith('run:'):
        return t['per_run'][group[4:]][member]['mean']
    if group.startswith('pair:'):
        return t['pairs'][group[5:]][member]
    (name, s), = t['per_run'][group].items()          # a control run: its only entry is its dev
    return s['mean']


def value(cal, group, member, measure, stat, pool):
    """(value, detail) of one member of one group ('main', 'run:<l>', 'pair:<l1>+<l2>' or a control run label)"""
    if stat == 'S_ens_flmax':
        v1, v2 = sens(cal, group, member, measure, 1), sens(cal, group, member, measure, 2)
        return max(v1, v2), 'nfe1 %+.3f nfe2 %+.3f' % (v1, v2)
    st = cal['strata'][measure][group][member]
    if pool == 'by_stratum':
        mp = cal['min_stratum_pts']
        ok = {k: v for k, v in st['nfb'].items() if v[1] >= mp}
        k, (v, n) = max(ok.items(), key=lambda kv: kv[1][0])
        return v, 'at %s (%d pts; strata %d, pts %d..%d)' % (k, n, len(st['nfb']), min(x[1] for x in st['nfb'].values()),
                                                              max(x[1] for x in st['nfb'].values()))
    if stat == 'S_ens_q95':
        return st['q95'], 'n %d' % st['n']
    if stat == 'S_ens_max':
        return st['max'], 'n %d' % st['n']
    if group in ('main',) or group[:4] == 'run:':
        t = cal['tables'][measure]
        s = t['pooled'][member] if group == 'main' else t['per_run'][group[4:]][member]
        return s['mean'], '[%+.3f,%+.3f] n %d' % (s['lo'], s['hi'], s['n'])
    if group[:5] == 'pair:':
        return sens(cal, group, member, measure), 'pair'
    (name, s), = cal['tables'][measure]['per_run'][group].items()
    return s['mean'], '[%+.3f,%+.3f] n %d' % (s['lo'], s['hi'], s['n'])


def ctl_label(r):
    return r['label'] + ('' if r['role'] == 'ctl' else ' (eps:%s)' % r['x'])


def t2_groups(cal):
    return sorted(g for g in cal['strata']['pTT'] if g[:4] == 'run:' or g[:5] == 'pair:')


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    outf = sys.argv[1]
    res = dict(gates=[dict(name=g[0], measure=g[1], statistic=g[2], pool=g[3], bound=g[4], kind=g[5]) for g in GATES], paths={},
               t2_pooling_assumed='T2 pools 1 or 2 main runs per path: every main run alone and every pair of main runs')
    lines, flav, sep, small, efl, spr = [], [], [], [], [], []
    mps = set()
    for path, f in (('def', sys.argv[2]), ('log', sys.argv[3])):
        cal = json.load(open(f))
        if cal['path'] != path:
            sys.exit('%s: path %s, expected %s' % (f, cal['path'], path))
        mps.add(cal['min_stratum_pts'])
        refs = cal['refs']
        mains = [r['label'] for r in cal['runs'] if r['role'] == 'main']
        ctls = [r for r in cal['runs'] if r['role'] in ('ctl', 'flt')]
        T2 = t2_groups(cal)
        lines.append('== path %s: %s; main runs %s; program %s; reader %s' % (path, f, ','.join(mains), cal['prog'][:16], cal.get('reader', '?')))
        lines.append('   T2 pooling groups (%d): %s' % (len(T2), ' '.join(T2)))
        P = res['paths'][path] = dict(cal=f, mains=mains, refs=refs, rows=[])
        # smallest stratum seen (every measure, every group, every member, both stratum levels)
        for w, st in cal['strata'].items():
            for gname, g in st.items():
                for m, d in g.items():
                    for lev in ('nf', 'nfb'):
                        small.append((min(v[1] for v in d[lev].values()), path, w, gname, m, lev))
        for name, meas, stat, pool, bound, kind in GATES:
            row = dict(gate=name, bound=bound, kind=kind)
            dv, dd = value(cal, 'main', 'dev', meas, stat, pool)
            cv = [(value(cal, 'main', m, meas, stat, pool), m) for m in refs]
            (cmax, cdet), cm = max(cv, key=lambda z: z[0][0])
            # T2 pooling: each main run and each pair of main runs
            t2c = [(value(cal, g, m, meas, stat, pool)[0], g, m) for g in T2 for m in refs]
            t2d = [(value(cal, g, 'dev', meas, stat, pool)[0], g) for g in T2]
            t2v = [(value(cal, g, 'var8', meas, stat, pool)[0], g) for g in T2]
            run_c = max(z for z in t2c if z[1][:4] == 'run:') if any(z[1][:4] == 'run:' for z in t2c) else None
            pair_c = max(z for z in t2c if z[1][:5] == 'pair:') if any(z[1][:5] == 'pair:' for z in t2c) else None
            row.update(device=dv, device_detail=dd, correct_max=cmax, correct_max_member=cm, correct_detail=cdet,
                       t2_correct_max_run=run_c, t2_correct_max_pair=pair_c, t2_device_range=[min(t2d)[0], max(t2d)[0]] if t2d else None,
                       t2_var8_min=min(t2v) if t2v else None)
            cts = [('var8 (QVAR=8, main runs pooled)',) + value(cal, 'main', 'var8', meas, stat, pool)]
            if t2v:
                v8, g8 = min(t2v)
                cts.append(('var8, smallest T2 pooling (%s)' % g8, v8, 'min over %d runs/pairs' % len(t2v)))
            cts += [(ctl_label(r),) + value(cal, r['label'], 'dev', meas, stat, pool) for r in ctls]
            row['controls'] = [dict(name=n, value=v, detail=d, margin=(v - bound if bound is not None else None),
                                    killed=(v > bound if bound is not None else None)) for n, v, d in cts]
            P['rows'].append(row)
            b = '%.2f' % bound if bound is not None else 'report'
            hdr = '  %-16s %-9s bound %-6s' % (name, kind, b)
            if bound is not None:
                lines.append(hdr + ' device %+.3f (margin %+.3f) %s | correct max %+.3f %s (margin %+.3f) %s' % (
                    dv, bound - dv, dd, cmax, cm, bound - cmax, cdet))
            else:
                lines.append(hdr + ' device %+.3f %s | correct max %+.3f %s %s' % (dv, dd, cmax, cm, cdet))
            if t2d:
                s = '  %-34s T2 pooling: device %+.3f..%+.3f' % ('', min(t2d)[0], max(t2d)[0])
                for lab, z in (('single run', run_c), ('pair', pair_c)):
                    if z:
                        s += ' | correct max %s %+.3f (%s, %s)%s' % (lab, z[0], z[1], z[2], '' if bound is None else ' margin %+.3f' % (bound - z[0]))
                s += ' | var8 min %+.3f (%s)' % min(t2v)
                lines.append(s)
            for c in row['controls']:
                if bound is not None:
                    lines.append('  %-34s control %-44s %+.3f  margin %+.3f  %s%s  %s' % (
                        '', c['name'], c['value'], c['margin'], 'KILLED' if c['killed'] else 'NOT killed',
                        ' (margin < %.2f)' % NEAR_MARGIN if c['killed'] and c['margin'] < NEAR_MARGIN and meas == 'pTT' else '', c['detail']))
                else:
                    lines.append('  %-34s control %-44s %+.3f  %s' % ('', c['name'], c['value'], c['detail']))
            # separation by decision 37 (pTT S_ens gates): correct variants in every T2 pooling vs every pTT control
            if meas == 'pTT' and stat in ('S_ens', 'S_ens_flmax') and pool == 'pooled':
                cw = max([z for z in (run_c, pair_c) if z] + [(cmax, 'main', cm)])
                cc = [c for c in row['controls'] if c['name'].split(' ')[0] not in PROBES]
                probes = [c['name'] for c in row['controls'] if c['name'].split(' ')[0] in PROBES]
                cl = min(cc, key=lambda c: c['value'])
                notk = [c['name'] for c in cc if not c['killed']]
                near = [c['name'] for c in cc if c['killed'] and c['margin'] < NEAR_MARGIN]
                ok = cw[0] < bound and not notk
                sep.append('  %s %-16s bound %.2f: correct max over main/runs/pairs %+.3f (%s, %s; margin %+.3f); smallest control %+.3f (%s; '
                           'margin %+.3f); not killed: %s; killed with margin < %.2f: %s; probes not counted: %s -> %s' % (
                               path, name, bound, cw[0], cw[1], cw[2], bound - cw[0], cl['value'], cl['name'], cl['margin'],
                               ', '.join(notk) or 'none', NEAR_MARGIN, ', '.join(near) or 'none', ', '.join(probes) or 'none',
                               'SEPARATES' if ok else 'DOES NOT SEPARATE'))
                row['separation'] = dict(correct_worst=cw, smallest_control=cl, not_killed=notk, near=near, separates=ok)
        # flavour 1 vs flavour 2 (task 0): S_ens per flavour, pooled main runs for device / correct / var8, each control run
        for meas in ('pTT', 'emax', 'eF'):
            def fl(g, m):
                return sens(cal, g, m, meas, 1), sens(cal, g, m, meas, 2)
            rows = [('device', fl('main', 'dev'))]
            cmx = max(refs, key=lambda m: max(fl('main', m)))
            rows.append(('correct max (%s)' % cmx, fl('main', cmx)))
            d12 = [fl('main', m)[1] - fl('main', m)[0] for m in refs]
            rows.append(('correct: nfe2 - nfe1 range', (min(d12), max(d12))))
            rows.append(('var8', fl('main', 'var8')))
            rows += [(ctl_label(r), fl(r['label'], 'dev')) for r in ctls]
            for lab, (a, b2) in rows:
                if lab.startswith('correct: '):
                    flav.append('  %s %-5s %-46s nfe2-nfe1 from %+.3f to %+.3f' % (path, meas, lab, a, b2))
                else:
                    flav.append('  %s %-5s %-46s nfe1 %+.3f  nfe2 %+.3f  nfe2-nfe1 %+.3f' % (path, meas, lab, a, b2, b2 - a))
        # spread of the graded-fault runs: S_ens(pTT) pooled and per flavour, grouped by fault (x and position); the run-to-run sd of the
        # device and of the correct variants over the main runs (pooled pTT) for comparison
        byx = {}
        for r in ctls:
            if r['role'] == 'flt':
                byx.setdefault(str(r['x']), []).append(r['label'])
        dsd = statistics.stdev([sens(cal, 'run:' + l, 'dev', 'pTT') for l in mains]) if len(mains) > 1 else float('nan')
        csd = max(statistics.stdev([sens(cal, 'run:' + l, m, 'pTT') for l in mains]) for m in refs) if len(mains) > 1 else float('nan')
        spr.append('  %s run-to-run sd of S_ens(pTT) over the %d main runs: device %.3f, correct variants up to %.3f' % (path, len(mains), dsd, csd))
        for x, labs in sorted(byx.items()):
            v = [sens(cal, l, 'dev', 'pTT') for l in labs]
            f2 = [(sens(cal, l, 'dev', 'pTT', 1), sens(cal, l, 'dev', 'pTT', 2)) for l in labs]
            spr.append('  %s eps:%-12s runs %d: %s | mean %+.3f sd %s range %.3f' % (
                path, x, len(labs), ' '.join('%s %+.3f (nfe1 %+.3f nfe2 %+.3f)' % (l, a, b1, b2) for l, a, (b1, b2) in zip(labs, v, f2)),
                statistics.fmean(v), '%.3f' % statistics.stdev(v) if len(v) > 1 else '-', max(v) - min(v)))
        # eF vs emax (why eF is not gated): pooled S_ens of the device, the correct maximum and every control, both measures
        kill = {}
        for meas in ('emax', 'eF'):
            cm = max(refs, key=lambda m: sens(cal, 'main', m, meas))
            kill[meas] = [ctl_label(r) for r in ctls if sens(cal, r['label'], 'dev', meas) > 0.15]
            efl.append('  %s %-4s device %+.3f correct max %+.3f (%s) | ' % (path, meas, sens(cal, 'main', 'dev', meas), sens(cal, 'main', cm, meas), cm)
                       + ' '.join('%s %+.3f' % (ctl_label(r), sens(cal, r['label'], 'dev', meas)) for r in ctls))
        efl.append('  %s controls above 0.15: emax %s | eF %s -> %s' % (path, ', '.join(kill['emax']) or 'none', ', '.join(kill['eF']) or 'none',
                                                                   'the same set' if kill['emax'] == kill['eF'] else 'DIFFERENT sets'))
    if len(mps) != 1:
        sys.exit('the two cal files use different minimum stratum size rules: %s' % mps)
    mp = mps.pop()
    sm = min(small)
    res.update(min_stratum_pts=mp, smallest_stratum=dict(points=sm[0], path=sm[1], measure=sm[2], group=sm[3], member=sm[4], level=sm[5]),
               flavour_table=flav, separation=sep, eF_vs_emax=efl, fault_spread=spr)
    json.dump(res, open(outf, 'w'), indent=1)
    print('207-3 proposed accuracy gates (development, dev seeds; generated by t207_3_summary.py from the cal JSON files)')
    print('minimum stratum size rule: a stratum enters by_stratum if it has >= %d point(s) (the approved gate\'s rule: every stratum); '
          'smallest stratum seen: %d points (path %s, measure %s, group %s, member %s, level %s)' % ((mp,) + sm))
    print('T2 pooling assumed: 1 or 2 main runs per path (each main run alone and every pair of main runs)')
    print('\n'.join(lines))
    print('== separation of the S_ens(pTT) gates by decision 37 (every correct variant passes in every pooling, every pTT control fails)')
    for k, v in sorted(PROBES.items()):
        print('  probe %s: %s' % (k, v))
    print('\n'.join(sep))
    print('== spread of the graded-fault runs (S_ens(pTT) pooled per run, and per flavour)')
    print('\n'.join(spr))
    print('== flavour 1 vs flavour 2 (S_ens per flavour; device / correct / var8: main runs pooled; controls: their run)')
    print('\n'.join(flav))
    print('== why eF is not gated: if eF puts the same controls above 0.15 as emax (pooled S_ens, last line per path), it adds no '
          'detection; emax also carries the by_stratum gate')
    print('\n'.join(efl))


if __name__ == '__main__':
    main()
