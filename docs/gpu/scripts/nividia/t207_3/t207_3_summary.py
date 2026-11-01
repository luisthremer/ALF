#!/usr/bin/env python3
"""Ticket 207-3 (development only, dev seeds): summary table of the PROPOSED accuracy gates, generated from the calibration JSON files
of t207_3_cal.py (one per path).  Nothing here is a pre-registration: the bounds below are the proposal the manager brings to the owner
(fixed only in 207-6).

usage: t207_3_summary.py OUT.json CAL_def.json CAL_log.json

Gates (statistic of the per-point values r(p) of t207_3_cal.py; definitions as in t207_d5_gate.py of gpu-port-207-2 4118a1de):
  decisive (proposed):  S_ens(pTT) pooled <= 0.15;  S_ens(emax) pooled <= 0.15;  S_ens(emax) by_stratum <= 0.30;  S_ens_q95(emax) pooled <= 0.60
  report only:          S_ens(pTT) by_stratum, S_ens_q95(pTT) pooled, S_ens_max(pTT) pooled, S_ens_max(emax) pooled
  pooled     : over all points of the runs of the role (device and correct variants: the main runs pooled; a control: its run)
  by_stratum : largest per-stratum mean, stratum = (nst, nfe, br), every stratum (no minimum point count, as the gate), pooled over the
               runs of the role
  q95        : 95th percentile of the per-point values (linear interpolation), pooled
Rows per path: device (dev of the main runs), the correct-variant maximum (largest value over the frozen reference members, leave-one-out,
pooled main runs; also the largest single-run value for the pooled S_ens), and every control: var8 (QVAR=8 stand-in in the main runs),
the production control (ctl) and every graded fault (flt).  Margin: bound - value (device, correct variants), value - bound (controls;
> 0 = killed).
"""
import json, sys

GATES = [  # name, measure, statistic, pool, bound, decisive
    ('pTT_pooled', 'pTT', 'S_ens', 'pooled', 0.15, True),
    ('emax_pooled', 'emax', 'S_ens', 'pooled', 0.15, True),
    ('emax_by_stratum', 'emax', 'S_ens', 'by_stratum', 0.30, True),
    ('emax_q95', 'emax', 'S_ens_q95', 'pooled', 0.60, True),
    ('pTT_by_stratum', 'pTT', 'S_ens', 'by_stratum', None, False),
    ('pTT_q95', 'pTT', 'S_ens_q95', 'pooled', None, False),
    ('pTT_max', 'pTT', 'S_ens_max', 'pooled', None, False),
    ('emax_max', 'emax', 'S_ens_max', 'pooled', None, False),
]


def value(cal, group, member, measure, stat, pool):
    """(value, detail) of one member of one group ('main' or a control run label) from the cal JSON"""
    st = cal['strata'][measure][group][member]
    if pool == 'by_stratum':
        k, (v, n) = max(st['nfb'].items(), key=lambda kv: kv[1][0])
        return v, 'at %s (%d pts; strata %d, pts %d..%d)' % (k, n, len(st['nfb']), min(x[1] for x in st['nfb'].values()),
                                                              max(x[1] for x in st['nfb'].values()))
    if stat == 'S_ens_q95':
        return st['q95'], 'n %d' % st['n']
    if stat == 'S_ens_max':
        return st['max'], 'n %d' % st['n']
    t = cal['tables'][measure]
    if group == 'main':
        s = t['pooled'][member]
    else:
        (name, s), = t['per_run'][group].items()
    return s['mean'], '[%+.3f,%+.3f] n %d' % (s['lo'], s['hi'], s['n'])


def ctl_label(r):
    return r['label'] + ('' if r['role'] == 'ctl' else ' (eps:%s)' % r['x'])


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    outf = sys.argv[1]
    res = dict(gates=[dict(name=g[0], measure=g[1], statistic=g[2], pool=g[3], bound=g[4], decisive=g[5]) for g in GATES], paths={})
    lines = []
    for path, f in (('def', sys.argv[2]), ('log', sys.argv[3])):
        cal = json.load(open(f))
        if cal['path'] != path:
            sys.exit('%s: path %s, expected %s' % (f, cal['path'], path))
        refs = cal['refs']
        mains = [r['label'] for r in cal['runs'] if r['role'] == 'main']
        ctls = [r for r in cal['runs'] if r['role'] in ('ctl', 'flt')]
        lines.append('== path %s: %s; main runs %s; program %s; reader %s' % (path, f, ','.join(mains), cal['prog'][:16], cal.get('reader', '?')))
        P = res['paths'][path] = dict(cal=f, mains=mains, refs=refs, rows=[])
        for name, meas, stat, pool, bound, dec in GATES:
            row = dict(gate=name, bound=bound, decisive=dec)
            dv, dd = value(cal, 'main', 'dev', meas, stat, pool)
            cv = [(value(cal, 'main', m, meas, stat, pool), m) for m in refs]
            (cmax, cdet), cm = max(cv, key=lambda z: z[0][0])
            row.update(device=dv, device_detail=dd, correct_max=cmax, correct_max_member=cm, correct_detail=cdet)
            if stat == 'S_ens' and pool == 'pooled':
                pr = [(cal['tables'][meas]['per_run'][l][m]['mean'], l, m) for l in mains for m in refs]
                row['correct_max_single_run'] = max(pr)
            cts = [('var8 (QVAR=8, main runs)',) + value(cal, 'main', 'var8', meas, stat, pool)]
            cts += [(ctl_label(r),) + value(cal, r['label'], 'dev', meas, stat, pool) for r in ctls]
            row['controls'] = [dict(name=n, value=v, detail=d, margin=(v - bound if bound is not None else None),
                                    killed=(v > bound if bound is not None else None)) for n, v, d in cts]
            P['rows'].append(row)
            b = '%.2f' % bound if bound is not None else 'report'
            hdr = '  %-16s bound %-6s' % (name, b)
            if bound is not None:
                lines.append(hdr + ' device %+.3f (margin %+.3f) %s | correct max %+.3f %s (margin %+.3f) %s' % (
                    dv, bound - dv, dd, cmax, cm, bound - cmax, cdet))
            else:
                lines.append(hdr + ' device %+.3f %s | correct max %+.3f %s %s' % (dv, dd, cmax, cm, cdet))
            if 'correct_max_single_run' in row:
                v, l, m = row['correct_max_single_run']
                lines.append('  %-30s correct max single run %+.3f (%s, %s)' % ('', v, l, m))
            for c in row['controls']:
                if bound is not None:
                    lines.append('  %-30s control %-34s %+.3f  margin %+.3f  %s  %s' % ('', c['name'], c['value'], c['margin'],
                                                                                     'KILLED' if c['killed'] else 'NOT killed', c['detail']))
                else:
                    lines.append('  %-30s control %-34s %+.3f  %s' % ('', c['name'], c['value'], c['detail']))
    json.dump(res, open(outf, 'w'), indent=1)
    print('207-3 proposed accuracy gates (development, dev seeds; generated by t207_3_summary.py from the cal JSON files)')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
