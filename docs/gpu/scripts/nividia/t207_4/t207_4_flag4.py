#!/usr/bin/env python3
"""Ticket 207-4 task 3 (development only): do the row-permuted offline null variants (flag 4 of the code n<c>: c mod 100 has bit 4, i.e.
'row-permuted QR input' in nul_cgr22) sit at the top of the correct band?  Generated from t207_4_cal.py / t207_3_cal.py JSON files.

usage: t207_4_flag4.py LABEL=CAL.json [LABEL=CAL.json ...]     (e.g. b10_log=logs/report_cal_log.json b4_log=../t207-3/logs/report_cal_log.json)
Per file and measure (pTT, emax): S_ens of every correct member (leave-one-out, reference list of the file) per main run and pooled, the
members split into flag-4 (n4, n100004, n5, n100005) and the others; the mean of each group; the rank of the flag-4 members.
"""
import json, sys


def flag4(m):
    return m.startswith('n') and m[1:].isdigit() and (int(m[1:]) % 100) & 4 != 0


def main():
    for arg in sys.argv[1:]:
        label, f = arg.split('=', 1)
        res = json.load(open(f))
        refs = res['refs']
        mains = [r['label'] for r in res['runs'] if r['role'] == 'main']
        f4 = [m for m in refs if flag4(m)]
        ot = [m for m in refs if not flag4(m)]
        print('=== %s (%s, path %s): flag-4 members %s; others %s' % (label, f, res['path'], ','.join(f4), ','.join(ot)))
        for meas in ('pTT', 'emax'):
            t = res['tables'][meas]
            print('  -- S_ens(%s)' % meas)
            print('  %-10s' % 'member' + ''.join('%14s' % l for l in mains) + '%14s' % 'pooled')
            for m in f4 + ot:
                print('  %-10s' % (m + (' *' if m in f4 else '')) + ''.join('%+14.3f' % t['per_run'][l][m]['mean'] for l in mains)
                      + '%+14.3f' % t['pooled'][m]['mean'])
            for name, grp in (('mean flag-4', f4), ('mean others', ot)):
                print('  %-10s' % name[:10] + ''.join('%+14.3f' % (sum(t['per_run'][l][m]['mean'] for m in grp) / len(grp)) for l in mains)
                      + '%+14.3f' % (sum(t['pooled'][m]['mean'] for m in grp) / len(grp)))
            top = []
            for l in mains + ['pooled']:
                v = t['pooled'] if l == 'pooled' else t['per_run'][l]
                order = sorted(refs, key=lambda m: -v[m]['mean'])
                top.append('%s: top-%d = %s' % (l, len(f4), ','.join(order[:len(f4)])))
            print('  ranks: ' + '; '.join(top))


if __name__ == '__main__':
    main()
