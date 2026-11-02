#!/usr/bin/env python3
"""Ticket 207-4 (development only): measured response of the gate statistics to the graded sweep output fault ALF_GPU_FAULT_TAU=eps:<x>
(global, not localized) and the x that would put a control a given margin above a bound.  Generated from t207_4_cal.py JSON files.

usage: t207_4_xneed.py BOUND_pTT BOUND_emax MARGIN CAL.json [CAL.json ...]
Model (an interpolation, not a law): the fault adds an error proportional to x, so S(x) = S_dev + log10(1 + x/x0); x0 is solved from each
measured graded run (S(x_i), S_dev = the main runs' pooled dev), and x_need = x0 (10^(bound + margin - S_dev) - 1) from the measured run
whose S is nearest to the target.  Printed per path and measure (pTT, emax; S_ens pooled): S_dev, each graded run (x, S, x0), x_need.
"""
import json, math, sys


def mean(v):
    return sum(v) / len(v)


def main():
    bp, be, mg = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    for f in sys.argv[4:]:
        res = json.load(open(f))
        mains = [r['label'] for r in res['runs'] if r['role'] == 'main']
        flt = [r for r in res['runs'] if r['role'] == 'flt' and isinstance(r['x'], (int, float))]
        print('=== path %s (%s)' % (res['path'], f))
        for meas, b in (('pTT', bp), ('emax', be)):
            sdev = mean([p[3] for l in mains for p in res['points'][meas][l]['dev']])
            tgt = b + mg
            rows = []
            for r in sorted(flt, key=lambda r: r['x']):
                s = mean([p[3] for p in res['points'][meas][r['label']]['dev']])
                d = s - sdev
                x0 = r['x'] / (10 ** d - 1) if d > 0 else float('inf')
                rows.append((r['x'], s, x0, r['label']))
            print('  %s: S_dev(main pooled) %+.3f, target bound %+.2f + margin %.2f = %+.3f' % (meas, sdev, b, mg, tgt))
            for x, s, x0, l in rows:
                print('    x %-8g S %+.3f  x0 %.3g  (%s)' % (x, s, x0, l))
            fin = [r for r in rows if math.isfinite(r[2])]
            if fin and tgt > sdev:
                x, s, x0, l = min(fin, key=lambda r: abs(r[1] - tgt))
                print('    x_need %.3g (from %s, x0 %.3g)' % (x0 * (10 ** (tgt - sdev) - 1), l, x0))


if __name__ == '__main__':
    main()
