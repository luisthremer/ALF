#!/usr/bin/env python3
"""Ticket 207-4 (development only, dev seeds): separation tables of the candidate accuracy gates at beta = 10, generated from the per-point
values that t207_4_cal.py stores in its JSON (res['points']).  No number is typed in.

usage: t207_4_sep.py CAL.json [CAL.json ...]      (one per path)

Statistics (gate semantics of t207_d5_gate.py, 4118a1de), on the per-point log ratios r(p) of a member over a pool of runs:
  S_ens = mean; q95 = 95th percentile (linear interpolation, as t207_3_cal.q); max = maximum; by_stratum = the largest stratum mean,
  stratum = (nst, nfe, br), EVERY stratum counts (no minimum size, as in the gate); the smallest stratum size is printed.
Measures: pTT (points with B), pTT_f1 / pTT_f2 (flavour 1 / 2 only), pTT_fmax (the larger of the two flavour means), emax, eF (all points).
Poolings of the main runs (correct variants, dev, var8): 'all' = all main runs pooled; 'pair' = the largest value over every pair of main
  runs (a T2 gate on two seed sets); 'run' = the largest value over the single runs.  Control runs: one run each (their dev).
For each candidate: dev, correct max (over the frozen reference members, leave-one-out, with member and pool), var8, every control,
margins against the beta = 4 bound proposed by 207-3 where one exists (bound - value for dev/correct, value - bound for defects), and the
range decision 37 allows per pooling: max(correct, dev) <= b < min over the defects above it (defects at or below max(correct, dev) are not
catchable by that gate at any bound).  A control killed by less than 0.10 over the beta = 4 bound is flagged NEAR-THIN (207-3 review).
Margin rule (manager, from the 207-3c/3d/3e review, going into 207-6; as corrected 2026-11-02): a correct variant passes with margin when
its worst single-run margin (bound - its largest single-run value) is >= 0.03 AND z = (bound - its mean over the main runs) / sd >= 2.5,
sd = the run-to-run sd (n-1) of its single-run values over the main runs; a control is killed with margin when value - bound >= 0.10 on
each of >= 2 fresh seed sets.  Printed per gate:
  L = max over the correct members of max(largest single-run value + 0.03, mean + 2.5 sd): the smallest bound every correct variant passes
      with margin (the member that sets it, and the dev's largest single-run value for comparison);
  per control: the largest bound that kills it with margin, value - 0.10 (here ONE seed set per control: the '>= 2 seed sets' part cannot
      be met by these data), and whether L <= value - 0.10 (catchable with margin by some bound of this gate);
  for the beta = 4 bound: whether it satisfies the rule on the correct side (b4 >= L) and which controls it kills with margin.
"""
import itertools, json, math, statistics, sys

B4 = {('pTT', 'S_ens'): 0.15, ('emax', 'S_ens'): 0.15, ('emax', 'by_stratum'): 0.30, ('emax', 'q95'): 0.60, ('pTT_fmax', 'S_ens'): 0.15}
B4_PATH = {('def', 'pTT_fmax', 'S_ens'): 0.20}     # 207-3d: def pTT_flavour_max proposed at 0.20 (LOG 0.15)
CANDIDATES = [('pTT', 'S_ens'), ('pTT_f1', 'S_ens'), ('pTT_f2', 'S_ens'), ('pTT_fmax', 'S_ens'), ('pTT', 'by_stratum'), ('pTT', 'q95'),
              ('pTT', 'max'), ('emax', 'S_ens'), ('emax', 'by_stratum'), ('emax', 'q95'), ('emax', 'max'), ('eF', 'S_ens'),
              ('eF', 'by_stratum'), ('eF', 'q95')]
NEAR = 0.10


def q(v, p):
    v = sorted(v)
    k = (len(v) - 1) * p
    i = int(math.floor(k)); j = min(i + 1, len(v) - 1)
    return v[i] + (v[j] - v[i]) * (k - i)


def pts(res, meas, labels, member):
    base = 'pTT' if meas.startswith('pTT') else meas
    out = []
    for l in labels:
        out.extend(res['points'][base][l][member])
    if meas in ('pTT_f1', 'pTT_f2'):
        out = [p for p in out if p[1] == int(meas[-1])]
    return out


def stat(res, meas, st, labels, member):
    """(value, detail) of statistic st of member pooled over the runs labels"""
    if meas == 'pTT_fmax':
        a = stat(res, 'pTT_f1', st, labels, member)[0]; b = stat(res, 'pTT_f2', st, labels, member)[0]
        return (a, 'f1') if a >= b else (b, 'f2')
    p = pts(res, meas, labels, member)
    v = [x[3] for x in p]
    if st == 'S_ens':
        return sum(v) / len(v), ''
    if st == 'q95':
        return q(v, 0.95), ''
    if st == 'max':
        return max(v), ''
    sd = {}
    for x in p:
        sd.setdefault((x[0], x[1], x[2]), []).append(x[3])
    k, lst = max(sd.items(), key=lambda kv: sum(kv[1]) / len(kv[1]))
    return sum(lst) / len(lst), 'stratum %d.%d.%d/%d of %d strata, smallest %d' % (k + (len(lst), len(sd), min(len(x) for x in sd.values())))


def main():
    for f in sys.argv[1:]:
        res = json.load(open(f))
        refs = res['refs']
        mains = [r['label'] for r in res['runs'] if r['role'] == 'main']
        ctls = [r['label'] for r in res['runs'] if r['role'] in ('ctl', 'flt')]
        pools = [('all', [mains])] + ([('pair', [list(c) for c in itertools.combinations(mains, 2)])] if len(mains) > 2 else []) + \
                ([('run', [[l] for l in mains])] if len(mains) > 1 else [])
        print('=== path %s: %s; main runs %s; controls %s; program %s' % (res['path'], f, ', '.join(mains), ', '.join(ctls), res['prog'][:16]))
        for meas in ('pTT', 'emax'):
            row = []
            for l in mains + ctls:
                sd = {}
                for x in res['points'][meas][l]['dev']:
                    sd[(x[0], x[1], x[2])] = sd.get((x[0], x[1], x[2]), 0) + 1
                row.append('%s %d strata, sizes %d..%d' % (l, len(sd), min(sd.values()), max(sd.values())))
            print('  strata (nst,nfe,br) of %s per run: %s' % (meas, '; '.join(row)))
        caught = {}
        rulecaught = {}
        for meas, st in CANDIDATES:
            b4 = B4_PATH.get((res['path'], meas, st), B4.get((meas, st)))
            print('-- %s %s   (beta = 4 bound of 207-3: %s)' % (st, meas, '%+.2f' % b4 if b4 is not None else 'none'))
            mg = lambda v, d: '' if b4 is None else ' (margin %+.3f)' % ((v - b4) if d else (b4 - v))
            ctlv = [(stat(res, meas, st, [c], 'dev'), c) for c in ctls]
            for pname, groups in pools:
                def worst(members):
                    best = None
                    for g in groups:
                        for m in members:
                            v, det = stat(res, meas, st, g, m)
                            if best is None or v > best[0]:
                                best = (v, m, '+'.join(g) if pname != 'all' else 'all', det)
                    return best
                dv, cv, v8 = worst(['dev']), worst(refs), worst(['var8'])
                print('   pool %-4s dev %+.3f%s [%s%s] | correct max %+.3f%s [%s, %s%s] | var8 %+.3f%s [%s]' % (
                    pname, dv[0], mg(dv[0], False), dv[2], ', ' + dv[3] if dv[3] else '', cv[0], mg(cv[0], False), cv[1], cv[2],
                    ', ' + cv[3] if cv[3] else '', v8[0], mg(v8[0], True), v8[2]))
                lo = max(dv[0], cv[0])
                planted = [(v8[0], 'var8(%s)' % pname)] + [(v, c) for (v, _), c in ctlv]
                above = sorted((v, c) for v, c in planted if v > lo)
                below = [(v, c) for v, c in planted if v <= lo]
                print('        range (decision 37): %+.3f <= b < %s; not catchable: %s' % (
                    lo, '%+.3f (%s)' % above[0] if above else 'none', ', '.join('%s %+.3f' % (c, v) for v, c in below) or '-'))
                if b4 is not None and pname == 'all':
                    for v, c in planted:
                        caught.setdefault(c, []).append(('%s %s' % (st, meas), v > b4))
            if len(mains) > 1:
                lm = []
                for m in refs:
                    rv = [stat(res, meas, st, [l], m)[0] for l in mains]
                    sdm = statistics.stdev(rv)
                    lm.append((max(max(rv) + 0.03, statistics.fmean(rv) + 2.5 * sdm), m, max(rv), sdm, statistics.fmean(rv)))
                L = max(lm)
                dmax = max(stat(res, meas, st, [l], 'dev')[0] for l in mains)
                print('   margin rule: L = %+.3f (set by %s: largest single run %+.3f, mean %+.3f, run-to-run sd %.3f); dev largest single run %+.3f' % (L[0], L[1], L[2], L[4], L[3], dmax))
                cm = [(v - NEAR, c) for (v, _), c in ctlv]
                print('        controls, largest bound killing with margin (value - 0.10; one seed set each): ' +
                      ', '.join('%s %+.3f%s' % (c, u, '' if u >= L[0] else ' (< L: not catchable with margin)') for u, c in cm))
                ok = [(u, c) for u, c in cm if u >= L[0]]
                print('        bounds meeting the rule for the catchable controls: %s' % (
                    '[%+.3f, %+.3f]' % (L[0], min(ok)[0]) if ok else 'NONE (no control is catchable with margin)'))
                if b4 is not None:
                    print('        beta = 4 bound %+.2f: correct side %s; kills with margin: %s' % (
                        b4, 'meets the rule' if b4 >= L[0] else 'FAILS the rule (b4 < L by %.3f)' % (L[0] - b4),
                        ', '.join(c for (v, _), c in ctlv if v - b4 >= NEAR) or 'none'))
                for (v, _), c in ctlv:
                    rulecaught.setdefault(c, []).append(('%s %s' % (st, meas), v - NEAR >= L[0]))
            for (v, det), c in ctlv:
                flag = ''
                if b4 is not None and 0 < v - b4 < NEAR:
                    flag = '  NEAR-THIN: killed by %+.3f < %.2f' % (v - b4, NEAR)
                print('   control %-28s %+.3f%s%s%s' % (c, v, mg(v, True), ' [' + det + ']' if det else '', flag))
        print('-- planted defects vs the gates with a beta = 4 bound (%s), pool all: caught by' % ', '.join('%s %s %+.2f' % (k[1], k[0], b) for k, b in B4.items()))
        for c, lst in caught.items():
            g = [n for n, ok in lst if ok]
            print('   %-28s %s' % (c, ', '.join(g) if g else 'NONE'))
        print('-- margin rule: candidate gates that could catch each control with margin at some bound >= L (single seed set)')
        for c, lst in rulecaught.items():
            g = [n for n, ok in lst if ok]
            print('   %-28s %s' % (c, ', '.join(g) if g else 'NONE'))


if __name__ == '__main__':
    main()
