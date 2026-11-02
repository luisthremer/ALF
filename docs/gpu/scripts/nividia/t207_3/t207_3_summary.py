#!/usr/bin/env python3
"""Ticket 207-3 (development only, dev seeds): summary table of the PROPOSED accuracy gates, generated from the calibration JSON files
of t207_3_cal.py (one per path).  Nothing here is a pre-registration: the bounds below are the proposal the manager brings to the owner
(fixed only in 207-6).

usage: t207_3_summary.py OUT.json CAL_def.json CAL_log.json

Gates (statistic of the per-point values r(p) of t207_3_cal.py; definitions as in t207_d5_gate.py of gpu-port-207-2 4118a1de):
  decisive (proposed):  S_ens(pTT) pooled <= 0.15;  S_ens(emax) pooled <= 0.15;  S_ens(emax) by_stratum <= 0.30;  S_ens_q95(emax) pooled <= 0.60
  per-path bounds (BOUND_PATH, 207-3d): def pTT_flavour_max <= 0.20 (LOG keeps the 0.15 of the table)
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
Separation rule (207-3e, as written in the ticket review): a correct variant passes WITH MARGIN when its single-run margin is >= 0.03 and
z = (bound - mean over the main runs) / sd >= 2.5 (sd = run-to-run sd of that member over the main runs; this z reproduces the review's
2.7 / 3.5 / 3.3 / 1.9 of 207-3d, the margin of the worst single run is not divided by sd); a control is killed WITH MARGIN when margin >= 0.10.  Per gate:
SEPARATES (correct side with margin, every control killed with margin) / THIN (everything passes or is killed, but a margin rule fails) /
DOES NOT SEPARATE (a correct variant fails the bound or a control is not killed).  Gate-set claims per path (GATE_SETS) are computed
here, not assembled by hand.  PROPOSED marks the controls proposed for pre-registration (a column; the counted runs are unchanged).
T207_3_TEST_BOUND='def:pTT_flavour_max=0.15' overrides one bound for a test of the rule (printed as TEST ONLY).
Margin: bound - value (device, correct variants), value - bound (controls; > 0 = killed).  Separation (decision 37; 207-3f: for every
gated statistic, pTT and emax): every correct variant (reference members AND the device) passes in every T2 pooling.  The pTT gates judge
every counted control (var8, the production control, every graded fault) and every control must fail; the emax gates judge the PROPOSED
controls only and list the other counted controls (killed / not killed) as information.  The near-threshold rule of the 207-3 review (margin >= 0.10 above the bound) is flagged per control.
Probes: graded-fault runs measured as response probes, known to lie below the detection threshold and not candidate controls (the 207-3
closing round: LOG eps 2e-15, +0.032), are shown in every table but not counted in the separation; they are listed in PROBES below with
their reason and printed with the separation lines.
"""
import json, math, os, re, statistics, sys

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
# per-path bound overrides (207-3d): (path, gate name) -> bound; replaces the bound of GATES for that path only
BOUND_PATH = {('def', 'pTT_flavour_max'): 0.20}
NEAR_MARGIN = 0.10
# 207-3e/3f: the written separation rule (207-6 "Margin rule").  A correct variant (every frozen reference member and the device, each over
# the main runs) passes WITH MARGIN when its worst single-run margin (bound - largest main-run value) is >= CORR_MARGIN and
# z = (bound - mean over the main runs) / sd >= CORR_Z, sd = run-to-run sd of that member over the main runs (the single-run margin is NOT
# divided by sd).  A control is killed WITH MARGIN when value - bound >= CTL_MARGIN on each of >= CTL_SETS fresh seed sets (kinds of
# PROPOSED controls, enforced below); a control present in every main run (var8) may instead show z_fix = (smallest main-run value - bound)
# / run-to-run sd over the main runs >= FIX_Z.  The correct side also includes the device: above the bound in any run, pair or the main
# pool it FAILS (DOES NOT SEPARATE).  The rule is generated for every gated statistic (pTT and emax).
CORR_MARGIN, CORR_Z, CTL_MARGIN = 0.03, 2.5, NEAR_MARGIN
CTL_SETS, FIX_Z = 2, 4.0
# gate sets proposed for 207-6 (S_ens(pTT) gates only): the claim computed below is that every control is killed with margin by >= 1 gate
GATE_SETS = {'def': ['pTT_pooled', 'pTT_flavour_max'], 'log': ['pTT_pooled', 'pTT_flavour_max']}
# controls proposed for pre-registration (label regexes; var8 only against pTT_pooled); a marking column, not a filter of the counted runs
PROPOSED = {'def': [r'f_def_4em15_dev\d+', r'f_def_4em15_at0x2_dev\d+'], 'log': [r'f_log_1em14_dev\d+', r'f_log_6em15_at0x2_dev\d+']}
# test-only bound override for checking that the rule can fail: T207_3_TEST_BOUND='def:pTT_flavour_max=0.15' (printed as TEST ONLY)
TEST_BOUND = os.environ.get('T207_3_TEST_BOUND')
if TEST_BOUND:
    try:
        _k, _v = TEST_BOUND.split('=')
        _p, _g = _k.split(':')
        _v = float(_v)
    except ValueError:
        sys.exit('T207_3_TEST_BOUND=%r: expected path:gate=value' % TEST_BOUND)
    if _p not in GATE_SETS or _g not in [g[0] for g in GATES if g[4] is not None]:
        sys.exit('T207_3_TEST_BOUND=%r: unknown path %r or gate %r (paths %s, gated %s)' % (
            TEST_BOUND, _p, _g, sorted(GATE_SETS), [g[0] for g in GATES if g[4] is not None]))
    BOUND_PATH[(_p, _g)] = _v
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


def ckey(n):
    return re.sub(r' \((run|pair):[^)]*\)$', '', n)      # the smallest-pooling var8 names its group, which differs per gate


def is_proposed(path, name, gate):
    lab = name.split(' ')[0]
    if lab == 'var8':
        return gate == 'pTT_pooled'
    return any(re.fullmatch(q, lab) for q in PROPOSED[path])


def ctl_label(r):
    return r['label'] + ('' if r['role'] == 'ctl' else ' (eps:%s)' % r['x'])


def t2_groups(cal):
    return sorted(g for g in cal['strata']['pTT'] if g[:4] == 'run:' or g[:5] == 'pair:')


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    outf = sys.argv[1]
    res = dict(gates=[dict(name=g[0], measure=g[1], statistic=g[2], pool=g[3], bound=g[4], kind=g[5],
                           bound_by_path={pp: BOUND_PATH.get((pp, g[0]), g[4]) for pp in ('def', 'log')}) for g in GATES],
               rule=dict(correct_margin=CORR_MARGIN, correct_z=CORR_Z, control_margin=CTL_MARGIN), test_only_bound=TEST_BOUND,
               gate_sets=GATE_SETS, proposed_controls=PROPOSED,
               bound_path={'%s:%s' % k: v for k, v in BOUND_PATH.items()}, paths={},
               t2_pooling_assumed='T2 pools 1 or 2 main runs per path: every main run alone and every pair of main runs')
    lines, flav, sep, small, efl, spr, gsl = [], [], [], [], [], [], []
    mps = set()
    for path, f in (('def', sys.argv[2]), ('log', sys.argv[3])):
        cal = json.load(open(f))
        if cal['path'] != path:
            sys.exit('%s: path %s, expected %s' % (f, cal['path'], path))
        mps.add(cal['min_stratum_pts'])
        refs = cal['refs']
        mains = [r['label'] for r in cal['runs'] if r['role'] == 'main']
        if len(mains) < 2:
            sys.exit('%s: %d main run(s); the run-to-run sd and z need >= 2: refusing' % (f, len(mains)))
        mseeds = [r['seeds'] for r in cal['runs'] if r['role'] == 'main']
        if len(set(mains)) != len(mains) or len(set(mseeds)) != len(mseeds):
            sys.exit('%s: duplicated main run (label or seeds sha256 repeated): refusing' % f)
        seedof = {r['label']: r['seeds'] for r in cal['runs']}
        role = {r['label']: r['role'] for r in cal['runs']}
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
            bound = BOUND_PATH.get((path, name), bound)
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
                                    killed=(v > bound if bound is not None else None),
                                    proposed=is_proposed(path, n, name)) for n, v, d in cts]
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
                        ' (margin < %.2f)' % NEAR_MARGIN if c['killed'] and c['margin'] < NEAR_MARGIN else '', c['detail'])
                        + ('  [PREREG-PROPOSED]' if c['proposed'] else ''))
                else:
                    lines.append('  %-34s control %-44s %+.3f  %s' % ('', c['name'], c['value'], c['detail']))
            # separation by decision 37 (pTT S_ens gates): correct variants in every T2 pooling vs every pTT control
            if bound is not None:
                # correct side: reference members AND the device over every run, pair and the main pool (207-3f)
                cw = max([z for z in (run_c, pair_c) if z] + [(cmax, 'main', cm)] + [(v, g, 'dev') for v, g in t2d] + [(dv, 'main', 'dev')])
                cc = [c for c in row['controls'] if c['name'].split(' ')[0] not in PROBES]
                # pTT gates (decision 37): every counted control is judged; emax gates: the PROPOSED controls are judged, the counted controls
                # not killed are listed for information (the emax gates are not claimed to detect every graded fault)
                jc = cc if meas == 'pTT' else [c for c in cc if c['proposed']]
                probes = [c['name'] for c in row['controls'] if c['name'].split(' ')[0] in PROBES]
                cl = min(jc or cc, key=lambda c: c['value'])
                notk = [c['name'] for c in jc if not c['killed']]
                notk_info = [c['name'] for c in cc if not c['killed'] and c not in jc]
                kil_info = ['%s %+.3f' % (c['name'], c['margin']) for c in cc if c['killed'] and c not in jc]
                # written rule (207-3e): correct side per member (reference members and device) over the main runs
                cs = []
                for m in list(refs) + ['dev']:
                    vs = [value(cal, 'run:' + l, m, meas, stat, pool)[0] for l in mains]
                    sd = statistics.stdev(vs)
                    if not (math.isfinite(sd) and sd > 0):
                        sys.exit('%s %s %s: run-to-run sd %r of %s over the main runs is not finite and > 0: refusing' % (path, name, m, sd, vs))
                    mg = bound - max(vs)                       # single-run margin: the worst main run
                    mm = bound - statistics.fmean(vs)          # margin of the member's mean over the main runs
                    cs.append(dict(member=m, worst_run=max(vs), margin=mg, mean_margin=mm, sd=sd, z=mm / sd))
                zw = min(cs, key=lambda c: c['z']); mw = min(cs, key=lambda c: c['margin'])
                corr_pass = cw[0] < bound
                corr_ok = corr_pass and all(c['margin'] >= CORR_MARGIN and c['z'] >= CORR_Z for c in cs)
                corr_side = 'PASSES WITH MARGIN' if corr_ok else ('THIN' if corr_pass else 'FAILS')
                v8 = [value(cal, 'run:' + l, 'var8', meas, stat, pool)[0] for l in mains]
                v8sd = statistics.stdev(v8)
                if not (math.isfinite(v8sd) and v8sd > 0):
                    sys.exit('%s %s var8: run-to-run sd %r over the main runs is not finite and > 0: refusing' % (path, name, v8sd))
                v8z = (min(v8) - bound) / v8sd              # z_fix of the 207-6 rule: worst single main run vs the bound
                near = [c['name'] for c in jc if c['killed'] and c['margin'] < CTL_MARGIN]
                if not corr_pass or notk:
                    ok = 'DOES NOT SEPARATE'
                elif corr_ok and not near:
                    ok = 'SEPARATES'
                else:
                    ok = 'THIN'
                sep.append('  %s %-16s [%s] bound %.2f%s: correct side %s (rule: single-run margin >= %.2f and z >= %.1f; worst margin %+.3f (%s), worst z %.1f '
                           '(%s: mean margin %+.3f, sd %.3f); main-pooled/pair/run max %+.3f (%s, %s)); controls: smallest %+.3f (%s; margin %+.3f), '
                           'not killed: %s%s; killed with margin < %.2f: %s; var8 z_fix = (smallest main run - bound)/sd = %.1f; probes not counted: %s -> %s' % (
                               path, name, meas, bound, ' (TEST ONLY bound)' if TEST_BOUND else '', corr_side, CORR_MARGIN, CORR_Z, mw['margin'], mw['member'],
                               zw['z'], zw['member'], zw['mean_margin'], zw['sd'], cw[0], cw[1], cw[2], cl['value'], cl['name'], cl['margin'],
                               ', '.join(notk) or 'none', (' (judged: %s; counted but not judged, information: not killed: %s; killed: %s)' % (
                                   'every counted control' if meas == 'pTT' else 'PROPOSED controls', ', '.join(notk_info) or 'none',
                                   ', '.join(kil_info) or 'none')),
                               CTL_MARGIN, ', '.join(near) or 'none', v8z, ', '.join(probes) or 'none', ok))
                row['separation'] = dict(correct_worst=cw, smallest_control=cl, not_killed=notk, near=near, correct_side=corr_side,
                                         correct_members=cs, var8_margin_over_sd=v8z, verdict=ok, separates=(ok == 'SEPARATES'))
        # gate-set claim (207-3e): every counted pTT control is killed with margin (>= CTL_MARGIN) by at least one gate of the set
        rowd = {r['gate']: r for r in P['rows']}
        gs = GATE_SETS[path]
        gsl.append('  gate set %s {%s}%s: correct side %s' % (path, ', '.join('%s %.2f' % (g, rowd[g]['bound']) for g in gs),
                                                         ' (TEST ONLY bound)' if TEST_BOUND else '',
                                                         '; '.join('%s %s' % (g, rowd[g]['separation']['correct_side']) for g in gs)))
        # claim B needs, per kind of PROPOSED control, >= CTL_SETS runs on distinct fresh seed sets (distinct from each other, from every
        # main run and from every run of another role); computed from the cal JSON, not assumed (207-3f)
        kinds = {}
        for q in PROPOSED[path]:
            runs = [r['label'] for r in ctls if re.fullmatch(q, r['label'])]
            sh = [seedof[l] for l in runs]
            other = {seedof[l] for l in seedof if role[l] != 'flt'} | {seedof[l] for l in mains}
            fresh = len(set(sh)) == len(sh) and not (set(sh) & other)
            kinds[q] = dict(runs=runs, seeds=[x[:16] for x in sh], distinct_fresh=fresh, enough=(len(runs) >= CTL_SETS and fresh))
            gsl.append('  gate set %s control kind %-28s runs %s seeds %s: %d runs, %s -> %s' % (
                path, q, ','.join(runs) or 'none', ','.join(x[:8] for x in sh) or '-', len(runs),
                'distinct and fresh (not a main-run or other-role seed set)' if fresh else 'seed sets NOT distinct or not fresh',
                'enough (>= %d)' % CTL_SETS if kinds[q]['enough'] else 'TOO FEW fresh seed sets (need >= %d)' % CTL_SETS))
        gcl = []
        for c0 in rowd[gs[0]]['controls']:
            if c0['name'].split(' ')[0] in PROBES:
                continue
            per = []
            for g in gs:
                c = [x for x in rowd[g]['controls'] if ckey(x['name']) == ckey(c0['name'])][0]
                per.append((c['margin'], g, c['proposed']))
            # NOTE (conservative): for var8 each gate's value is the smallest over the main runs / pairs (ckey drops the group), so this is
            # max over gates of (min over groups), which is <= min over groups of (max over gates): the var8 gate-set max-min is conservative
            bm, bg, _ = max(per)
            prop = any(z[2] for z in per)
            # a proposed control counts against the gates it is proposed for (var8: pTT_pooled only)
            pm, pg, _ = max([z for z in per if z[2]] or per)
            # which rule the control meets: margin >= CTL_MARGIN, or (var8, present in every main run) z_fix >= FIX_Z on the gate it is proposed for
            isv8 = c0['name'].split(' ')[0] == 'var8'
            zf = rowd[pg]['separation']['var8_margin_over_sd'] if isv8 else None
            rules = (['margin>=%.2f' % CTL_MARGIN] if pm >= CTL_MARGIN else []) + (['z_fix>=%.0f' % FIX_Z] if isv8 and zf >= FIX_Z and pm > 0 else [])
            gcl.append(dict(control=c0['name'], gate=bg, margin=bm, with_margin=bm >= CTL_MARGIN, proposed=prop,
                            prop_gate=pg, prop_margin=pm, prop_with_margin=pm >= CTL_MARGIN, z_fix=zf, rules_met=rules,
                            prop_ok=bool(rules)))
            gsl.append('    %-58s killed by %-16s margin %+.3f  %s  %s' % (c0['name'], bg if bm > 0 else '-- none --', bm,
                       'WITH MARGIN' if bm >= CTL_MARGIN else ('thin' if bm > 0 else 'NOT KILLED'),
                       ('[PREREG-PROPOSED: gate %s margin %+.3f; rule met: %s%s]' % (
                           pg, pm, ' + '.join(rules) or 'NONE', '' if zf is None else '; z_fix = %.1f' % zf)) if prop else ''))
        na = sum(g['with_margin'] for g in gcl)
        pr = [g for g in gcl if g['proposed']]
        npr = sum(g['prop_ok'] for g in pr)
        prr = [g for g in pr if not g['control'].startswith('var8')]
        nrr = sum(g['prop_with_margin'] for g in prr)
        kinds_ok = all(k['enough'] for k in kinds.values())
        P['gate_set'] = dict(kinds=kinds, kinds_ok=kinds_ok, gates=gs, bounds={g: rowd[g]['bound'] for g in gs}, controls=gcl,
                             all_with_margin=(na == len(gcl)), proposed_with_margin=(npr == len(pr) and kinds_ok),
                             proposed_resizable_with_margin=(nrr == len(prr) and kinds_ok))
        gsl.append('  gate set %s claim A (every counted control killed with margin >= %.2f by >= 1 gate): %d of %d -> %s' % (
            path, CTL_MARGIN, na, len(gcl), 'HOLDS' if na == len(gcl) else 'DOES NOT HOLD'))
        gsl.append('  gate set %s claim B (every control proposed for pre-registration killed with margin >= %.2f, var8 by pTT_pooled alternatively '
                   'z_fix >= %.0f; every kind of control on >= %d fresh seed sets): %d of %d, kinds %s -> %s' % (
            path, CTL_MARGIN, FIX_Z, CTL_SETS, npr, len(pr), 'ok' if kinds_ok else 'NOT ok',
            'HOLDS' if (npr == len(pr) and kinds_ok) else 'DOES NOT HOLD'))
        gsl.append('  gate set %s claim B-resizable (the proposed controls without the fixed var8; kinds on >= %d fresh seed sets): %d of %d, kinds %s -> %s; not with margin: %s' % (
            path, CTL_SETS, nrr, len(prr), 'ok' if kinds_ok else 'NOT ok', 'HOLDS' if (nrr == len(prr) and kinds_ok) else 'DOES NOT HOLD',
            ', '.join('%s %+.3f' % (g['control'], g['prop_margin']) for g in pr if not g['prop_with_margin']) or 'none'))
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
               flavour_table=flav, separation=sep, gate_set_lines=gsl, eF_vs_emax=efl, fault_spread=spr)
    json.dump(res, open(outf, 'w'), indent=1)
    print('207-3 proposed accuracy gates (development, dev seeds; generated by t207_3_summary.py from the cal JSON files)')
    print('minimum stratum size rule: a stratum enters by_stratum if it has >= %d point(s) (the approved gate\'s rule: every stratum); '
          'smallest stratum seen: %d points (path %s, measure %s, group %s, member %s, level %s)' % ((mp,) + sm))
    print('T2 pooling assumed: 1 or 2 main runs per path (each main run alone and every pair of main runs)')
    print('\n'.join(lines))
    print('== separation of every gate (decision 37 + margin rule): every correct variant incl. the device passes in every pooling; pTT gates: every counted control fails; emax gates: the PROPOSED controls fail, the rest is information')
    for k, v in sorted(PROBES.items()):
        print('  probe %s: %s' % (k, v))
    print('\n'.join(sep))
    print('== gate sets (pTT gates; computed from the rows above; rule: control killed with margin >= %.2f by >= 1 gate of the set)%s' % (
        CTL_MARGIN, '; TEST ONLY bound override ' + TEST_BOUND if TEST_BOUND else ''))
    print('\n'.join(gsl))
    print('== spread of the graded-fault runs (S_ens(pTT) pooled per run, and per flavour)')
    print('\n'.join(spr))
    print('== flavour 1 vs flavour 2 (S_ens per flavour; device / correct / var8: main runs pooled; controls: their run)')
    print('\n'.join(flav))
    print('== why eF is not gated: if eF puts the same controls above 0.15 as emax (pooled S_ens, last line per path), it adds no '
          'detection; emax also carries the by_stratum gate')
    print('\n'.join(efl))


if __name__ == '__main__':
    main()
