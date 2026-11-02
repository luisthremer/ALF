#!/usr/bin/env python3
"""Ticket 302c-2 round 3 (development only, dev seeds): sizing tables of the candidate 207 accuracy gates on the B580, generated from
processed runs.  Reader, per-point values and statistics are round 2's (t302c2b_gates.py: approved 207-2 code of 4118a1de, t207_d5_qparse
strict reader, t207_d5_agg.per_run, r_m(p) = log10 M_m(p) - mean of the other reference members, S_ens / S_ens_q95 / S_ens_max, pools
'pooled' and 'by_stratum' with strata (nst, nfe, br), no minimum stratum size = the approved gate's rule).  Added in round 3:

 (a) groups of main runs: every main run alone ('run:<l>'), every pair of main runs ('pair:<l1>+<l2>'; T2 pools 2 main runs per path and
     beta, 207-3c manager decision) and all main runs ('all').  The decisive value of a must-pass member (device 'dev', every frozen
     reference member = correct variant, leave-one-out) is its largest single-run value (the margin rule uses the single run); the pair
     maximum is reported next to it.  var8 (fixed control, a member of the main runs) is given at its smallest single-run / pair value.
 (b) S_ens_flmax = max(S_ens of the flavour-1 points, S_ens of the flavour-2 points) of the group (207-3c 'pTT_flavour_max'), and per
     flavour S_ens(pTT) of every member in every main run and of the device in every control run.
 (c) per candidate gate (measure / statistic / pool), decision 37's allowed range of the bound b, restricted:
        lower L = max(correct max single run, device max single run) + MARGIN_CORRECT (0.05)  (the device gets the same margin; flagged
                  if it is the binding one)
        control k: upper U_k = value_k - MARGIN_CONTROL (0.10);  the range for a chosen set of controls is L <= b <= min U_k
     printed per control, and per control class with the smallest size of the class that admits the range (the 'chosen' control of the
     class for that gate); 'NO BOUND' where U < L.
 (d) per control class (ladder of fault sizes) and gate: the smallest size killed with margin >= 0.10 at b = L, from which every larger
     size of the ladder is killed with that margin too; several runs of one size: the smallest value counts.

 (c-rule) manager rule (2026-11-02, from the 207-3c/3d review, for 207-6; with the manager's correction of the same day): a must-pass
     member passes with margin iff its worst single-run margin >= 0.03 and z = (b - its mean over the main runs) / its run-to-run sd >= 2.5,
     so L_rule = max over members of max(max_single + 0.03, mean + 2.5 sd); a control is killed with margin iff its margin >= 0.10 on each of
     >= 2 fresh seed sets ('one set' where only one was measured).  Fields L_rule, need, classes_rule, controls_rule.

Control classes (from the run's fault env): 'all-flavour' = eps:x without AT (judged by pTT gates: class 'pTT all-flavour'; by emax gates:
class 'emax'), 'flavour-only' = AT 0:<nfe>, 'localized' = AT <nst>:<nfe> with nst > 0.  Runs of another work dir (round 2) enter with
their label; their binaries are listed in the JSON.

usage: t302c2c_gates.py OUT.json TAG def|log PROG_SHA256 ROLE:LABEL=RUNDIR ...      ROLE: main | flt   (checks as t302c2b_gates.py)
"""
import importlib.util, itertools, json, math, os, re, statistics, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location('t302c2b_gates', os.path.join(HERE, 't302c2b_gates.py'))
G2 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(G2)
AG, QP, QError = G2.AG, G2.QP, G2.QError
MARGIN_CORRECT, MARGIN_CONTROL = 0.05, 0.10
RULE_MIN_MARGIN, RULE_Z = 0.03, 2.5          # manager rule 2026-11-02 (207-3c/3d review, for 207-6); see the docstring
MEASURES = ('pTT', 'emax', 'eF')
CANDS = [(w, st, pool) for w in MEASURES for st in ('S_ens', 'S_ens_q95', 'S_ens_max') for pool in ('pooled', 'by_stratum')] + \
        [(w, 'S_ens_flmax', 'pooled') for w in MEASURES]
FAMILY = {'pTT': ('pTT all-flavour', 'pTT flavour-only', 'localized'), 'emax': ('emax', 'localized'), 'eF': ('emax', 'localized')}


def ckey(c):
    return '%s/%s/%s' % c


def pts_of(group, member, measure):
    out = []
    for x in group:
        d = x['res'][measure][member]
        out += [(tuple(p), v) for p, v in zip(d['pos'], d['re'])]
    return out


def value(group, member, cand):
    """value of a candidate statistic for one member over a group of runs (pooled over the group's points)"""
    w, st, pool = cand
    pts = pts_of(group, member, w)
    if st == 'S_ens_flmax':
        f = [statistics.fmean([v for p, v in pts if p[1] == k]) for k in (1, 2)]
        v = max(f)
    elif pool == 'pooled':
        v = G2.stat(st, [v for _, v in pts])
    else:
        strata = {}
        for p, x in pts:
            strata.setdefault(p, []).append(x)
        v = max(G2.stat(st, s) for s in strata.values())
    if not math.isfinite(v):
        raise QError('%s %s: not finite' % (member, ckey(cand)))
    return v


def flav(group, member):
    pts = pts_of(group, member, 'pTT')
    return [statistics.fmean([v for p, v in pts if p[1] == k]) for k in (1, 2)]


def fclass(x):
    """(ladder, size, position) of a fault run from its x field 'eps' or 'eps@nst:nfe'"""
    if '@' in x:
        e, at = x.split('@')
        nst = int(at.split(':')[0])
        return ('flavour-only' if nst == 0 else 'localized'), float(e), at
    return 'all-flavour', float(x), ''


def cls_for(ladder, measure):
    if ladder == 'all-flavour':
        return 'pTT all-flavour' if measure == 'pTT' else 'emax'
    return 'pTT flavour-only' if ladder == 'flavour-only' else ladder


def main():
    a = sys.argv[1:]
    if len(a) < 5 or a[2] not in ('def', 'log') or not re.match(r'^[0-9a-f]{64}$', a[3]):
        sys.exit(__doc__)
    outf, tag, path, progsha = a[:4]
    runs, dup, hosts = [], {}, set()
    for arg in a[4:]:
        m = re.match(r'^(main|flt):([A-Za-z0-9_.@:/-]+)=(.+)$', arg)
        if not m:
            raise QError('bad argument %r' % arg)
        qf = os.path.join(m.group(3), 'qref.txt')
        run = QP.load(qf)
        role, x = G2.role_of(run.spec['run']['env'])
        if role != m.group(1):
            raise QError('%s: given role %s, env says %s' % (m.group(2), m.group(1), role))
        if run.spec['run']['variant'] != path or run.spec['args']['path'] != path or run.spec['args']['cover'] != 1:
            raise QError('%s: not a cover-1 %s run' % (m.group(2), path))
        prog = AG.proc_record(qf, run.spec, progsha)
        r = run.spec['run']
        for what, key in (('qref.txt', QP.sha256_file(qf)), ('run.log', r['records_sha256']['run.log']), ('seeds+role+fault', (role, x, r['seeds_sha256']))):
            if (what, key) in dup:
                raise QError('%s is the same run as %s (%s)' % (m.group(2), dup[(what, key)], what))
            dup[(what, key)] = m.group(2)
        hosts.add(run.host)
        runs.append(dict(label=m.group(2), dir=m.group(3), role=role, x=x, run=run, prog=prog, seeds=r['seeds_sha256'], keep=run.spec['args']['keep'],
                         binary=r.get('binary_sha256'), npts=len(run.points), nb=run.n_with_b))
    if len(hosts) != 1 or None in hosts:
        raise QError('host member labels %s: one expected' % hosts)
    host = hosts.pop()
    refs = QP.ref_names(path, host)
    for x in runs:
        x['res'] = {w: AG.per_run(x['run'], w, refs)[0] for w in MEASURES}
    mains = [x for x in runs if x['role'] == 'main']
    flts = [x for x in runs if x['role'] == 'flt']
    if len(mains) < 2:
        raise QError('at least two main runs needed (pairs)')
    groups = [('run:' + x['label'], [x]) for x in mains] + [('pair:%s+%s' % (p['label'], q['label']), [p, q]) for p, q in itertools.combinations(mains, 2)] + \
             [('all', mains)]
    out = dict(tag=tag, path=path, host_label=host, refs=refs, prog_sha256=progsha, alf_binaries=sorted(set(x['binary'] for x in runs)),
               margin_correct=MARGIN_CORRECT, margin_control=MARGIN_CONTROL, groups=[g for g, _ in groups],
               runs=[{k: v for k, v in x.items() if k not in ('run', 'res')} for x in runs], candidates={}, flavour_pTT={}, members_pTT={})
    for f in flts:
        f['ladder'], f['size'], f['at'] = fclass(f['x'])
    for cand in CANDS:
        k = ckey(cand)
        vals = {g: {m: value(rs, m, cand) for m in refs + ['dev', 'var8']} for g, rs in groups}
        single = [g for g, _ in groups if g.startswith('run:')]
        pairs = [g for g, _ in groups if g.startswith('pair:')]

        def worst(gs, members, hi=True):
            z = [(vals[g][m], g, m) for g in gs for m in members]
            return list(max(z) if hi else min(z))
        cs, cp, ca = worst(single, refs), worst(pairs, refs), worst(['all'], refs)
        ds, dp = worst(single, ['dev']), worst(pairs, ['dev'])
        v8s, v8p = worst(single, ['var8'], False), worst(pairs, ['var8'], False)
        L = max(cs[0], ds[0]) + MARGIN_CORRECT
        # manager rule (2026-11-02, from the 207-3c/3d review, going into 207-6): a must-pass member passes with margin iff its single-run
        # margin (bound - largest single-run value) >= RULE_MIN_MARGIN and z = (bound - mean of its single-run values over the main runs) / sd
        # >= RULE_Z, sd = run-to-run sd of those values (manager correction 2026-11-02); so the lowest admissible bound for member m is
        # max(max_single(m) + RULE_MIN_MARGIN, mean(m) + RULE_Z * sd(m)); L_rule = max over members (device included)
        need = {}
        for m in refs + ['dev']:
            sv = [vals[g][m] for g in single]
            sd = statistics.stdev(sv)
            mu = statistics.fmean(sv)
            need[m] = dict(max_single=max(sv), mean=mu, sd=sd, lowest_bound=max(max(sv) + RULE_MIN_MARGIN, mu + RULE_Z * sd))
        bind = max(need, key=lambda m: need[m]['lowest_bound'])
        L_rule = need[bind]['lowest_bound']
        ctl = {}
        for f in flts:
            v = value([f], 'dev', cand)
            U = v - MARGIN_CONTROL
            ctl[f['label']] = dict(x=f['x'], ladder=f['ladder'], cls=cls_for(f['ladder'], cand[0]), size=f['size'], at=f['at'], value=v, upper=U,
                                   range=[L, U] if U >= L else None, killed_at_L=v > L, margin_at_L=v - L)
        # classes: the smallest size admitting the range (every larger size of the class too), and the range it gives
        classes = {}
        for cl in sorted(set(c['cls'] for c in ctl.values())):
            sizes = sorted(set(c['size'] for c in ctl.values() if c['cls'] == cl))
            per = {s: min(c['upper'] for c in ctl.values() if c['cls'] == cl and c['size'] == s) for s in sizes}
            ok = [s for s in sizes if all(per[t] >= L for t in sizes if t >= s)]
            sm = ok[0] if ok else None
            classes[cl] = dict(sizes=sizes, upper_by_size={'%g' % s: per[s] for s in sizes}, smallest_killed=sm,
                               range=[L, per[sm]] if sm is not None else None)
        # the same under the manager rule: lower end L_rule; per (class, size, position) the number of seed sets (a kill 'with margin'
        # needs margin >= 0.10 on each of >= 2 fresh seed sets; fewer: 'one set')
        classes_rule = {}
        for cl in sorted(set(c['cls'] for c in ctl.values())):
            sizes = sorted(set(c['size'] for c in ctl.values() if c['cls'] == cl))
            per = {s: min(c['upper'] for c in ctl.values() if c['cls'] == cl and c['size'] == s) for s in sizes}
            nsets = {s: len(set(next(f['seeds'] for f in flts if f['label'] == l) for l, c in ctl.items() if c['cls'] == cl and c['size'] == s))
                     for s in sizes}
            ok = [s for s in sizes if all(per[t] >= L_rule for t in sizes if t >= s)]
            sm = ok[0] if ok else None
            classes_rule[cl] = dict(sizes=sizes, upper_by_size={'%g' % s: per[s] for s in sizes}, n_sets_by_size={'%g' % s: nsets[s] for s in sizes},
                                    smallest_killed=sm, smallest_killed_n_sets=nsets[sm] if sm is not None else None,
                                    meets_two_sets=(nsets[sm] >= 2) if sm is not None else None,
                                    range=[L_rule, per[sm]] if sm is not None else None)
        out['candidates'][k] = dict(
            rule=dict(min_margin=RULE_MIN_MARGIN, z=RULE_Z, ctl_margin=MARGIN_CONTROL, ctl_min_sets=2), L_rule=L_rule, L_rule_member=bind,
            need=need, classes_rule=classes_rule,
            controls_rule={l: dict(margin_at_L_rule=c['value'] - L_rule, killed_with_margin_one_set=c['value'] - L_rule >= MARGIN_CONTROL) for l, c in ctl.items()},
            values=vals, correct_max_single=cs, correct_max_pair=cp, correct_max_all=ca, device_max_single=ds, device_max_pair=dp,
            var8_min_single=v8s, var8_min_pair=v8p, lower=L, lower_binding='device' if ds[0] > cs[0] else 'correct',
            var8_margin_at_L=v8s[0] - L, controls=ctl, classes=classes)
    for g, rs in groups:
        if g.startswith('run:'):
            out['flavour_pTT'][g] = {m: flav(rs, m) for m in refs + ['dev', 'var8']}
            out['members_pTT'][g] = {m: value(rs, m, ('pTT', 'S_ens', 'pooled')) for m in refs + ['dev', 'var8']}
    for f in flts:
        out['flavour_pTT'][f['label']] = {'dev': flav([f], 'dev')}
    with open(outf, 'w') as fo:
        json.dump(out, fo, indent=1)
    print('wrote %s: %s path %s, %d main runs, %d controls, %d candidates' % (outf, tag, path, len(mains), len(flts), len(CANDS)))


if __name__ == '__main__':
    try:
        main()
    except QError as e:
        sys.stderr.write('t302c2c_gates ERROR: %s\n' % e)
        sys.exit(1)
