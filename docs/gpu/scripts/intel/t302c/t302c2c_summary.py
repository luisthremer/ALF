#!/usr/bin/env python3
"""Ticket 302c-2 round 3 (development only): markdown tables generated from tab/gates_<tag>.json (t302c2c_gates.py).  No number is typed.
  usage: t302c2c_summary.py TABDIR                 all tables (a)-(d) per beta and path
         t302c2c_summary.py TABDIR check TAG GATE B   bound check of one gate at bound B on one tag (b4_def, b10_log, ...): every must-pass
                                                    member (device, correct variants) in every single main run and pair passes iff value <= B
                                                    (and is flagged if its margin < margin_correct); every control fails iff value > B
                                                    (flagged if its margin < margin_control); exit 0 iff all pass and all controls fail.
Tables per tag (beta, path):
  (a) candidates: device max single run (pair), correct max single run (member, run) and pair, var8 min single run, L = lower end of the
      allowed range (correct and device single-run max + margin_correct), every control's value with U = value - margin_control, the
      decision-37 range [L, U] or NO BOUND
  (b) per-flavour S_ens(pTT): every main run (device, var8, every correct member: nfe1 / nfe2) and every control run (device), with
      pTT_flavour_max = max(nfe1, nfe2)
  (c) per gate and control class: smallest size of the class killed with margin >= margin_control at b = L (every larger size too) and
      the range it allows, or NO BOUND
  (d) per control class at beta = 10: smallest size killed with margin >= margin_control by any gate (at its L), with the gate(s)"""
import json, os, sys

D = sys.argv[1]
TAGS = ('b4_def', 'b4_log', 'b10_def', 'b10_log')
GATES_A = ['pTT/S_ens/pooled', 'pTT/S_ens_flmax/pooled', 'pTT/S_ens/by_stratum', 'pTT/S_ens_q95/pooled', 'pTT/S_ens_max/pooled',
           'emax/S_ens/pooled', 'emax/S_ens/by_stratum', 'emax/S_ens_q95/pooled', 'emax/S_ens_max/pooled', 'emax/S_ens_flmax/pooled',
           'eF/S_ens/pooled', 'eF/S_ens/by_stratum', 'eF/S_ens_q95/pooled']


def load(tag):
    f = os.path.join(D, 'gates_%s.json' % tag)
    return json.load(open(f)) if os.path.exists(f) else None


def cname(c):
    return '%s %s%s' % (c['ladder'][:4], c['x'].split('@')[0], (' @' + c['at']) if c['at'] else '')


def short(g):
    return g.replace('run:', '').replace('pair:', '')


def check(tag, gate, b):
    o = load(tag)
    c = o['candidates'][gate]
    mc, mk = o['margin_correct'], o['margin_control']
    bad = 0
    print('bound check %s %s <= %g (from %s/gates_%s.json; every control of the tag counts as must-fail, sub-threshold probes included; manager rule for must-pass members: worst-run margin >= %.2f and z >= %.1f; controls: margin >= %.2f per set)' % (tag, gate, b, D, tag, c['rule']['min_margin'], c['rule']['z'], mk))
    for g, vals in c['values'].items():
        if g == 'all':
            continue
        for m in o['refs'] + ['dev']:
            v = vals[m]
            ok = v <= b
            bad += not ok
            sd, mu, mx = c['need'][m]['sd'], c['need'][m]['mean'], c['need'][m]['max_single']
            z = (b - mu) / sd if sd > 0 else float('inf')
            rule = b - mx >= c['rule']['min_margin'] and z >= c['rule']['z']
            print('  must-pass %-8s %-40s %+.3f  %s  margin %+.3f  member: worst-run margin %+.3f, z = (b - mean %+.3f) / sd %.3f = %.1f%s' % (
                m, g, v, 'PASS' if ok else 'FAIL', b - v, b - mx, mu, sd, z,
                '' if not g.startswith('run:') else ('  manager rule: with margin' if rule else '  manager rule: NOT with margin')))
    for l, k in c['controls'].items():
        ok = k['value'] > b
        bad += not ok
        print('  control   %-50s %+.3f  %s  margin %+.3f%s' % (l, k['value'], 'FAILS (killed)' if ok else 'PASSES (NOT killed)', k['value'] - b,
                                                              ('  (margin >= %.2f, one set)' % mk) if k['value'] - b >= mk else
                                                              ('  (margin < %.2f)' % mk) if ok else ''))
    v8 = c['var8_min_single']
    print('  fixed control var8 smallest single run %+.3f (%s): %s' % (v8[0], v8[1], 'killed' if v8[0] > b else 'NOT killed'))
    print('  -> %s' % ('every must-pass member passes and every control fails' if bad == 0 else '%d wrong outcome(s)' % bad))
    return bad


if len(sys.argv) > 2 and sys.argv[2] == 'check':
    sys.exit(1 if check(sys.argv[3], sys.argv[4], float(sys.argv[5])) else 0)

best, best_r = {}, {}
for tag in TAGS:
    o = load(tag)
    if o is None:
        continue
    C = o['candidates']
    mains = [r for r in o['runs'] if r['role'] == 'main']
    flts = [r for r in o['runs'] if r['role'] == 'flt']
    print('\n## %s\n' % tag)
    print('runs: main %s; controls %s.  ALF binaries %s; binary128 program %s; refs %s; margins: correct %.2f, control %.2f\n' % (
        ', '.join('%s (%s)' % (r['label'], 't302c-2b' if 't302c-2b' in r['dir'] else 't302c-2c') for r in mains),
        ', '.join(r['label'] for r in flts), ' '.join(b[:16] for b in o['alf_binaries']), o['prog_sha256'][:16], ','.join(o['refs']),
        o['margin_correct'], o['margin_control']))
    ctl_l = [r['label'] for r in flts]
    print('### (a) candidates: must-pass maxima (single run / pair), var8, controls (value, allowed range [L, value - %.2f])\n' % o['margin_control'])
    print('| candidate | device max run (pair) | correct max run (member, run) | correct max pair | var8 min run | L | ' +
          ' | '.join(cname(C[GATES_A[0]]['controls'][l]) for l in ctl_l) + ' |')
    print('|' + '---|' * (6 + len(ctl_l)))
    for g in GATES_A:
        c = C[g]
        cells = []
        for l in ctl_l:
            k = c['controls'][l]
            cells.append('%+.3f (%s)' % (k['value'], ('%.3f..%.3f' % tuple(k['range'])) if k['range'] else 'NO BOUND'))
        print('| %s | %+.3f (%+.3f) | %+.3f (%s, %s) | %+.3f | %+.3f | %.3f%s | %s |' % (
            g, c['device_max_single'][0], c['device_max_pair'][0], c['correct_max_single'][0], c['correct_max_single'][2],
            short(c['correct_max_single'][1]), c['correct_max_pair'][0], c['var8_min_single'][0], c['lower'],
            ' (device)' if c['lower_binding'] == 'device' else '', ' | '.join(cells)))
    print('\n### (b) per-flavour S_ens(pTT) (nfe1 / nfe2; pTT_flavour_max = the larger)\n')
    print('| run | device | var8 | correct members (nfe1/nfe2), sorted by flavour max |')
    print('|---|---|---|---|')
    for g, d in o['flavour_pTT'].items():
        if 'var8' in d:
            cm = sorted(((m, v) for m, v in d.items() if m in o['refs']), key=lambda kv: -max(kv[1]))
            print('| %s | %+.3f/%+.3f | %+.3f/%+.3f | %s |' % (short(g), d['dev'][0], d['dev'][1], d['var8'][0], d['var8'][1],
                                                          ', '.join('%s %+.3f/%+.3f' % (m, a, b) for m, (a, b) in cm)))
        else:
            print('| %s (control) | %+.3f/%+.3f (max %+.3f) | | |' % (g, d['dev'][0], d['dev'][1], max(d['dev'])))
    print('\n### (c) allowed range per gate and control class (L = must-pass single-run max + %.2f; class: smallest size killed with '
          'margin >= %.2f at b = L, every larger size too)\n' % (o['margin_correct'], o['margin_control']))
    print('| gate | L | var8 (fixed) margin at L | class: smallest killed size -> range [L, U] | all classes of the gate with a killed size: range |')
    print('|---|---|---|---|---|')
    for g in GATES_A:
        c = C[g]
        fam = {'pTT': ('pTT all-flavour', 'pTT flavour-only', 'localized'), 'emax': ('emax', 'localized'), 'eF': ('emax', 'localized')}[g.split('/')[0]]
        cells, ups_k = [], []
        for cl in fam:
            k = c['classes'].get(cl)
            if k is None:
                continue
            if k['smallest_killed'] is not None:
                ups_k.append((k['range'][1], cl))
            ups = ', '.join('%s: U %.3f' % (s, u) for s, u in k['upper_by_size'].items())
            cells.append('**%s** %s (%s)' % (cl, ('%g -> %.3f..%.3f' % (k['smallest_killed'], k['range'][0], k['range'][1])) if k['smallest_killed'] is not None
                                            else 'NO BOUND (no size of the class killed with margin)', ups))
            best.setdefault((tag, cl), []).append((k['smallest_killed'] if k['smallest_killed'] is not None else float('inf'), g))
        comb = ('%.3f..%.3f (%s)' % (c['lower'], min(ups_k)[0], ', '.join(cl for _, cl in ups_k))) if ups_k else 'NO BOUND'
        print('| %s | %.3f | %+.3f | %s | %s |' % (g, c['lower'], c['var8_margin_at_L'], '; '.join(cells), comb))
    r0 = C[GATES_A[0]]['rule']
    print('\n### (c-rule) manager rule (2026-11-02): must-pass member passes with margin iff its worst single-run margin >= %.2f and z = (b - its mean over '
          'the main runs) / its run-to-run sd >= %.1f, so L_rule = max over members of max(max_single + %.2f, mean + %.1f sd); control killed with margin iff margin >= %.2f '
          'on each of >= %d fresh seed sets ("one set" = only one set measured)\n' % (r0['min_margin'], r0['z'], r0['min_margin'], r0['z'],
                                                                                    r0['ctl_margin'], r0['ctl_min_sets']))
    print('| gate | L_rule (binding member: max single, mean, sd) | top-3 members by lowest bound | class: smallest killed size (sets) -> range | combined |')
    print('|---|---|---|---|---|')
    for g in GATES_A:
        c = C[g]
        n = c['need'][c['L_rule_member']]
        top = sorted(c['need'].items(), key=lambda kv: -kv[1]['lowest_bound'])[:3]
        fam = {'pTT': ('pTT all-flavour', 'pTT flavour-only', 'localized'), 'emax': ('emax', 'localized'), 'eF': ('emax', 'localized')}[g.split('/')[0]]
        cells, ups_k = [], []
        for cl in fam:
            k = c['classes_rule'].get(cl)
            if k is None:
                continue
            if k['smallest_killed'] is not None:
                ups_k.append((k['range'][1], cl))
                best_r.setdefault((tag, cl), []).append((k['smallest_killed'], g, k['smallest_killed_n_sets']))
            else:
                best_r.setdefault((tag, cl), []).append((float('inf'), g, 0))
            ups = ', '.join('%s: U %.3f (%d set%s)' % (s, u, k['n_sets_by_size'][s], '' if k['n_sets_by_size'][s] == 1 else 's')
                            for s, u in k['upper_by_size'].items())
            cells.append('**%s** %s (%s)' % (cl, ('%g (%s) -> %.3f..%.3f' % (
                k['smallest_killed'], 'one set' if k['smallest_killed_n_sets'] < 2 else '%d sets' % k['smallest_killed_n_sets'],
                k['range'][0], k['range'][1])) if k['smallest_killed'] is not None else 'NO BOUND', ups))
        comb = ('%.3f..%.3f (%s)' % (c['L_rule'], min(ups_k)[0], ', '.join(cl for _, cl in ups_k))) if ups_k else 'NO BOUND'
        print('| %s | %.3f (%s: max %+.3f, mean %+.3f, sd %.3f) | %s | %s | %s |' % (
            g, c['L_rule'], c['L_rule_member'], n['max_single'], n['mean'], n['sd'],
            ', '.join('%s %.3f (max %+.3f mean %+.3f sd %.3f)' % (m, d['lowest_bound'], d['max_single'], d['mean'], d['sd']) for m, d in top), '; '.join(cells), comb))
print('\n## (d) smallest size per control class killed with margin >= 0.10 by a gate (at its L), beta = 10 (beta = 4 for reference)\n')
print('| tag | class | smallest size | gate(s) |')
print('|---|---|---|---|')
for (tag, cl), lst in sorted(best.items(), key=lambda kv: (kv[0][0][:3] != 'b10', kv[0])):
    s = min(x[0] for x in lst)
    gs = [g for x, g in lst if x == s]
    print('| %s | %s | %s | %s |' % (tag, cl, ('%g' % s) if s != float('inf') else 'NONE (no gate kills any size of the class with margin)',
                                     ', '.join(gs) if s != float('inf') else ''))
print('\n## (d-rule) the same under the manager rule (L_rule; control margin >= 0.10; "one set" where only one fresh seed set was measured)\n')
print('| tag | class | smallest size | sets | gate(s) |')
print('|---|---|---|---|---|')
for (tag, cl), lst in sorted(best_r.items(), key=lambda kv: (kv[0][0][:3] != 'b10', kv[0])):
    s = min(x[0] for x in lst)
    gs = [g for x, g, _ in lst if x == s]
    ns = min(n for x, _, n in lst if x == s)
    print('| %s | %s | %s | %s | %s |' % (tag, cl, ('%g' % s) if s != float('inf') else 'NONE (no gate kills any size of the class with margin)',
                                          ('one set' if ns < 2 else '%d sets' % ns) if s != float('inf') else '', ', '.join(gs) if s != float('inf') else ''))
