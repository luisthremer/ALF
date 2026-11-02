#!/usr/bin/env python3
"""Ticket 302c-2 round 2 (development only, dev seeds): separation tables of the candidate 207 accuracy gates on the B580, generated from
the processed runs with the approved 207-2 code (commit 4118a1de: t207_d5_qparse.load = strict reader, t207_d5_agg.per_run = the per-point
values, agg.q = the gate's quantile), computed as t207_d5_gate.py (4118a1de) computes a gate:

  per point p and member m:  r_m(p) = log10 M_m(p) - mean_{k in refs, k != m} log10 M_k(p)       (M = pTT on the points with B; emax, eF on
                             all points; refs = the frozen list of the path with the in-run host member under the label the runs carry)
  statistic   S_ens (mean of r), S_ens_q95 (95th percentile, linear interpolation), S_ens_max (maximum)
  pool        pooled     : the statistic over the points of all runs of the group
              each_run   : per run; decisive value for a member that must pass = largest run value, for a control = smallest run value
              by_stratum : per stratum (nst, nfe, br) over the points of the group's runs; decisive value = largest stratum value

Groups: 'main' = all main runs (members dev, var8 and every reference member, leave-one-out); every fault run = its member dev.
For every candidate (measure, statistic, pool):  device value, correct-variant maximum (max over the reference members of their must-pass
value), var8, every fault control's value, and the separation: a bound b separates (decision 37: every correct variant and the device
pass, the planted defect fails) iff max(device, correct max) <= b < control value; the range that allows is printed, or 'no bound
separates'.  The NVIDIA-proposed gates (207-3, beta = 4) are evaluated at their bounds; this script chooses and moves no bound.
Also: per run the S_ens(pTT) of every correct member (the beta = 10 LOG n5 / n4 question).

usage: t302c2b_gates.py OUT.json TAG def|log PROG_SHA256 ROLE:LABEL=RUNDIR ...      ROLE: main | flt
  main: env exactly ALF_GPU_{DEBUG=1, TAUDIAG_QHOST=1, TAUDIAG_QVAR=8, TAUDIAG_QREF, TAUDIAG_FROM}; flt: the same without QVAR plus
  ALF_GPU_FAULT_TAU=eps:<x> [ALF_GPU_FAULT_TAU_AT=<nst>:<nfe>].  Every run: qref_proc.json passing for the arguments of its spec and
  program sha256 = PROG_SHA256 (agg.proc_record); one host member label; no run twice (qref.txt / run.log / seeds+role+fault size:
  fault runs of different x on one seed set are allowed, they have the same trajectory and differ only by the fault).
"""
import importlib.util, json, math, os, re, statistics, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location('t207_d5_agg', os.path.join(HERE, 't207_d5_agg.py'))
AG = importlib.util.module_from_spec(_s)
_s.loader.exec_module(AG)
QP = AG.QP
QError = QP.QError
MEASURES = ('pTT', 'emax', 'eF')
STATS = ('S_ens', 'S_ens_q95', 'S_ens_max')
POOLS = ('pooled', 'each_run', 'by_stratum')
# NVIDIA proposals of 207-3 (beta = 4, decision 37 rule), with the controls 207-3 used to show each is sharp; report-only candidates
PROPOSED = [('pTT', 'S_ens', 'pooled', 0.15, 'near eps (2e-15 / 1e-14)'), ('emax', 'S_ens', 'pooled', 0.15, 'eps 1e-12'),
            ('emax', 'S_ens', 'by_stratum', 0.30, 'localized eps 2e-12 @ (2,2)'), ('emax', 'S_ens_q95', 'pooled', 0.60, '(207-3: no own control named)')]
REPORT_ONLY = [('pTT', 'S_ens', 'by_stratum'), ('pTT', 'S_ens_q95', 'pooled'), ('pTT', 'S_ens_max', 'pooled'), ('emax', 'S_ens_max', 'pooled')]


def stat(name, v):
    if not v:
        raise QError('statistic of an empty set')
    if name == 'S_ens':
        return statistics.fmean(v)
    if name == 'S_ens_q95':
        return AG.q(v, 0.95)
    return max(v)


def role_of(env):
    base = {'ALF_GPU_DEBUG': '1', 'ALF_GPU_TAUDIAG_QHOST': '1'}
    for k in ('ALF_GPU_TAUDIAG_QREF', 'ALF_GPU_TAUDIAG_FROM'):
        if not re.match(r'^[0-9]+$', env.get(k, 'x')):
            raise QError('%s missing or not a number' % k)
        base[k] = env[k]
    if env == dict(base, ALF_GPU_TAUDIAG_QVAR='8'):
        return 'main', None
    f = env.get('ALF_GPU_FAULT_TAU', '')
    m = re.match(r'^eps:([0-9.eE+-]+)$', f)
    if m and env == dict(base, ALF_GPU_FAULT_TAU=f):
        return 'flt', m.group(1)
    at = env.get('ALF_GPU_FAULT_TAU_AT', '')
    if m and re.match(r'^[0-9]+:[0-9]+$', at) and env == dict(base, ALF_GPU_FAULT_TAU=f, ALF_GPU_FAULT_TAU_AT=at):
        return 'flt', '%s@%s' % (m.group(1), at)
    raise QError('env %s is neither main nor flt' % ' '.join('%s=%s' % kv for kv in sorted(env.items())))


def values(group, member, measure, st, pool, as_control):
    """decisive value and details of a member over the runs of a group (as t207_d5_gate.values / judge)"""
    pts = []
    for i, x in enumerate(group):
        d = x['res'][measure][member]
        pts += [(i, tuple(p), v) for p, v in zip(d['pos'], d['re'])]
    allv = [v for _, _, v in pts]
    per = [stat(st, [v for i, _, v in pts if i == j]) for j in range(len(group))]
    strata = {}
    for _, p, v in pts:
        strata.setdefault(p, []).append(v)
    strata = {p: (len(v), stat(st, v)) for p, v in sorted(strata.items())}
    pooled = stat(st, allv)
    for v in [pooled] + per + [t[1] for t in strata.values()]:
        if not math.isfinite(v):
            raise QError('%s %s %s: not finite' % (member, measure, st))
    if pool == 'pooled':
        dv = pooled
    elif pool == 'each_run':
        dv = min(per) if as_control else max(per)
    else:
        dv = max(t[1] for t in strata.values())
    top = max(strata.items(), key=lambda kv: kv[1][1])
    return dict(value=dv, pooled=pooled, per_run=per, n=len(allv), n_strata=len(strata),
                top_stratum=dict(pos=list(top[0]), n=top[1][0], value=top[1][1]), stratum_n=[min(t[0] for t in strata.values()), max(t[0] for t in strata.values())])


def main():
    a = sys.argv[1:]
    if len(a) < 5 or a[2] not in ('def', 'log') or not re.match(r'^[0-9a-f]{64}$', a[3]):
        sys.exit(__doc__)
    outf, tag, path, progsha = a[:4]
    runs, dup, hosts = [], {}, set()
    for arg in a[4:]:
        m = re.match(r'^(main|flt):([A-Za-z0-9_.@:-]+)=(.+)$', arg)
        if not m:
            raise QError('bad argument %r' % arg)
        qf = os.path.join(m.group(3), 'qref.txt')
        run = QP.load(qf)
        role, x = role_of(run.spec['run']['env'])
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
    main_runs = [x for x in runs if x['role'] == 'main']
    flts = [x for x in runs if x['role'] == 'flt']
    if not main_runs:
        raise QError('no main run')
    bins = sorted(set(x['binary'] for x in runs))
    out = dict(tag=tag, path=path, host_label=host, refs=refs, prog_sha256=progsha, alf_binaries=bins,
               runs=[{k: v for k, v in x.items() if k not in ('run', 'res')} for x in runs], candidates={}, members_pTT_per_run={})
    for w in MEASURES:
        for st in STATS:
            for pool in POOLS:
                key = '%s/%s/%s' % (w, st, pool)
                dev = values(main_runs, 'dev', w, st, pool, False)
                var8 = values(main_runs, 'var8', w, st, pool, True)
                corr = {mm: values(main_runs, mm, w, st, pool, False) for mm in refs}
                cmax_m = max(corr, key=lambda mm: corr[mm]['value'])
                lower = max(dev['value'], corr[cmax_m]['value'])
                ctl = {}
                for f in flts:
                    c = values([f], 'dev', w, st, pool, True)
                    c['separates'] = c['value'] > lower
                    c['range'] = [lower, c['value']] if c['separates'] else None
                    c['margin_vs_correct_max'] = c['value'] - corr[cmax_m]['value']
                    ctl[f['label']] = dict(c, x=f['x'])
                out['candidates'][key] = dict(device=dev, var8=var8, correct=corr, correct_max=dict(member=cmax_m, value=corr[cmax_m]['value']),
                                              lower=lower, controls=ctl)
    for x in main_runs:
        out['members_pTT_per_run'][x['label']] = {mm: statistics.fmean(x['res']['pTT'][mm]['re']) for mm in refs + ['dev', 'var8']}
    with open(outf, 'w') as f:
        json.dump(out, f, indent=1)
    report(out)


def report(o):
    C = o['candidates']
    print('302c-2 round 2 gate tables %s (path %s), generated from %s; host member label %s; refs %s' % (o['tag'], o['path'], ', '.join(
        x['label'] for x in o['runs']), o['host_label'], ','.join(o['refs'])))
    print('  binary128 program %s; ALF binaries %s' % (o['prog_sha256'][:16], ' '.join(b[:16] for b in o['alf_binaries'])))
    for x in o['runs']:
        print('  run %-26s %-4s x %-14s seeds %s points %d (with B %d) keep %d  %s' % (x['label'], x['role'], x['x'], x['seeds'][:16], x['npts'], x['nb'], x['keep'], x['dir']))
    mains = [x['label'] for x in o['runs'] if x['role'] == 'main']
    print('== per main run S_ens(pTT) of every member (correct members leave-one-out; sorted)')
    for l, d in o['members_pTT_per_run'].items():
        print('  %-26s ' % l + ' '.join('%s %+.3f' % kv for kv in sorted(d.items(), key=lambda kv: -kv[1])))
    print('== candidates: value per (measure / statistic / pool); device and correct members judged as must-pass, controls as must-fail')
    print('   correct max = max over the reference members; range = bounds b with max(device, correct max) <= b < control (decision 37)')
    for key, c in C.items():
        dv, cm = c['device'], c['correct_max']
        print('-- %s: device %+.3f (runs %s), correct max %+.3f (%s; per run of it %s), var8 %+.3f' % (
            key, dv['value'], ' '.join('%+.3f' % v for v in dv['per_run']), cm['value'], cm['member'],
            ' '.join('%+.3f' % v for v in c['correct'][cm['member']]['per_run']), c['var8']['value']))
        if key.endswith('by_stratum'):
            print('     strata: device %d strata (%d..%d points), top %s; correct max member top %s' % (
                dv['n_strata'], dv['stratum_n'][0], dv['stratum_n'][1], fs(dv['top_stratum']), fs(c['correct'][cm['member']]['top_stratum'])))
        for l, k in c['controls'].items():
            print('     control %-28s %+.3f  margin to correct max %+.3f  %s%s' % (
                l, k['value'], k['margin_vs_correct_max'], ('separable: %+.3f <= b < %+.3f' % tuple(k['range'])) if k['separates'] else 'no bound separates',
                ('  top stratum %s' % fs(k['top_stratum'])) if key.endswith('by_stratum') else ''))
    print('== NVIDIA-proposed gates (207-3, beta = 4) at their bounds on these runs')
    for w, st, pool, b, what in PROPOSED:
        c = C['%s/%s/%s' % (w, st, pool)]
        dv, cm = c['device']['value'], c['correct_max']['value']
        ok_pass = dv <= b and cm <= b
        print('  %s %s %s <= %.2f [control in 207-3: %s]: device %+.3f %s, correct max %+.3f (%s) %s; bound - max(device, correct) = %+.3f' % (
            w, st, pool, b, what, dv, 'pass' if dv <= b else 'FAIL', cm, c['correct_max']['member'], 'pass' if cm <= b else 'FAIL', b - max(dv, cm)))
        for l, k in c['controls'].items():
            print('      %-28s %+.3f  %s (control - bound %+.3f)' % (l, k['value'], 'killed' if k['value'] > b else 'NOT KILLED', k['value'] - b))
        print('      all must-pass members below the bound: %s' % ('yes' if ok_pass else 'NO'))
    print('== report-only candidates of 207-3: ' + '; '.join('%s/%s/%s device %+.3f correct max %+.3f' % (
        w, st, pool, C['%s/%s/%s' % (w, st, pool)]['device']['value'], C['%s/%s/%s' % (w, st, pool)]['correct_max']['value']) for w, st, pool in REPORT_ONLY))


def fs(t):
    return '(%s):%d:%+.3f' % (','.join(str(v) for v in t['pos']), t['n'], t['value'])


if __name__ == '__main__':
    try:
        main()
    except QError as e:
        sys.stderr.write('t302c2b_gates ERROR: %s\n' % e)
        sys.exit(1)
