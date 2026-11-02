#!/usr/bin/env python3
"""Ticket 207 dev5 (development only): calibration table of the accuracy test against the quad reference.

Per run (seed set) and pooled over the runs of one path, for every method m and a measure (default pTT: max|B dGTT Bi|, the error
of the CGR2_2 output as the next precision check receives it; emax: pointwise max-relative error):
  S_cpu(m) = mean over points of log10(e_m / e_cpu)            (reference: the host CGR2_2 of the oracle algorithm, same inputs)
  S_ens(m) = mean over points of log10(e_m / g)                 (g: geometric mean of e over the members of the --ref list at the
                                                                  point; leave-one-out for a member)
with a bootstrap 95 % interval over points, and the share of bootstrap means above each candidate bound.

usage: t207_d5_agg.py OUT.json MEASURE --ref def|log|LIST [--host LABEL] [--allow-v1] [--prog-sha SHA256] ROLE:LABEL=path/qref.txt ...
  ROLE: main (device = 'dev', stand-in control = 'var8') | ctl (production control run: its 'dev' is the defect)
  --ref: the reference ensemble, an explicit list (ticket 207-2): 'def' = cpu,H,n4,n100004,n1,n100001,n2,n1608,n4808,n5,n100005,
         'log' = cpu,H,n4,n100004,n1,n100001,n2,n5,n100005 (H: the host member label), or a comma list of names cpu, the host label,
         n<code> (code != 0).  An unknown name, a duplicate, or a listed name missing from any run is an error (no method is ever treated
         as correct by default).
  --host: the label of the in-run host CGR2_2 member (it depends on the build of the offline program, e.g. nvh NVIDIA, ifh Intel);
         default: the label the runs carry at the host position (t207_d5_qparse.load), which must be the same in every run.
  ROLE is tied to the run's own env (ticket 207-2, review): a main run has exactly ALF_GPU_{DEBUG=1, TAUDIAG_QREF, TAUDIAG_FROM,
         TAUDIAG_QHOST=1, TAUDIAG_QVAR=8}, a ctl run the same plus TAU_NS=0 (default path) or TAUSWAP=8 (LOG path) (QVAR optional); any
         other ALF_GPU_* variable, or a role that does not match the env, is an error.
  The same run twice (same qref.txt sha256, same run.log sha256, or the same role + seeds + parameters) is an error.
  Every run needs a qref_proc.json next to its qref.txt (written by t207_d5_proc.sh, which checks the binary against its .sha256 file),
         checked by t207_d5_qparse.proc_record: rc 0, check_rc 0, the arguments of the spec, valid program and source sha256, and the
         sha256 of qref.txt as it is now; a run without it is an error.  --prog-sha (optional here): every run's program hash must equal it.
This script is a development and report tool: the accuracy gate (t207_d5_gate.py) does not read its output; it imports ref_list, per_run,
pool_values and summ from here and recomputes everything from the run directories (one implementation of the statistics).
Input is read through the strict parser t207_d5_qparse.py (every point, every member, every number checked); any deviation ends the
run with an error and no output file.  Measures pTT / pmax use the points with B (their number must equal the number derived from
the run parameters); emax / eF use all points.  An error number that is zero or not finite is an error (no clamping).
"""
import importlib.util, json, math, os, random, re, statistics, sys
sys.dont_write_bytecode = True      # the loaders below import sibling scripts by path: no __pycache__ in the scripts directory

_spec = importlib.util.spec_from_file_location('t207_d5_qparse', os.path.join(os.path.dirname(os.path.abspath(__file__)), 't207_d5_qparse.py'))
QP = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(QP)
QError = QP.QError


def ref_list(arg, host):
    """the explicit reference ensemble: validated names, in the given order; 'def' / 'log' are the frozen lists with the host member
    labelled host (t207_d5_qparse.ref_names)"""
    QP.check_host(host)
    names = QP.ref_names(arg, host) if arg in ('def', 'log') else arg.split(',')
    seen = set()
    for n in names:
        if not (n in ('cpu', host) or re.match(r'^n[1-9][0-9]*$', n)):
            raise QError('--ref: %r is not a reference name (allowed: cpu, the host label %s, n<code> with code != 0)' % (n, host))
        if n in seen:
            raise QError('--ref: duplicate name %s' % n)
        seen.add(n)
    if 'cpu' not in seen or len(names) < 2:
        raise QError('--ref: the list must contain cpu and at least one more member')
    return names


def parse(fn, allow_v1=False):
    """the validated run of one qref.txt: list of (Point), n points with B"""
    run = QP.load(fn, allow_v1)
    return run


def meas(e, which):
    if which == 'emax':
        return max(e[0:4])
    if which == 'eF':
        return e[4]
    if which == 'eTT':
        return e[2]
    if len(e) != 8 or e[5] < 0:
        raise QError('measure %s needs the propagated numbers (point without B)' % which)
    return e[5] if which == 'pTT' else max(e[5:8])


def q(v, p):
    v = sorted(v)
    k = (len(v) - 1) * p
    i = int(math.floor(k)); j = min(i + 1, len(v) - 1)
    return v[i] + (v[j] - v[i]) * (k - i)


def boot_means(v, n=4000, seed=20261030):
    rng = random.Random(seed)
    L = len(v)
    return [sum(v[rng.randrange(L)] for _ in range(L)) / L for _ in range(n)]


def logs(run, which):
    """per point: dict method -> log10 e, for the points the measure exists for (all, or those with B); returns (list, derived count)"""
    out = []
    for pt in run.points:
        if which in ('pTT', 'pmax') and not pt.hasb:
            continue
        r = {}
        for k, e in pt.m.items():
            if k in QP.DIAG_NAMES:
                continue
            x = meas(e, which)
            if not (math.isfinite(x) and x > 0):
                raise QError('dump %d member %s: error number %r is not positive and finite (log undefined; no clamping)' % (pt.hd['file'], k, x))
            r[k] = math.log10(x)
        out.append(r)
    n_exp = run.n_with_b if which in ('pTT', 'pmax') else len(run.points)
    if len(out) != n_exp or not out:
        raise QError('%d points used for measure %s, %d derived from the run parameters' % (len(out), which, n_exp))
    return out, n_exp


def proc_record(qref_txt, spec, progsha):
    """qref_proc.json of the processing (always required): must say rc 0, check_rc 0 and the arguments of the spec and carry a valid program
    sha256 (equal to --prog-sha if given); returns the recorded program sha256"""
    rundir = os.path.dirname(os.path.abspath(qref_txt))
    pj = os.path.join(rundir, 'qref_proc.json')
    if not os.path.isfile(pj):
        raise QError('%s is missing: every run needs the record of its processing (t207_d5_proc.sh) with the program hash' % pj)
    try:
        pr = QP.proc_record(rundir, spec)              # the one implementation of the record check (qparse; also used by thin and the gate)
    except QError as e:
        raise QError('%s: not a passing processing of the arguments in qref_spec.json (rc, check_rc, args, hashes): %s' % (pj, e))
    if progsha is not None and pr['prog'] != progsha:
        raise QError('%s: program sha256 %s differs from --prog-sha %s' % (pj, pr['prog'], progsha))
    return pr['prog']


def positions(run, which):
    """per point used by measure which (same order as logs()): its position (nst, nfe, br) in the call cycle (fix round 4: strata)"""
    return [(pt.hd['nst'], pt.hd['nfe'], pt.hd['br']) for pt in run.points if not (which in ('pTT', 'pmax') and not pt.hasb)]


def per_run(run, which, refs):
    """per member: n, rc (S_cpu per point), re (S_ens per point) and pos (the (nst, nfe, br) of every point, for per-stratum statistics)"""
    L, n_exp = logs(run, which)
    pos = positions(run, which)
    if len(pos) != len(L):
        raise QError('positions: %d points, %d values' % (len(pos), len(L)))
    methods = list(L[0])
    for p in L:
        if list(p) != methods:
            raise QError('the member list differs between points')
    for r in refs:
        if r not in methods:
            raise QError('reference member %s is missing in the run (members: %s)' % (r, ' '.join(methods)))
    res = {}
    for m in methods:
        members = [k for k in refs if k != m]
        rc = [p[m] - p['cpu'] for p in L]
        re_ = [p[m] - sum(p[k] for k in members) / len(members) for p in L]
        res[m] = dict(n=len(rc), rc=rc, re=re_, pos=list(pos))
    return res, n_exp


STATISTICS = {'S_cpu': 'rc', 'S_ens': 're'}       # statistic name -> per-point list of per_run()


def pool_values(results, member, statistic):
    """the per-point values of a statistic of one member, concatenated over the per_run() results of several runs (the pooled statistic is
    their mean: every point weighs the same)"""
    key = STATISTICS[statistic]
    out = []
    for res in results:
        out += res[member][key]
    return out


def summ(v):
    b = boot_means(v)
    return dict(n=len(v), mean=statistics.fmean(v), lo=q(b, 0.025), hi=q(b, 0.975), sd=statistics.pstdev(v), boot=b)


def main():
    args = sys.argv[1:]
    allow_v1 = False
    progsha = None
    refarg = None
    host = None
    rest = []
    i = 0
    while i < len(args):
        if args[i] == '--allow-v1':
            allow_v1 = True
        elif args[i] == '--prog-sha':
            if i + 1 >= len(args) or not re.match(r'^[0-9a-f]{64}$', args[i + 1]):
                raise QError('--prog-sha needs a sha256 (64 lower-case hex digits)')
            progsha = args[i + 1]; i += 1
        elif args[i] == '--ref':
            if i + 1 >= len(args):
                raise QError('--ref needs a value')
            refarg = args[i + 1]; i += 1
        elif args[i] == '--host':
            if i + 1 >= len(args):
                raise QError('--host needs a value')
            host = QP.check_host(args[i + 1], '--host'); i += 1
        else:
            rest.append(args[i])
        i += 1
    if refarg is None or len(rest) < 3 or rest[1] not in ('pTT', 'pmax', 'emax', 'eF', 'eTT'):
        raise QError('usage: t207_d5_agg.py OUT.json pTT|pmax|emax|eF|eTT --ref def|log|LIST [--allow-v1] ROLE:LABEL=path/qref.txt ...')
    outf, which = rest[0], rest[1]
    runs = []
    seen = set()
    for a in rest[2:]:
        m = re.match(r'^(main|ctl):([A-Za-z0-9_.-]+)=(.+)$', a)
        if not m:
            raise QError('bad run argument %r (ROLE:LABEL=path)' % a)
        if m.group(2) in seen:
            raise QError('duplicate label %s' % m.group(2))
        seen.add(m.group(2))
        runs.append((m.group(1), m.group(2), m.group(3), parse(m.group(3), allow_v1)))
        run = runs[-1][3]
        role_env = QP.env_role(run.spec['run']['env'], run.spec['run']['variant'])
        if role_env != m.group(1):
            raise QError('run %s: role %s given, but its env is that of a %s run' % (m.group(2), m.group(1), role_env))
    labels = sorted(set(r[3].host for r in runs if r[3].host is not None))
    if host is None:
        if len(labels) != 1:
            raise QError('host member label: the runs carry %s; one label expected (or give --host)' % (labels or 'none'))
        host = labels[0]
    elif labels and labels != [host]:
        raise QError('host member label: the runs carry %s, --host says %s' % (labels, host))
    refs = ref_list(refarg, host)
    table = {}
    pooled = {}
    inputs = {}
    dup = {}
    for role, lab, fn, run in runs:
        res, n_exp = per_run(run, which, refs)
        r = run.spec['run']
        prog = proc_record(fn, run.spec, progsha)
        keys = [('qref.txt sha256', QP.sha256_file(fn)), ('run.log sha256', r['records_sha256']['run.log']),
                ('run record (role, seeds, parameters)', (role, r['seeds_sha256'], r['parameters_sha256']))]
        for what, key in keys:
            if (what, key) in dup:
                raise QError('run %s is the same run as %s (same %s): a run may enter once' % (lab, dup[(what, key)], what))
            dup[(what, key)] = lab
        inputs[lab] = dict(role=role, path=fn, qref_sha256=keys[0][1], run_log_sha256=keys[1][1], seeds_sha256=r['seeds_sha256'],
                           parameters_sha256=r['parameters_sha256'], prog_sha256=prog, n_points=len(run.points), n_used=res['cpu']['n'],
                           n_expected=n_exp, nstm=r['nstm'], nfl=r['nfl'], stride=r['stride'], v1=any(p.hd['ver'] == 1 for p in run.points), variant=r['variant'], args=run.spec['args'],
                           env=r['env'], members=sorted(res))
        table[lab] = {}
        for m, d in res.items():
            name = m
            if role == 'ctl' and m == 'dev':
                name = 'ctl_dev'
            elif role == 'ctl':
                continue                            # the other methods of a control run are not part of the calibration
            sc, se = summ(d['rc']), summ(d['re'])
            table[lab][name] = dict(S_cpu=sc, S_ens=se)
            pooled.setdefault(name, []).append({name: d})
    corr_all = set(refs)
    pool = {m: dict(S_cpu=summ(pool_values(rs, m, 'S_cpu')), S_ens=summ(pool_values(rs, m, 'S_ens'))) for m, rs in pooled.items()}
    # print
    labs = [lab for _, lab, _, _ in runs]
    order = ['dev', 'var8', 'ctl_dev', 'cpu', host] + sorted(m for m in pool if re.match(r'^n[1-9][0-9]*$', m)) + ['u1', 'u2']
    for ref in ('S_cpu', 'S_ens'):
        print('== measure %s, statistic %s (mean log10 ratio; [95%% bootstrap interval over points])' % (which, ref))
        print('  %-9s' % 'method' + ''.join('%24s' % l for l in labs) + '%24s' % 'pooled')
        for m in order:
            if m not in pool:
                continue
            row = '  %-9s' % m
            for l in labs:
                s = table[l].get(m, {}).get(ref)
                row += '%24s' % ('%+.3f [%+.3f,%+.3f]' % (s['mean'], s['lo'], s['hi']) if s else '-')
            s = pool[m][ref]
            row += '%24s' % ('%+.3f [%+.3f,%+.3f]' % (s['mean'], s['lo'], s['hi']))
            print(row)
        # separation and false-alarm view
        nulls = [m for m in pool if m in corr_all and m != 'cpu']
        per_set_null_max = {l: max(table[l][m][ref]['mean'] for m in nulls if m in table[l]) for l in labs if any(m in table[l] for m in nulls)}
        print('  max over correct variants per set: ' + ', '.join('%s %+.3f' % (l, v) for l, v in per_set_null_max.items()))
        for cand in (0.10, 0.15, 0.20, 0.25):
            fa = []
            for l in labs:
                for m in nulls:
                    if m in table[l]:
                        b = table[l][m][ref]['boot']
                        fa.append(sum(1 for x in b if x > cand) / len(b))
            dv = [table[l]['dev'][ref]['mean'] for l in labs if 'dev' in table[l]]
            cv = [table[l][c][ref]['mean'] for l in labs for c in ('var8', 'ctl_dev') if c in table[l]]
            print('  bound %.2f: correct variants: mean bootstrap exceedance %.4f, max %.4f, (set, variant) pairs %d, of them with mean above the bound: %d;'
                  ' device per set max %+.3f (k_min %.2f); controls per set min %+.3f' % (
                      cand, statistics.fmean(fa) if fa else float('nan'), max(fa) if fa else float('nan'), len(fa),
                      sum(1 for l in labs for m in nulls if m in table[l] and table[l][m][ref]['mean'] > cand),
                      max(dv) if dv else float('nan'), 10 ** (cand - max(dv)) if dv else float('nan'), min(cv) if cv else float('nan')))
    for l in table:
        for m in table[l]:
            for ref in ('S_cpu', 'S_ens'):
                table[l][m][ref].pop('boot')
    for m in pool:
        for ref in ('S_cpu', 'S_ens'):
            pool[m][ref].pop('boot')
    json.dump(dict(measure=which, refs=refs, host_member=host, allow_v1=allow_v1, prog_sha_required=progsha, inputs=inputs, runs=table, pooled=pool), open(outf, 'w'), indent=1)


if __name__ == '__main__':
    try:
        main()
    except QError as e:
        sys.stderr.write('t207_d5_agg ERROR: %s\n' % e)
        sys.exit(1)
