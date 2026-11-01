#!/usr/bin/env python3
"""Ticket 207-2: accuracy gate against the quad reference, recomputed from the evidence (third design, after review 207-2c).

usage: t207_d5_gate.py --spec SPEC.json --report OUT.json RUNDIR [RUNDIR ...]

The gate decides only on numbers it computes itself from the run directories; no aggregate is ever an input.  For every RUNDIR (the
directory holding the run's records, qref.txt and the processing records of t207_d5_proc.sh):
  1. qref_spec.json re-derived from the run's own records (parameters, info, run.log, stamp.txt, provenance.json); the run's path
     (from info) is the spec's path;
  2. qref_proc.json (t207_d5_qparse.proc_record): rc 0, check_rc 0, the arguments of qref_spec.json, the program sha256 = spec prog_sha,
     the program source sha256 = spec src_sha, and qref_sha256 = the sha256 of qref.txt as it is now;
  3. the dump evidence (t207_d5_qparse.verify_dumps): the manifest lists exactly the covered dumps, the kept dumps match it;
  4. its processing arguments qs, nul, cover are the spec's; its role is the row of the spec's env table (env_roles) that equals its
     ALF_GPU_* env exactly (t207_d5_qparse.role_of_env);
  5. qref.txt is loaded by the strict parser (t207_d5_qparse.load, version 1 dumps never accepted, every number checked, eF within its
     by-construction bounds) with the spec's member list of the role (exactly those members, in that order);
  6. per gate measure, the per-point log ratios S_ens (and S_cpu, report only) by t207_d5_agg.per_run (the same function the aggregation
     uses); the number of points used must equal the spec's expect_points for that measure (pTT: points with B, emax / eF: all points).
Across the runs: the same parameters; no run twice (directory, qref.txt sha256, run.log sha256, seeds, a dump sha256 or a QREF data line
shared by two runs; control seed sets are distinct from main seed sets, decision 33); the number of runs per role equals the spec's n_runs.

Decision, per gate {name, measure, statistic S_ens, pool, bound}:
  pool 'pooled'   : the mean of the per-point values over all runs of the role (every point weighs the same, t207_d5_agg.pool_values);
  pool 'each_run' : every run on its own (device: the largest run mean must be <= bound; control: the smallest run mean must be > bound).
  device (member 'dev' of the main runs): FAIL if its value > bound.
  control (every control of the spec, always required): 'variant' = a member of the main runs (e.g. var8), 'run' = member 'dev' of the
  runs whose env is the control's role.  It is judged on the statistic, pool and bound of the gate it names in 'trips' and must exceed
  the bound (killed); for every gate in its optional 'must_not_trip' list it must stay <= that gate's bound.  Otherwise FAIL.
Verdicts and exit codes (printed as the last line 'GATE <VERDICT>...'; the report carries the same verdict and exit code):
  PASS 0        every gate passes, every control behaves as declared, every gate is tripped by at least one control, spec not test_only
  FAIL 1        a gate's own criterion caught something: device above a bound, a control not killed, or tripping a gate it must not trip
  INPUT_ERROR 2 the gate could not be evaluated: a bad or incomplete spec, a missing or inconsistent file, a count differing from the spec,
                any refusal above, or any unexpected exception (never a FAIL, never a detection)
  DEV-ONLY 3    the spec declares no control, or a gate no control trips: the numbers are computed and reported, never a PASS
  TEST-ONLY 4   all criteria met, but the spec is marked test_only (e.g. the kept dev5 subset, cover != 1): not an acceptance verdict
The report (--report, must not exist) holds the spec, the sha256 of the spec and of this script, t207_d5_agg.py, t207_d5_qparse.py and
t207_d5_qref.F90, per run the identities and per-member statistics, every decisive number, the verdict and the exit code.  It is written
in every case where the report path is usable (also for INPUT_ERROR) and is never read by the gate.

Spec file (JSON object, no duplicate keys, no NaN / Infinity; every key below required unless marked optional, no other key):
  spec_version 1; path 'def' | 'log'; host_member: the label the offline program prints for the in-run host CGR2_2 of QHOST=1
  (t207_d5_qref.F90 prints 'nvh' whatever the run's compiler; on the B580 it is the ifx/MKL host); members: {role: [the exact ordered
  statistical members of that role's runs]} for main and the roles of the run controls (the parser derives the list from the run's
  records with the host member labelled host_member and refuses any difference; diagnostic members qself, qbr, din, cpuh, cpuhx follow
  from qs and hostin); refs: the frozen reference ensemble, at least two distinct members of every role's list, not dev, not n0 (the
  oracle cpu by construction), not a variant control member;
  qs: one integer > 0 (one qs per gate run); nul: the null codes of the processing (list of distinct integers >= 0, with 0 and the codes
  of refs); cover: 1 (another value only with test_only); test_only (optional): a non-empty reason; expect_points: {measure: integer >= 1}
  for exactly the measures of the gates; prog_sha, src_sha: sha256 of the offline program and of t207_d5_qref.F90 (src_sha must also be
  the sha256 of the t207_d5_qref.F90 next to this script); env_roles: {'main': {VAR: value, ...}, ROLE: {'extra_env': {...},
  'omit_env': [...]}, ...} (exact envs, see t207_d5_qparse.env_table); n_runs: {'main': N, ROLE: N, ...} for main and exactly the roles
  of the run controls; gates: [{name, measure (pTT | emax | eF), statistic ('S_ens'), pool ('pooled' | 'each_run'), bound}, ...];
  controls: [{name, kind 'variant', member (a member of main, not dev, not a reference), trips, must_not_trip (optional)} |
             {name, kind 'run', role, trips, must_not_trip (optional)}, ...].
  allow_v1 is refused in any form.
"""
import hashlib, importlib.util, json, math, os, re, statistics, sys
sys.dont_write_bytecode = True      # the loader below imports a sibling script by path: no __pycache__ in the scripts directory
HERE = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location('t207_d5_agg', os.path.join(HERE, 't207_d5_agg.py'))
AGG = importlib.util.module_from_spec(_s)
_s.loader.exec_module(AGG)
QP = AGG.QP                                  # one instance of the parser module, the one the aggregation uses
QError = QP.QError

EXIT = {'PASS': 0, 'FAIL': 1, 'INPUT_ERROR': 2, 'DEV-ONLY': 3, 'TEST-ONLY': 4}
MEASURES = ('pTT', 'emax', 'eF')             # measures of qref.txt the gate accepts (t207_d5_agg.meas); anything else is refused
STATISTICS = ('S_ens',)                      # decision statistics; S_cpu is reported only
POOLS = ('pooled', 'each_run')
SPEC_REQUIRED = {'spec_version', 'path', 'host_member', 'members', 'refs', 'qs', 'nul', 'cover', 'expect_points', 'prog_sha', 'src_sha',
                 'env_roles', 'n_runs', 'gates', 'controls'}
SPEC_OPTIONAL = {'test_only'}
GATE_KEYS = {'name', 'measure', 'statistic', 'pool', 'bound'}
NAME = re.compile(r'^[A-Za-z0-9_.-]+$')


class InputError(Exception):
    pass


def isint(x):
    return isinstance(x, int) and not isinstance(x, bool)


def isnum(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def sha(p):
    return QP.sha256_file(p)


# ------------------------------------------------------------------------------------------------------------------- spec

def load_spec(path):
    try:
        s = QP.read_json(path, 'gate spec')
    except QError as e:
        raise InputError(str(e))
    if not isinstance(s, dict):
        raise InputError('gate spec: not a JSON object')
    if 'allow_v1' in s:
        raise InputError('gate spec: allow_v1 is refused (version 1 dumps have no B, Bi; the gate never accepts them), value %r' % (s['allow_v1'],))
    unknown = sorted(set(s) - SPEC_REQUIRED - SPEC_OPTIONAL)
    missing = sorted(SPEC_REQUIRED - set(s))
    if unknown or missing:
        raise InputError('gate spec: unknown keys %s, missing keys %s' % (unknown, missing))
    if not (isint(s['spec_version']) and s['spec_version'] == 1):
        raise InputError('gate spec: spec_version %r (1 expected)' % (s['spec_version'],))
    if s['path'] not in ('def', 'log'):
        raise InputError('gate spec: path %r (def or log)' % (s['path'],))
    if not (isinstance(s['host_member'], str) and NAME.match(s['host_member'])):
        raise InputError('gate spec: host_member %r: the label the offline program prints for the in-run host CGR2_2 (QHOST) expected'
                         % (s['host_member'],))
    refs = s['refs']
    if not (isinstance(refs, list) and len(refs) >= 2 and all(isinstance(r, str) for r in refs)):
        raise InputError('gate spec: refs %r is not a list of at least two names' % (refs,))
    if len(set(refs)) != len(refs):
        raise InputError('gate spec: refs %s: duplicate name' % refs)
    for r in refs:
        if r == 'dev' or r == 'n0' or r in QP.DIAG_NAMES:
            raise InputError('gate spec: refs: %r is not a reference member (dev is the device under test; n0 is the oracle cpu by '
                             'construction, checked by the parser, and would count twice; %s are diagnostics)' % (r, ', '.join(QP.DIAG_NAMES)))
    if not (isint(s['qs']) and s['qs'] > 0):
        raise InputError('gate spec: qs %r: one integer > 0 expected (one qs per gate run)' % (s['qs'],))
    nul = s['nul']
    if not (isinstance(nul, list) and all(isint(c) and c >= 0 for c in nul) and len(set(nul)) == len(nul) and 0 in nul):
        raise InputError('gate spec: nul %r: a list of distinct integers >= 0 containing 0 expected' % (nul,))
    if not (isint(s['cover']) and s['cover'] >= 1):
        raise InputError('gate spec: cover %r: an integer >= 1 expected' % (s['cover'],))
    if 'test_only' in s and not (isinstance(s['test_only'], str) and s['test_only'].strip()):
        raise InputError('gate spec: test_only must be a non-empty reason')
    if s['cover'] != 1 and 'test_only' not in s:
        raise InputError('gate spec: cover %d: qref.txt must cover every dump of the run (cover 1); another cover only in a spec marked '
                         'test_only' % s['cover'])
    for k in ('prog_sha', 'src_sha'):
        if not (isinstance(s[k], str) and QP.HEX64.match(s[k])):
            raise InputError('gate spec: %s %r is not a sha256 (64 lower-case hex digits)' % (k, s[k]))
    here_src = sha(os.path.join(HERE, 't207_d5_qref.F90'))
    if s['src_sha'] != here_src:
        raise InputError('gate spec: src_sha %s differs from the sha256 %s of the t207_d5_qref.F90 next to the gate' % (s['src_sha'], here_src))
    try:
        table = QP.env_table(s['env_roles'])
    except QError as e:
        raise InputError('gate spec: %s' % e)
    # gates
    gates = s['gates']
    if not (isinstance(gates, list) and gates):
        raise InputError('gate spec: gates must be a non-empty list')
    gnames = []
    for i, g in enumerate(gates):
        if not isinstance(g, dict) or set(g) != GATE_KEYS:
            raise InputError('gate spec: gates[%d]: a dictionary with exactly the keys %s expected, got %r' % (i, sorted(GATE_KEYS), g))
        if not (isinstance(g['name'], str) and NAME.match(g['name'])) or g['name'] in gnames:
            raise InputError('gate spec: gates[%d]: name %r is not a unique name' % (i, g['name']))
        gnames.append(g['name'])
        if g['measure'] not in MEASURES:
            raise InputError('gate spec: gate %s: unknown measure %r (the gate supports %s)' % (g['name'], g['measure'], ', '.join(MEASURES)))
        if g['statistic'] not in STATISTICS:
            raise InputError('gate spec: gate %s: statistic %r is not a decision statistic (%s)' % (g['name'], g['statistic'], ', '.join(STATISTICS)))
        if g['pool'] not in POOLS:
            raise InputError('gate spec: gate %s: pool %r (%s)' % (g['name'], g['pool'], ' | '.join(POOLS)))
        if not isnum(g['bound']):
            raise InputError('gate spec: gate %s: bound %r is not a finite number' % (g['name'], g['bound']))
    ep = s['expect_points']
    want_m = sorted(set(g['measure'] for g in gates))
    if not isinstance(ep, dict) or sorted(ep) != want_m:
        raise InputError('gate spec: expect_points %r must have exactly the measures of the gates %s' % (ep, want_m))
    for m, n in ep.items():
        if not (isint(n) and n >= 1):
            raise InputError('gate spec: expect_points %s: %r is not an integer >= 1' % (m, n))
    # member lists: the exact ordered statistical members of the runs of each role (names come from the spec, not from the scripts)
    mem = s['members']
    if not isinstance(mem, dict) or 'main' not in mem:
        raise InputError('gate spec: members must be a dictionary {role: [member, ...]} with main')
    for role, lst in mem.items():
        if not (isinstance(lst, list) and lst and all(isinstance(x, str) and NAME.match(x) for x in lst) and len(set(lst)) == len(lst)):
            raise InputError('gate spec: members.%s %r: a list of distinct names expected' % (role, lst))
        if 'dev' not in lst or any(r not in lst for r in refs) or any(x in QP.DIAG_NAMES for x in lst):
            raise InputError('gate spec: members.%s must contain dev and every reference member %s and no diagnostic member %s'
                             % (role, refs, QP.DIAG_NAMES))
        if s['host_member'] not in lst:
            raise InputError('gate spec: members.%s does not contain host_member %s' % (role, s['host_member']))
    # controls
    ctl = s['controls']
    if not isinstance(ctl, list):
        raise InputError('gate spec: controls must be a list')
    cnames = []
    run_roles = set()
    try:                                         # a malformed QVAR of main is refused here (the parser would refuse every run)
        QP.parse_int_list(table['main'].get('ALF_GPU_TAUDIAG_QVAR', ''), 'env_roles.main ALF_GPU_TAUDIAG_QVAR')
    except QError as e:
        raise InputError('gate spec: %s' % e)
    for i, c in enumerate(ctl):
        if not isinstance(c, dict):
            raise InputError('gate spec: controls[%d] is not a dictionary' % i)
        kind = c.get('kind')
        keys = {'name', 'kind', 'trips'} | ({'member'} if kind == 'variant' else {'role'} if kind == 'run' else set())
        if kind not in ('variant', 'run') or not (keys <= set(c) <= keys | {'must_not_trip'}):
            raise InputError('gate spec: controls[%d]: kind variant (keys name, kind, member, trips) or run (keys name, kind, role, trips), '
                             'optional must_not_trip; got %r' % (i, c))
        if not (isinstance(c['name'], str) and NAME.match(c['name'])) or c['name'] in cnames:
            raise InputError('gate spec: controls[%d]: name %r is not a unique name' % (i, c['name']))
        cnames.append(c['name'])
        if c['trips'] not in gnames:
            raise InputError('gate spec: control %s: trips %r is not a gate of the spec %s' % (c['name'], c['trips'], gnames))
        mnt = c.get('must_not_trip', [])
        if not (isinstance(mnt, list) and all(x in gnames for x in mnt) and len(set(mnt)) == len(mnt) and c['trips'] not in mnt):
            raise InputError('gate spec: control %s: must_not_trip %r must be a list of distinct gates other than %s' % (c['name'], mnt, c['trips']))
        if kind == 'variant':
            if c['member'] == 'dev' or c['member'] in refs or c['member'] not in mem['main']:
                raise InputError('gate spec: control %s: member %r must be a member of main in members, not dev and not a reference member'
                                 % (c['name'], c['member']))
        else:
            if c['role'] not in table or c['role'] == 'main':
                raise InputError('gate spec: control %s: role %r is not a control role of env_roles %s' % (c['name'], c['role'], sorted(table)))
            run_roles.add(c['role'])
    unused = sorted(set(table) - {'main'} - run_roles)
    if unused:
        raise InputError('gate spec: env_roles %s are used by no control' % unused)
    if set(mem) != {'main'} | run_roles:
        raise InputError('gate spec: members %s must have exactly the roles %s' % (sorted(mem), sorted({'main'} | run_roles)))
    nr = s['n_runs']
    if not isinstance(nr, dict) or set(nr) != {'main'} | run_roles:
        raise InputError('gate spec: n_runs %r must have exactly the roles %s' % (nr, sorted({'main'} | run_roles)))
    for r, n in nr.items():
        if not (isint(n) and n >= 1):
            raise InputError('gate spec: n_runs %s: %r is not an integer >= 1' % (r, n))
    return s, table


# ------------------------------------------------------------------------------------------------------------------- runs

def load_run(rundir, s, table):
    d = os.path.abspath(rundir)
    what = 'run %s' % d
    if not os.path.isdir(d):
        raise InputError('%s: not a directory' % what)
    if not os.path.isfile(os.path.join(d, 'qref.txt')):
        raise InputError('%s: qref.txt missing' % what)
    try:
        rspec = QP.load_spec(d)                      # re-derived from the run's records and compared with the stored one
        if rspec['run']['variant'] != s['path'] or rspec['args']['path'] != s['path']:
            raise InputError('%s: a %s-path run (info), processed as %s; the spec is for the %s path' % (
                what, rspec['run']['variant'], rspec['args']['path'], s['path']))
        pr = QP.proc_record(d, rspec, prog=s['prog_sha'], src=s['src_sha'])
        manifest = QP.verify_dumps(d, rspec)
        for k in ('qs', 'nul', 'cover'):
            if rspec['args'][k] != s[k]:
                raise InputError('%s: processed with %s %r, the spec says %r' % (what, k, rspec['args'][k], s[k]))
        role = QP.role_of_env(rspec['run']['env'], table)
        run = QP.load(os.path.join(d, 'qref.txt'), allow_v1=False, members=s['members'][role], host=s['host_member'])
        r = run.spec['run']
        res = {}
        for m in sorted(set(g['measure'] for g in s['gates'])):
            res[m], n_exp = AGG.per_run(run, m, s['refs'])
            n = res[m]['cpu']['n']
            if not (isint(n) and n == n_exp == s['expect_points'][m]):
                raise InputError('%s: %r points used for measure %s (derived from the run %r), the spec expects %d' % (what, n, m, n_exp, s['expect_points'][m]))
    except QError as e:
        raise InputError('%s: %s' % (what, e))
    lines = [ln for ln in open(os.path.join(d, 'qref.txt')).read().split('\n') if ln.startswith('QREF file ')]
    return dict(dir=d, real=os.path.realpath(d), role=role, run=run, res=res, proc=pr, manifest=manifest, lines=lines,
                qref_sha256=pr['qref_sha256'], run_log_sha256=r['records_sha256']['run.log'], seeds_sha256=r['seeds_sha256'],
                parameters_sha256=r['parameters_sha256'], manifest_sha256=sha(os.path.join(d, 'qref_manifest.txt')))


def cross_checks(runs, s):
    if len(set(x['parameters_sha256'] for x in runs)) != 1:
        raise InputError('the runs have different parameters (parameters sha256 %s)' % sorted(set(x['parameters_sha256'][:16] for x in runs)))
    seen = {}
    for x in runs:
        keys = [('directory', x['real']), ('qref.txt sha256', x['qref_sha256']), ('run.log sha256', x['run_log_sha256']),
                ('seeds (control seed sets must differ from main seed sets, decision 33)', x['seeds_sha256'])]
        keys += [('dump sha256', h) for h in x['manifest'].values()]
        keys += [('QREF data line', ln.split(' | ', 1)[1]) for ln in x['lines']]
        own = set()
        for what, k in keys:
            if (what, k) in own:
                continue
            own.add((what, k))
            if (what, k) in seen:
                raise InputError('run %s is the same run as %s (same %s): a run may enter once' % (x['dir'], seen[(what, k)], what))
        for wk in own:
            seen[wk] = x['dir']
    cnt = {}
    for x in runs:
        cnt[x['role']] = cnt.get(x['role'], 0) + 1
    for c in s['controls']:
        if c['kind'] == 'run' and c['role'] not in cnt:
            raise InputError('control %s is declared in the spec, but no run has its role %s (controls are always required)' % (c['name'], c['role']))
    if 'main' not in cnt:
        raise InputError('no main run (the device under test)')
    if cnt != s['n_runs']:
        raise InputError('runs per role %s, the spec expects %s' % (dict(sorted(cnt.items())), dict(sorted(s['n_runs'].items()))))


# ------------------------------------------------------------------------------------------------------------------- decision

def values(group, member, g):
    """(pooled mean over the points of all runs in group, [per-run means]) of member under gate g's measure and statistic"""
    for x in group:
        if member not in x['res'][g['measure']]:
            raise InputError('run %s (%s): member %s missing (members: %s)' % (x['dir'], x['role'], member, ' '.join(sorted(x['res'][g['measure']]))))
    allv = AGG.pool_values([x['res'][g['measure']] for x in group], member, g['statistic'])
    per = [statistics.fmean(AGG.pool_values([x['res'][g['measure']]], member, g['statistic'])) for x in group]
    for v in [statistics.fmean(allv)] + per:
        if not math.isfinite(v):
            raise InputError('%s %s: statistic not finite' % (member, g['name']))
    return statistics.fmean(allv), per, allv


def judge(v, per, g, as_control):
    """the decisive value: pooled mean, or per run the largest (device / must-not-trip) or smallest (control must trip) run mean"""
    if g['pool'] == 'pooled':
        return v
    return min(per) if as_control else max(per)


def decide(runs, s):
    main = [x for x in runs if x['role'] == 'main']
    gates = {g['name']: g for g in s['gates']}
    out, reasons, rep_g, rep_c = [], [], [], []
    verdict = 'PASS'
    for g in s['gates']:
        v, per, allv = values(main, 'dev', g)
        dv = judge(v, per, g, False)
        b = AGG.summ(allv); b.pop('boot')
        ok = dv <= g['bound']
        out.append('gate %s: device %s(%s) %s = %+.4f (pooled over %d points %+.4f [%+.4f, %+.4f]; runs %s), bound %+.4f: %s' % (
            g['name'], g['statistic'], g['measure'], 'pooled' if g['pool'] == 'pooled' else 'largest run mean', dv, b['n'], v, b['lo'], b['hi'],
            ' '.join('%+.4f' % p for p in per), g['bound'], 'pass' if ok else 'FAIL'))
        rep_g.append(dict(name=g['name'], measure=g['measure'], statistic=g['statistic'], pool=g['pool'], bound=g['bound'], device_value=dv,
                          device_pooled=b, device_run_means=per, passed=ok))
        if not ok:
            verdict = 'FAIL'
            reasons.append('gate %s: device %+.4f > bound %+.4f' % (g['name'], dv, g['bound']))
    for c in s['controls']:
        if c['kind'] == 'variant':
            group, member = main, c['member']
        else:
            group, member = [x for x in runs if x['role'] == c['role']], 'dev'
        rc = dict(name=c['name'], kind=c['kind'], member=member, role=c.get('role', 'main'), checks=[])
        for gname, must in [(c['trips'], True)] + [(x, False) for x in c.get('must_not_trip', [])]:
            g = gates[gname]
            v, per, allv = values(group, member, g)
            cv = judge(v, per, g, must)
            hit = cv > g['bound']
            ok = hit if must else not hit
            out.append('control %s (%s %s): gate %s %s(%s) %s = %+.4f (runs %s), must be %s %+.4f: %s' % (
                c['name'], c['kind'], member if c['kind'] == 'variant' else c['role'], gname, g['statistic'], g['measure'],
                'pooled' if g['pool'] == 'pooled' else ('smallest' if must else 'largest') + ' run mean', cv, ' '.join('%+.4f' % p for p in per),
                '>' if must else '<=', g['bound'], ('killed' if must else 'not tripped') if ok else ('NOT KILLED' if must else 'TRIPPED')))
            rc['checks'].append(dict(gate=gname, must_trip=must, value=cv, pooled=v, run_means=per, n=len(allv), ok=ok))
            if not ok:
                verdict = 'FAIL'
                reasons.append('control %s %s gate %s (%+.4f %s bound %+.4f)' % (c['name'], 'not killed by' if must else 'trips', gname, cv,
                                                                                '<=' if must else '>', g['bound']))
        rep_c.append(rc)
    uncovered = [g['name'] for g in s['gates'] if not any(c['trips'] == g['name'] for c in s['controls'])]
    if verdict == 'PASS':
        if not s['controls'] or uncovered:
            verdict = 'DEV-ONLY'
            reasons.append('no control declared' if not s['controls'] else 'gates without a control that trips them: %s' % uncovered)
        elif 'test_only' in s:
            verdict = 'TEST-ONLY'
            reasons.append('spec marked test_only (%s): not an acceptance verdict' % s['test_only'])
    return verdict, reasons, out, rep_g, rep_c


def run_report(x, s):
    rr = dict(dir=x['dir'], role=x['role'], qref_sha256=x['qref_sha256'], run_log_sha256=x['run_log_sha256'], seeds_sha256=x['seeds_sha256'],
              parameters_sha256=x['parameters_sha256'], manifest_sha256=x['manifest_sha256'], prog=x['proc']['prog'], prog_src=x['proc']['prog_src'],
              args=x['run'].spec['args'], env=x['run'].spec['run']['env'], n_points=len(x['run'].points), measures={})
    for m, res in x['res'].items():
        rr['measures'][m] = {mem: dict(n=d['n'], S_ens_mean=statistics.fmean(d['re']), S_ens_sd=statistics.pstdev(d['re']),
                                       S_cpu_mean=statistics.fmean(d['rc'])) for mem, d in res.items()}
    return rr


def main(argv):
    rep = dict(tool='t207_d5_gate.py (207-2 design 3)', argv=argv,
               scripts={f: sha(os.path.join(HERE, f)) for f in ('t207_d5_gate.py', 't207_d5_agg.py', 't207_d5_qparse.py', 't207_d5_qref.F90')})
    report = None
    try:
        opts, dirs, i = {}, [], 0
        while i < len(argv):
            a = argv[i]
            if a in ('--spec', '--report'):
                if a in opts or i + 1 >= len(argv):
                    raise InputError('%s given twice or without a value' % a)
                opts[a] = argv[i + 1]
                i += 2
            elif a.startswith('-'):
                raise InputError('unknown option %r' % a)
            else:
                dirs.append(a)
                i += 1
        if '--spec' not in opts or '--report' not in opts or not dirs:
            raise InputError('usage: t207_d5_gate.py --spec SPEC.json --report OUT.json RUNDIR [RUNDIR ...]')
        if os.path.lexists(opts['--report']) or not os.path.isdir(os.path.dirname(os.path.abspath(opts['--report']))):
            raise InputError('report %s exists or its directory does not (reports are never overwritten)' % opts['--report'])
        report = opts['--report']
        rep['spec_path'] = os.path.abspath(opts['--spec'])
        s, table = load_spec(opts['--spec'])
        rep['spec_sha256'] = sha(opts['--spec'])
        rep['spec'] = s
        runs = [load_run(d, s, table) for d in dirs]
        cross_checks(runs, s)
        rep['runs'] = [run_report(x, s) for x in runs]
        for x in runs:
            print('  run %s: role %s, %s' % (x['dir'], x['role'], ', '.join('%s %d points' % (m, x['res'][m]['cpu']['n']) for m in sorted(x['res']))))
        verdict, reasons, lines, rep['gates'], rep['controls'] = decide(runs, s)
        for ln in lines:
            print('  ' + ln)
    except InputError as e:
        verdict, reasons = 'INPUT_ERROR', [str(e)]
    except Exception as e:                       # anything unforeseen (a bug, an unreadable file): not a verdict, never a FAIL
        verdict, reasons = 'INPUT_ERROR', ['%s: %s' % (type(e).__name__, e)]
    code = EXIT[verdict]
    rep.update(verdict=verdict, reasons=reasons, exit_code=code)
    if report is not None:
        try:
            with open(report, 'x') as f:
                json.dump(rep, f, indent=1, allow_nan=False)
        except Exception as e:
            print('GATE INPUT_ERROR: report %s could not be written: %s: %s' % (report, type(e).__name__, e))
            return EXIT['INPUT_ERROR']
    print('GATE %s%s' % (verdict, (': ' + '; '.join(reasons)) if reasons else ''))
    return code


if __name__ == '__main__':
    try:
        rc = main(sys.argv[1:])
    except BaseException as e:                   # nothing escapes as a traceback with another exit status
        print('GATE INPUT_ERROR: %s: %s' % (type(e).__name__, e))
        rc = EXIT['INPUT_ERROR']
    sys.exit(rc)
