#!/usr/bin/env python3
"""Ticket 207-2: accuracy gate against the quad reference, recomputed from the evidence (design 3, fix round 4 after review 207-2d).

usage: t207_d5_gate.py --spec SPEC.json --spec-sha SHA256 --prog QREF_BINARY --report OUT.json RUNDIR [RUNDIR ...]

  --spec-sha  the sha256 of the pre-registered spec (from the acceptance file): a spec file with another sha256 is refused before it is parsed;
  --prog      the offline program t207_d5_qref.F90 built for the spec's path (qref_pp.out / qref_log.out); its sha256 on disk must be the
              spec's prog_sha.  The gate re-runs it on every kept dump of every run (authenticity, step 9 below).

The gate decides only on numbers it computes itself from the run directories; no aggregate is ever an input.  For every RUNDIR (the
directory holding the run's records, qref.txt and the processing records of t207_d5_proc.sh):
  1. qref_spec.json re-derived from the run's own records (parameters, info, run.log, stamp.txt, provenance.json); the run's path
     (from info) is the spec's path;
  2. qref_proc.json (t207_d5_qparse.proc_record): rc 0, check_rc 0, the arguments of qref_spec.json, the program sha256 = spec prog_sha,
     the program source sha256 = spec src_sha, qref_sha256 = the sha256 of the qref.txt bytes the gate parses, manifest_sha256 = the
     sha256 of the manifest bytes the gate parses;
  3. the dump evidence (t207_d5_qparse.verify_dumps): the manifest lists exactly the covered dumps, the kept dumps match it; every kept dump
     is a dump of this run at its place (t207_d5_qparse.check_dump_headers: magic, version 2, file / pass / nst / nfe / br of run.log, N,
     path, QVAR masks, QHOST, TAUDIAG, hasb, exact size as the program checks it; a text file posing as a dump is refused);
  4. its processing arguments qs, nul, cover are the spec's; its role is the row of the env table that equals its ALF_GPU_* env exactly;
  5. qref.txt is loaded by the strict parser (version 1 dumps never accepted, every number checked; the eF invariant is a sanity check
     only) with the spec's member list of the role (exactly those members, in that order);
  6. per gate measure, the per-point log ratios S_ens (and S_cpu, report only) by t207_d5_agg.per_run; the number of points used must equal
     the spec's expect_points for that measure (pTT: points with B, emax / eF: all points).
Across the runs: 7. the same parameters; no run twice (directory, qref.txt sha256, run.log sha256, seeds, a dump sha256 or a QREF data
line shared by two runs); the number of runs per role equals the spec's n_runs.  8. Pinning: every run's ALF binary (provenance.json
binary_sha256 = stamp.txt) is the spec's binary_sha256, its parameters sha256 the spec's parameters_sha256, and the seed sets per role are
exactly the spec's seeds (set equality; no run chosen out of a larger set, no other seed set).  9. Authenticity: the program (--prog) is run
on the kept dumps of every run with the run's processing arguments; its output must equal the run's qref.txt lines of those dumps (headers
except the timing tq exactly, member names and order exactly, every number within +-1 in the last printed digit; the dev5 re-processing
gave only equal numbers or 1-digit differences of the format change).
WHAT THIS GATE DOES NOT PROVE (left to the harness, ticket 207-5): the run directory is under the control of whoever produced it.  The gate
proves that the numbers of the kept dumps are what the pinned program computes from those dumps, that those dumps are consistent with the
run's records, and that the records are internally consistent.  It does not prove that the dumps, run.log, provenance.json, qref.txt lines of
deleted dumps (cover 1, keep 40: 39 of 40 lines are not re-computed), qref_proc.json or the manifest were produced by the stamped run: that
needs the stamps of the records, of qref.txt, of the processing record and of the manifest in the attempt log outside the run directory
(207-5), checked against the files the gate reads.  This gate alone is not forgery-proof.

Code-enforced spec content (not left to the spec writer):
  - env_roles.main is exactly {ALF_GPU_DEBUG=1, ALF_GPU_TAUDIAG_QHOST=1, ALF_GPU_TAUDIAG_QVAR=8, ALF_GPU_TAUDIAG_QREF=<int>=1..,
    ALF_GPU_TAUDIAG_FROM=<int>}; a control role may omit only ALF_GPU_TAUDIAG_QVAR and adds at least one fault from the table FAULTS
    (ALF_GPU_TAU_NS=0 default path, ALF_GPU_TAUSWAP=8 LOG path, ALF_GPU_FAULT_TAU=eps:<x>[@<nst>[:<nfe>]] both paths);
  - refs equal the frozen list of the path (FROZEN_REFS, host member = host_member) unless the spec is marked test_only; var<m>, u<k>,
    dev, n0 and the diagnostic members are never references.
Decision, per gate {name, measure, statistic, pool, bound}:
  statistic: S_ens (mean of the per-point values), S_ens_q95 (their 95th percentile, linear interpolation), S_ens_max (their maximum);
  pool 'pooled'    : the statistic over the per-point values of all runs of the role (every point weighs the same);
  pool 'each_run'  : the statistic of every run on its own (device: the largest run value must be <= bound; control: the smallest run value
                     must be > bound, i.e. every run detects it);
  pool 'by_stratum': the statistic over the points of every stratum = call-cycle position (nst, nfe, br), pooled over the runs of the role
                     (device: every stratum must be <= bound, i.e. the largest stratum value; control: killed if the largest stratum value is
                     > bound, i.e. the gate rejects it; must_not_trip: the largest stratum value <= bound).  The per-stratum point counts are
                     reported.  Which of these become gates and with which bounds is decided in 207-3 (decision 37), not here.
  device (member 'dev' of the main runs): FAIL if its decisive value > bound.
  control (every control of the spec, always required): 'variant' = a member of the main runs (e.g. var8), 'run' = member 'dev' of the
  runs whose env is the control's role; judged on the statistic, pool and bound of the gate named in 'trips' and must exceed the bound
  (killed); for every gate in its optional 'must_not_trip' list it must stay <= that gate's bound.  Otherwise FAIL.
  Always reported (not decisive unless a gate uses them): per measure and per judged member the pooled mean / q95 / max, the per-run values
  and the per-stratum means with their point counts (section 'localized' of the report).
Verdicts and exit codes (printed as the last line 'GATE <VERDICT>...'; the report carries the same verdict and exit code):
  PASS 0        every gate passes, every control behaves as declared, every gate is tripped by at least one control, spec not test_only
  FAIL 1        a gate's own criterion caught something: device above a bound, a control not killed, or tripping a gate it must not trip
  INPUT_ERROR 2 the gate could not be evaluated: a bad or incomplete spec, a spec hash or program hash other than required, a missing or
                inconsistent file, a count differing from the spec, a pin not met, a reprocessing mismatch, any refusal above, or any
                unexpected exception (never a FAIL, never a detection)
  DEV-ONLY 3    the spec declares no control, or a gate no control trips: the numbers are computed and reported, never a PASS
  TEST-ONLY 4   all criteria met, but the spec is marked test_only (e.g. the kept dev5 subset, cover != 1): not an acceptance verdict
The report (--report, must not exist) holds the spec and its sha256 (of the bytes parsed), the sha256 of this script, t207_d5_agg.py,
t207_d5_qparse.py and t207_d5_qref.F90, per run the identities (ALF binary, parameters, seeds, qref.txt, manifest), the authenticity check,
per-member statistics, every decisive number, the verdict and the exit code.  It is written in every case where the report path is usable
(also for INPUT_ERROR) and is never read by the gate.

Spec file (JSON object, no duplicate keys, no NaN / Infinity; every key below required unless marked optional, no other key):
  spec_version 1; path 'def' | 'log'; host_member: the label the offline program prints for the in-run host CGR2_2 of QHOST=1
  (t207_d5_qref.F90 prints 'nvh' whatever the run's compiler); members: {role: [the exact ordered statistical members of that role's runs]}
  for main and the roles of the run controls; refs: the frozen reference list of the path (another list only with test_only);
  qs: one integer > 0; nul: the null codes of the processing (distinct integers >= 0, with 0 and the codes of refs); cover: 1 (another
  value only with test_only); test_only (optional): a non-empty reason; expect_points: {measure: integer >= 1} for exactly the measures of
  the gates; prog_sha, src_sha: sha256 of the offline program and of t207_d5_qref.F90 (src_sha = the file next to this script);
  binary_sha256: the ALF binary of every run (provenance.json); parameters_sha256: the parameters file of every run;
  seeds: {role: [seeds_sha256, ...]} for main and exactly the roles of the run controls, distinct across all roles (control seed sets are
  not main seed sets); env_roles: {'main': {VAR: value, ...}, ROLE: {'extra_env': {...}, 'omit_env': [...]}, ...} (see above);
  n_runs: {'main': N, ROLE: N, ...}; gates: [{name, measure (pTT | emax | eF), statistic (S_ens | S_ens_q95 | S_ens_max),
  pool ('pooled' | 'each_run' | 'by_stratum'), bound}, ...];
  controls: [{name, kind 'variant', member (a member of main, not dev, not a reference), trips, must_not_trip (optional)} |
             {name, kind 'run', role, trips, must_not_trip (optional)}, ...].
  allow_v1 is refused in any form.
"""
import hashlib, importlib.util, json, math, os, re, statistics, subprocess, sys
sys.dont_write_bytecode = True      # the loader below imports a sibling script by path: no __pycache__ in the scripts directory
HERE = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location('t207_d5_agg', os.path.join(HERE, 't207_d5_agg.py'))
AGG = importlib.util.module_from_spec(_s)
_s.loader.exec_module(AGG)
QP = AGG.QP                                  # one instance of the parser module, the one the aggregation uses
QError = QP.QError

EXIT = {'PASS': 0, 'FAIL': 1, 'INPUT_ERROR': 2, 'DEV-ONLY': 3, 'TEST-ONLY': 4}
MEASURES = ('pTT', 'emax', 'eF')             # measures of qref.txt the gate accepts (t207_d5_agg.meas); anything else is refused
STATISTICS = ('S_ens', 'S_ens_q95', 'S_ens_max')     # decision statistics of the per-point S_ens values; S_cpu is reported only
POOLS = ('pooled', 'each_run', 'by_stratum')
SPEC_REQUIRED = {'spec_version', 'path', 'host_member', 'members', 'refs', 'qs', 'nul', 'cover', 'expect_points', 'prog_sha', 'src_sha',
                 'binary_sha256', 'parameters_sha256', 'seeds', 'env_roles', 'n_runs', 'gates', 'controls'}
SPEC_OPTIONAL = {'test_only'}
GATE_KEYS = {'name', 'measure', 'statistic', 'pool', 'bound'}
NAME = re.compile(r'^[A-Za-z0-9_.-]+$')
HEX64 = QP.HEX64
# The frozen reference ensembles (decision 33, ticket 207-2 task 3); 'H' stands for the spec's host_member.  A spec with another list is
# refused unless marked test_only (it can then never PASS).
FROZEN_REFS = {'def': ['cpu', 'H', 'n4', 'n100004', 'n1', 'n100001', 'n2', 'n1608', 'n4808', 'n5', 'n100005'],
               'log': ['cpu', 'H', 'n4', 'n100004', 'n1', 'n100001', 'n2', 'n5', 'n100005']}
NEVER_REF = re.compile(r'^(dev|n0|var[0-9]+|u[0-9]+)$')
# The env of the main runs (the device under test) and the faults a control role may add: fixed here, not in the spec.
MAIN_ENV_FIXED = {'ALF_GPU_DEBUG': '1', 'ALF_GPU_TAUDIAG_QHOST': '1', 'ALF_GPU_TAUDIAG_QVAR': '8'}
MAIN_ENV_INT = {'ALF_GPU_TAUDIAG_QREF': 1, 'ALF_GPU_TAUDIAG_FROM': 0}          # name -> smallest value
CONTROL_MAY_OMIT = ('ALF_GPU_TAUDIAG_QVAR',)
_EPS = r'(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eEdD][-+]?[0-9]+)?'
# fault variable -> (paths, value regex, meaning).  ALF_GPU_FAULT_TAU=eps:<x> is the graded fault of 207-3 (x > 0, < 1e300, as
# gpu_state_mod validates it); '@<nst>[:<nfe>]' restricts it to one block (and flavour) position.  207-3 defines the position syntax: if it
# differs from this, this table is the one place to change (QUESTION in the 207-2 report).
FAULTS = {'ALF_GPU_TAU_NS': (('def',), re.compile(r'^0$'), 'default path forced (no pre-pivot): NS=0 control'),
          'ALF_GPU_TAUSWAP': (('log',), re.compile(r'^8$'), 'LOG path: device CGR2_2 step swapped (mask 8)'),
          'ALF_GPU_FAULT_TAU': (('def', 'log'), re.compile(r'^eps:(' + _EPS + r')(?:@[1-9][0-9]*(?::[1-9][0-9]*)?)?$'),
                                'graded fault of 207-3 on the device CGR2_2 output')}
MEMO = {}       # reprocessing results within ONE process (the self-test calls main() in-process); see reprocess()


class InputError(Exception):
    pass


def isint(x):
    return isinstance(x, int) and not isinstance(x, bool)


def isnum(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def sha(p):
    return QP.sha256_file(p)


def frozen_refs(path, host):
    return [host if r == 'H' else r for r in FROZEN_REFS[path]]


# ------------------------------------------------------------------------------------------------------------------- spec

def load_spec(path, want_sha):
    """the spec from the bytes of the file, refused unless their sha256 is want_sha; returns (spec, env table, sha256)"""
    if not os.path.isfile(path):
        raise InputError('gate spec missing: %s' % path)
    with open(path, 'rb') as f:
        raw = f.read()
    got = hashlib.sha256(raw).hexdigest()
    if got != want_sha:
        raise InputError('gate spec %s: sha256 %s differs from --spec-sha %s (only the pre-registered spec is evaluated)' % (path, got, want_sha))
    try:
        s = QP.loads_strict(raw, 'gate spec', path)
    except QError as e:
        raise InputError(str(e))
    s, table = check_spec(s)
    return s, table, got


def check_spec(s):
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
    if 'test_only' in s and not (isinstance(s['test_only'], str) and s['test_only'].strip()):
        raise InputError('gate spec: test_only must be a non-empty reason')
    refs = s['refs']
    if not (isinstance(refs, list) and len(refs) >= 2 and all(isinstance(r, str) for r in refs)):
        raise InputError('gate spec: refs %r is not a list of at least two names' % (refs,))
    if len(set(refs)) != len(refs):
        raise InputError('gate spec: refs %s: duplicate name' % refs)
    for r in refs:
        if NEVER_REF.match(r) or r in QP.DIAG_NAMES:
            raise InputError('gate spec: refs: %r is not a reference member (dev is the device under test; n0 is the oracle cpu by '
                             'construction and would count twice; var<m> are device variants and u<k> input-sensitivity members, never '
                             'references; %s are diagnostics)' % (r, ', '.join(QP.DIAG_NAMES)))
    fz = frozen_refs(s['path'], s['host_member'])
    if refs != fz and 'test_only' not in s:
        raise InputError('gate spec: refs %s differ from the frozen reference list of the %s path %s (decision 33; host member %s); another '
                         'list only in a spec marked test_only' % (refs, s['path'], fz, s['host_member']))
    if not (isint(s['qs']) and s['qs'] > 0):
        raise InputError('gate spec: qs %r: one integer > 0 expected (one qs per gate run)' % (s['qs'],))
    nul = s['nul']
    if not (isinstance(nul, list) and all(isint(c) and c >= 0 for c in nul) and len(set(nul)) == len(nul) and 0 in nul):
        raise InputError('gate spec: nul %r: a list of distinct integers >= 0 containing 0 expected' % (nul,))
    if not (isint(s['cover']) and s['cover'] >= 1):
        raise InputError('gate spec: cover %r: an integer >= 1 expected' % (s['cover'],))
    if s['cover'] != 1 and 'test_only' not in s:
        raise InputError('gate spec: cover %d: qref.txt must cover every dump of the run (cover 1); another cover only in a spec marked '
                         'test_only' % s['cover'])
    for k in ('prog_sha', 'src_sha', 'binary_sha256', 'parameters_sha256'):
        if not (isinstance(s[k], str) and HEX64.match(s[k])):
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
    # member lists: the exact ordered statistical members of the runs of each role
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
    roles = {'main'} | run_roles
    if set(mem) != roles:
        raise InputError('gate spec: members %s must have exactly the roles %s' % (sorted(mem), sorted(roles)))
    nr = s['n_runs']
    if not isinstance(nr, dict) or set(nr) != roles:
        raise InputError('gate spec: n_runs %r must have exactly the roles %s' % (nr, sorted(roles)))
    for r, n in nr.items():
        if not (isint(n) and n >= 1):
            raise InputError('gate spec: n_runs %s: %r is not an integer >= 1' % (r, n))
    # seeds: the exact seed sets per role, distinct across roles (decision 33: control seed sets are not main seed sets)
    sd = s['seeds']
    if not isinstance(sd, dict) or set(sd) != roles:
        raise InputError('gate spec: seeds %r must be a dictionary with exactly the roles %s' % (sd, sorted(roles)))
    allseeds = {}
    for r, lst in sd.items():
        if not (isinstance(lst, list) and lst and all(isinstance(x, str) and HEX64.match(x) for x in lst) and len(set(lst)) == len(lst)):
            raise InputError('gate spec: seeds.%s %r: a non-empty list of distinct seeds_sha256 values (64 lower-case hex digits) expected' % (r, lst))
        for x in lst:
            if x in allseeds:
                raise InputError('gate spec: seeds: %s is listed for role %s and role %s (control seed sets must differ from main seed sets '
                                 'and from each other, decision 33)' % (x[:16], allseeds[x], r))
            allseeds[x] = r
    check_env_rules(s, table, run_roles)
    return s, table


def check_env_rules(s, table, run_roles):
    """the code-enforced env: main exactly MAIN_ENV_FIXED plus QREF / FROM integers; a control role = main minus at most CONTROL_MAY_OMIT
    plus at least one fault of FAULTS for the spec's path, nothing else"""
    main = s['env_roles']['main']
    want = set(MAIN_ENV_FIXED) | set(MAIN_ENV_INT)
    if set(main) != want or any(main[k] != v for k, v in MAIN_ENV_FIXED.items()) or \
            any(not re.match(r'^[0-9]+$', main[k]) or int(main[k]) < lo for k, lo in MAIN_ENV_INT.items()):
        raise InputError('gate spec: env_roles.main %s must be exactly %s plus ALF_GPU_TAUDIAG_QREF=<integer >= 1> and ALF_GPU_TAUDIAG_FROM='
                         '<integer >= 0> (fixed by the gate: any other variable could change what the device result is)' % (
                             ' '.join('%s=%s' % kv for kv in sorted(main.items())), ' '.join('%s=%s' % kv for kv in sorted(MAIN_ENV_FIXED.items()))))
    for role in sorted(run_roles):
        r = s['env_roles'][role]
        bad_omit = [k for k in r['omit_env'] if k not in CONTROL_MAY_OMIT]
        if bad_omit:
            raise InputError('gate spec: env_roles.%s.omit_env %s: a control may omit only %s' % (role, bad_omit, list(CONTROL_MAY_OMIT)))
        if not r['extra_env']:
            raise InputError('gate spec: env_roles.%s adds no fault (extra_env empty): a control role adds at least one fault of %s' % (role, sorted(FAULTS)))
        for k, v in sorted(r['extra_env'].items()):
            f = FAULTS.get(k)
            ok = f is not None and s['path'] in f[0] and f[1].match(v)
            if ok and k == 'ALF_GPU_FAULT_TAU':
                x = float(f[1].match(v).group(1).replace('d', 'e').replace('D', 'e'))
                ok = 0 < x < 1e300
            if not ok:
                raise InputError('gate spec: env_roles.%s.extra_env %s=%s is not a fault of the fixed list for the %s path (%s)' % (
                    role, k, v, s['path'], '; '.join('%s %s (%s)' % (n, t[1].pattern, '/'.join(t[0])) for n, t in sorted(FAULTS.items()))))


# ------------------------------------------------------------------------------------------------------------------- runs

def load_run(rundir, s, table):
    d = os.path.abspath(rundir)
    what = 'run %s' % d
    if not os.path.isdir(d):
        raise InputError('%s: not a directory' % what)
    qt = os.path.join(d, 'qref.txt')
    if not os.path.isfile(qt):
        raise InputError('%s: qref.txt missing' % what)
    try:
        rspec = QP.load_spec(d)                      # re-derived from the run's records and compared with the stored one
        if rspec['run']['variant'] != s['path'] or rspec['args']['path'] != s['path']:
            raise InputError('%s: a %s-path run (info), processed as %s; the spec is for the %s path' % (
                what, rspec['run']['variant'], rspec['args']['path'], s['path']))
        with open(qt, 'rb') as f:
            raw = f.read()                           # the bytes hashed here are the bytes parsed below (qref.txt is not read again)
        qsha = hashlib.sha256(raw).hexdigest()
        mparsed = QP.read_manifest_sha(d)            # likewise the manifest
        pr = QP.proc_record(d, rspec, prog=s['prog_sha'], src=s['src_sha'], qref_sha256=qsha, manifest_sha256=mparsed[1])
        manifest, msha = QP.verify_dumps(d, rspec, parsed=mparsed)
        ndh = QP.check_dump_headers(d, rspec)
        for k in ('qs', 'nul', 'cover'):
            if rspec['args'][k] != s[k]:
                raise InputError('%s: processed with %s %r, the spec says %r' % (what, k, rspec['args'][k], s[k]))
        role = QP.role_of_env(rspec['run']['env'], table)
        run = QP.load(qt, allow_v1=False, members=s['members'][role], host=s['host_member'], raw=raw)
        r = run.spec['run']
        res = {}
        for m in sorted(set(g['measure'] for g in s['gates'])):
            res[m], n_exp = AGG.per_run(run, m, s['refs'])
            n = res[m]['cpu']['n']
            if not (isint(n) and n == n_exp == s['expect_points'][m]):
                raise InputError('%s: %r points used for measure %s (derived from the run %r), the spec expects %d' % (what, n, m, n_exp, s['expect_points'][m]))
    except QError as e:
        raise InputError('%s: %s' % (what, e))
    return dict(dir=d, real=os.path.realpath(d), role=role, run=run, res=res, proc=pr, manifest=manifest, lines=run.lines, kept=QP.kept_dumps(rspec),
                n_headers=ndh, qref_sha256=qsha, run_log_sha256=r['records_sha256']['run.log'], seeds_sha256=r['seeds_sha256'],
                parameters_sha256=r['parameters_sha256'], binary_sha256=r['binary_sha256'], manifest_sha256=msha)


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


def pins(runs, s):
    """the runs are the pre-registered ones: ALF binary, parameters, and per role exactly the spec's seed sets"""
    for x in runs:
        if x['binary_sha256'] != s['binary_sha256']:
            raise InputError('run %s: ALF binary sha256 %s (provenance.json), the spec requires binary_sha256 %s' % (x['dir'], x['binary_sha256'], s['binary_sha256']))
        if x['parameters_sha256'] != s['parameters_sha256']:
            raise InputError('run %s: parameters sha256 %s, the spec requires parameters_sha256 %s' % (x['dir'], x['parameters_sha256'], s['parameters_sha256']))
    for role, want in sorted(s['seeds'].items()):
        have = sorted(x['seeds_sha256'] for x in runs if x['role'] == role)
        if set(have) != set(want) or len(have) != len(want):
            raise InputError('role %s: the runs have the seed sets %s, the spec requires exactly %s (no run may be left out, added or replaced)'
                             % (role, [h[:16] for h in have], sorted(w[:16] for w in want)))


# ------------------------------------------------------------------------------------------------------------------- authenticity

NUMTOK = re.compile(r'^(-?)([0-9])\.([0-9]{4})E([+-][0-9]{3})$')


def last_digit_diff(a, b):
    """|a - b| in units of the last printed digit (ES12.4E3) of the finer of the two numbers; None if a token is not such a number"""
    ma, mb = NUMTOK.match(a), NUMTOK.match(b)
    if not ma or not mb:
        return None
    va = (-1 if ma.group(1) else 1) * int(ma.group(2) + ma.group(3))
    vb = (-1 if mb.group(1) else 1) * int(mb.group(2) + mb.group(3))
    ea, eb = int(ma.group(4)), int(mb.group(4))
    e0 = min(ea, eb)
    if abs(ea - eb) > 2:
        return float('inf')
    return abs(va * 10 ** (ea - e0) - vb * 10 ** (eb - e0))


def compare_line(old, new, what):
    """(numbers, equal, one_digit) of the run's qref.txt line old against the re-computed line new; InputError on any other difference"""
    so, sn = old.split(' | '), new.split(' | ')
    ho, hn = so[0].split(), sn[0].split()
    if len(ho) != len(hn) or 'tq' not in ho or ho.index('tq') != len(ho) - 2 or hn.index('tq') != len(hn) - 2 or ho[:-1] != hn[:-1]:
        raise InputError('%s: header %r, re-computed %r (tq aside, they must be equal)' % (what, so[0], sn[0]))
    if len(so) != len(sn):
        raise InputError('%s: %d members in qref.txt, %d re-computed' % (what, len(so) - 1, len(sn) - 1))
    nnum = neq = n1 = 0
    for a, b in zip(so[1:], sn[1:]):
        ta, tb = a.split(), b.split()
        if not ta or not tb or ta[0] != tb[0] or len(ta) != len(tb):
            raise InputError('%s: member %r in qref.txt where the program prints %r' % (what, a[:40], b[:40]))
        for x, y in zip(ta[1:], tb[1:]):
            nnum += 1
            if x == y:
                neq += 1
                continue
            dd = last_digit_diff(x, y)
            if dd is None or dd > 1:
                raise InputError('%s: member %s: %s in qref.txt, %s re-computed by the program (allowed: +-1 in the last printed digit)'
                                 % (what, ta[0], x, y))
            n1 += 1
    return nnum, neq, n1


def program_env():
    e = dict(os.environ)
    e.update(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='1')      # as t207_d5_proc.sh (no reduction depends on the thread count)
    return e


def reprocess(x, s, prog):
    """the program's output lines for the kept dumps of run x.  The MEMO key is the program sha256 verified in this call, the processing
    arguments and the (name, sha256) of every kept dump, verified against the manifest in this call: a result is re-used only for the same
    program on the same bytes with the same arguments, and only within one process (one CLI call never meets the same dumps twice; the
    self-test's in-process calls do)."""
    a = x['run'].spec['args']
    names = [QP.dump_name(k) for (k, p, nst, nfe, br) in x['kept']]
    key = (s['prog_sha'], a['qs'], tuple(a['nul']), a['hostin'], tuple((nm, x['manifest'][nm]) for nm in names))
    if key in MEMO:
        return MEMO[key], 'memo'
    nul = ','.join(str(c) for c in a['nul']) if a['nul'] else '-'
    try:
        p = subprocess.run([prog, str(a['qs']), nul, str(a['hostin'])], input=''.join(nm + '\n' for nm in names), cwd=x['dir'],
                           capture_output=True, text=True, env=program_env(), timeout=24 * 3600)
    except (OSError, subprocess.SubprocessError) as e:
        raise InputError('run %s: the program %s could not be run: %s: %s' % (x['dir'], prog, type(e).__name__, e))
    lines = p.stdout.split('\n')
    if p.returncode != 0 or not p.stdout.endswith('\n'):
        raise InputError('run %s: re-processing of the kept dumps: program exit status %d (%s)' % (
            x['dir'], p.returncode, ' / '.join(l[:200] for l in lines if l.startswith('QREFERR'))[:600] or p.stderr.strip()[-300:]))
    lines = lines[:-1]
    want_end = 'QREFEND dumps %d ok %d err 0' % (len(names), len(names))
    if len(lines) != len(names) + 1 or lines[-1] != want_end:
        raise InputError('run %s: re-processing printed %d lines ending %r, expected %d QREF lines and %r' % (
            x['dir'], len(lines), lines[-1][:80] if lines else '', len(names), want_end))
    MEMO[key] = lines[:-1]
    return lines[:-1], 'program'


def authenticate(runs, s, prog):
    """step 9: the program re-run on the kept dumps of every run; the run's qref.txt lines of those dumps must be what it prints"""
    out = []
    for x in runs:
        new, src = reprocess(x, s, prog)
        bynum = {}
        for ln in x['lines']:
            bynum[int(ln.split()[2])] = ln
        tot = [0, 0, 0]
        for (k, p, nst, nfe, br), ln in zip(x['kept'], new):
            if not ln.startswith('QREF file %d ' % k):
                raise InputError('run %s: re-processing line for dump %d is %r' % (x['dir'], k, ln[:60]))
            if k not in bynum:
                raise InputError('run %s: kept dump %d has no line in qref.txt' % (x['dir'], k))
            c = compare_line(bynum[k], ln, 'run %s dump %d (re-processing check)' % (x['dir'], k))
            tot = [t + v for t, v in zip(tot, c)]
        out.append(dict(dir=x['dir'], dumps=len(new), numbers=tot[0], equal=tot[1], one_digit=tot[2], source=src))
    return out


# ------------------------------------------------------------------------------------------------------------------- decision

def stat_of(name, v):
    if not v:
        raise InputError('statistic %s of an empty set' % name)
    if name == 'S_ens':
        return statistics.fmean(v)
    if name == 'S_ens_q95':
        return AGG.q(v, 0.95)
    if name == 'S_ens_max':
        return max(v)
    raise InputError('unknown statistic %s' % name)


def points_of(group, member, measure):
    """[(run index, (nst, nfe, br), S_ens value)] of member over the runs of group"""
    out = []
    for i, x in enumerate(group):
        if member not in x['res'][measure]:
            raise InputError('run %s (%s): member %s missing (members: %s)' % (x['dir'], x['role'], member, ' '.join(sorted(x['res'][measure]))))
        d = x['res'][measure][member]
        out += [(i, tuple(p), v) for p, v in zip(d['pos'], d['re'])]
    return out


def values(group, member, g):
    """pooled value, per-run values, per-stratum {pos: (n, value)} of the gate's statistic, and all values"""
    pts = points_of(group, member, g['measure'])
    allv = [v for i, p, v in pts]
    st = g['statistic']
    pooled = stat_of(st, allv)
    per = [stat_of(st, [v for i, p, v in pts if i == j]) for j in range(len(group))]
    strata = {}
    for i, p, v in pts:
        strata.setdefault(p, []).append(v)
    strata = {p: (len(v), stat_of(st, v)) for p, v in sorted(strata.items())}
    for v in [pooled] + per + [t[1] for t in strata.values()]:
        if not math.isfinite(v):
            raise InputError('%s %s: statistic not finite' % (member, g['name']))
    return pooled, per, strata, allv


def judge(pooled, per, strata, g, as_control):
    """the decisive value (see the doc string): pooled; each_run: smallest run value for a control that must trip, else the largest;
    by_stratum: the largest stratum value"""
    if g['pool'] == 'pooled':
        return pooled
    if g['pool'] == 'each_run':
        return min(per) if as_control else max(per)
    return max(t[1] for t in strata.values())


def pool_text(g, must):
    if g['pool'] == 'pooled':
        return 'pooled'
    if g['pool'] == 'each_run':
        return ('smallest' if must else 'largest') + ' run value'
    return 'largest stratum value'


def strata_text(strata):
    return ' '.join('(%d,%d,%d):%d:%+.3f' % (p + (n, v)) for p, (n, v) in strata.items())


def decide(runs, s):
    main = [x for x in runs if x['role'] == 'main']
    gates = {g['name']: g for g in s['gates']}
    out, reasons, rep_g, rep_c = [], [], [], []
    verdict = 'PASS'
    for g in s['gates']:
        v, per, strata, allv = values(main, 'dev', g)
        dv = judge(v, per, strata, g, False)
        b = AGG.summ(allv); b.pop('boot')
        ok = dv <= g['bound']
        out.append('gate %s: device %s(%s) %s = %+.4f (pooled %+.4f over %d points, mean %+.4f [%+.4f, %+.4f]; runs %s; %d strata, points '
                   'per stratum %d..%d), bound %+.4f: %s' % (
                       g['name'], g['statistic'], g['measure'], pool_text(g, False), dv, v, b['n'], b['mean'], b['lo'], b['hi'],
                       ' '.join('%+.4f' % p for p in per), len(strata), min(t[0] for t in strata.values()), max(t[0] for t in strata.values()),
                       g['bound'], 'pass' if ok else 'FAIL'))
        if g['pool'] == 'by_stratum':
            out.append('  strata (nst,nfe,br):points:value  ' + strata_text(strata))
        rep_g.append(dict(name=g['name'], measure=g['measure'], statistic=g['statistic'], pool=g['pool'], bound=g['bound'], device_value=dv,
                          device_pooled_value=v, device_pooled=b, device_run_values=per, device_run_means=per if g['statistic'] == 'S_ens' else None,
                          device_strata=[dict(nst=p[0], nfe=p[1], br=p[2], n=n, value=val) for p, (n, val) in strata.items()], passed=ok))
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
            v, per, strata, allv = values(group, member, g)
            cv = judge(v, per, strata, g, must)
            hit = cv > g['bound']
            ok = hit if must else not hit
            out.append('control %s (%s %s): gate %s %s(%s) %s = %+.4f (runs %s), must be %s %+.4f: %s' % (
                c['name'], c['kind'], member if c['kind'] == 'variant' else c['role'], gname, g['statistic'], g['measure'],
                pool_text(g, must), cv, ' '.join('%+.4f' % p for p in per),
                '>' if must else '<=', g['bound'], ('killed' if must else 'not tripped') if ok else ('NOT KILLED' if must else 'TRIPPED')))
            rc['checks'].append(dict(gate=gname, must_trip=must, value=cv, pooled=v, run_values=per, run_means=per if g['statistic'] == 'S_ens' else None,
                                     strata=[dict(nst=p[0], nfe=p[1], br=p[2], n=n, value=val) for p, (n, val) in strata.items()], n=len(allv), ok=ok))
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


def localized(runs, s):
    """the localized statistics (fix round 4), always reported, decisive only where a gate uses them: per measure of the gates and per judged
    member (dev of main, every variant control member, dev of every run-control role): pooled mean / q95 / max, per-run mean / q95 / max, and
    per stratum (nst, nfe, br) the number of points, mean, q95, max"""
    main = [x for x in runs if x['role'] == 'main']
    who = [('main', 'dev', main)] + [('main', c['member'], main) for c in s['controls'] if c['kind'] == 'variant'] + \
          [(c['role'], 'dev', [x for x in runs if x['role'] == c['role']]) for c in s['controls'] if c['kind'] == 'run']
    out, lines = {}, []
    for m in sorted(set(g['measure'] for g in s['gates'])):
        out[m] = []
        for role, member, group in who:
            pts = points_of(group, member, m)
            allv = [v for i, p, v in pts]
            three = lambda v: dict(n=len(v), mean=statistics.fmean(v), q95=AGG.q(v, 0.95), max=max(v))       # noqa: E731
            st = {}
            for i, p, v in pts:
                st.setdefault(p, []).append(v)
            e = dict(role=role, member=member, pooled=three(allv), runs=[three([v for i, p, v in pts if i == j]) for j in range(len(group))],
                     strata=[dict(nst=p[0], nfe=p[1], br=p[2], **three(v)) for p, v in sorted(st.items())])
            out[m].append(e)
            sm = max(e['strata'], key=lambda z: z['mean'])
            lines.append('localized %s %s %s: pooled mean %+.4f q95 %+.4f max %+.4f; %d strata (points %d..%d), largest stratum mean %+.4f at '
                         '(nst %d, nfe %d, br %d, %d points)' % (m, role, member, e['pooled']['mean'], e['pooled']['q95'], e['pooled']['max'],
                                                                len(e['strata']), min(z['n'] for z in e['strata']), max(z['n'] for z in e['strata']),
                                                                sm['mean'], sm['nst'], sm['nfe'], sm['br'], sm['n']))
    return out, lines


def run_report(x, s):
    rr = dict(dir=x['dir'], role=x['role'], binary_sha256=x['binary_sha256'], parameters_sha256=x['parameters_sha256'], seeds_sha256=x['seeds_sha256'],
              qref_sha256=x['qref_sha256'], run_log_sha256=x['run_log_sha256'], manifest_sha256=x['manifest_sha256'], prog=x['proc']['prog'],
              prog_src=x['proc']['prog_src'], args=x['run'].spec['args'], env=x['run'].spec['run']['env'], n_points=len(x['run'].points),
              kept_dumps_header_checked=x['n_headers'], measures={})
    for m, res in x['res'].items():
        rr['measures'][m] = {mem: dict(n=d['n'], S_ens_mean=statistics.fmean(d['re']), S_ens_sd=statistics.pstdev(d['re']),
                                       S_cpu_mean=statistics.fmean(d['rc'])) for mem, d in res.items()}
    return rr


LIMITS = ('This gate re-computes the numbers of the kept dumps with the pinned program and checks the run records for consistency and against '
          'the pinned binary, parameters and seed sets.  It does not prove that the run directory (dumps, records, qref.txt lines of deleted '
          'dumps, qref_proc.json, manifest) was produced by the stamped run: that is the harness (207-5: stamps of the records, qref.txt, '
          'qref_proc.json and qref_manifest.txt in the attempt log outside the run directory).  Not forgery-proof on its own.')


def main(argv):
    rep = dict(tool='t207_d5_gate.py (207-2 design 3, fix round 4)', argv=argv, limits=LIMITS,
               scripts={f: sha(os.path.join(HERE, f)) for f in ('t207_d5_gate.py', 't207_d5_agg.py', 't207_d5_qparse.py', 't207_d5_qref.F90')})
    report = None
    try:
        opts, dirs, i = {}, [], 0
        while i < len(argv):
            a = argv[i]
            if a in ('--spec', '--spec-sha', '--prog', '--report'):
                if a in opts or i + 1 >= len(argv):
                    raise InputError('%s given twice or without a value' % a)
                opts[a] = argv[i + 1]
                i += 2
            elif a.startswith('-'):
                raise InputError('unknown option %r' % a)
            else:
                dirs.append(a)
                i += 1
        if any(k not in opts for k in ('--spec', '--spec-sha', '--prog', '--report')) or not dirs:
            raise InputError('usage: t207_d5_gate.py --spec SPEC.json --spec-sha SHA256 --prog QREF_BINARY --report OUT.json RUNDIR [RUNDIR ...]')
        if os.path.lexists(opts['--report']) or not os.path.isdir(os.path.dirname(os.path.abspath(opts['--report']))):
            raise InputError('report %s exists or its directory does not (reports are never overwritten)' % opts['--report'])
        report = opts['--report']
        rep['spec_path'] = os.path.abspath(opts['--spec'])
        if not HEX64.match(opts['--spec-sha']):
            raise InputError('--spec-sha %r is not a sha256 (64 lower-case hex digits)' % opts['--spec-sha'])
        s, table, rep['spec_sha256'] = load_spec(opts['--spec'], opts['--spec-sha'])
        rep['spec'] = s
        prog = os.path.abspath(opts['--prog'])
        if not (os.path.isfile(prog) and os.access(prog, os.X_OK)):
            raise InputError('program %s: not an executable file' % prog)
        psha = sha(prog)
        rep['prog'] = dict(path=prog, sha256=psha)
        if psha != s['prog_sha']:
            raise InputError('program %s: sha256 %s differs from the required prog_sha %s of the spec' % (prog, psha, s['prog_sha']))
        runs = [load_run(d, s, table) for d in dirs]
        cross_checks(runs, s)
        pins(runs, s)
        rep['runs'] = [run_report(x, s) for x in runs]
        rep['authenticity'] = authenticate(runs, s, prog)
        for x, au in zip(runs, rep['authenticity']):
            print('  run %s: role %s, binary %s, %s; re-processing: %d kept dumps, %d numbers, %d equal, %d within 1 in the last digit (%s)' % (
                x['dir'], x['role'], x['binary_sha256'][:16], ', '.join('%s %d points' % (m, x['res'][m]['cpu']['n']) for m in sorted(x['res'])),
                au['dumps'], au['numbers'], au['equal'], au['one_digit'], au['source']))
        rep['localized'], llines = localized(runs, s)
        verdict, reasons, lines, rep['gates'], rep['controls'] = decide(runs, s)
        for ln in lines + llines:
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
