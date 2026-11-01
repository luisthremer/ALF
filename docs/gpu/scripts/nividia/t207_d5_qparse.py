#!/usr/bin/env python3
"""Ticket 207-2: strict reader of the output of the offline quad-reference program (t207_d5_qref.F90) and of the records of
the run that produced the dumps.  Nothing is skipped silently: every deviation from the expected content is a QError.

Library (callers load it by file path and run python3 -I):
    load(qref_txt, allow_v1=False, raw=None) -> Run   the validated points of RUNDIR/qref.txt (RUNDIR = its directory); raw: the bytes
                                               of qref.txt already read by the caller (the bytes that are hashed are the bytes parsed);
                                               Run.qref_sha256 = sha256 of the parsed bytes, Run.lines = its 'QREF file' lines
    proc_record(rundir, spec, prog, src, qref_sha256=None, manifest_sha256=None)
                                               the processing record qref_proc.json checked (rc, check_rc, args, program and source hash,
                                               qref_sha256 = the hash of qref.txt (given: the hash of the bytes the caller parsed, else the
                                               file now), manifest_sha256 a sha256 (= the given hash of the manifest bytes the caller parsed));
                                               used by thin, t207_d5_agg.py and t207_d5_gate.py
    verify_dumps(rundir, spec)                 manifest = the covered dumps, the kept dumps match it (returns manifest, its sha256)
    check_dump_headers(rundir, spec)           every kept dump: header and exact size as t207_d5_qref.F90 checks them, and the header's file /
                                               pass / nst / nfe / br against run.log, N, path, QVAR masks, QHOST, TAUDIAG, hasb against the
                                               run's records (refuses text files or other data posing as dumps)
    env_table(env_roles), role_of_env(env, table)   the env allow-list as a table of exact envs (driven by the gate spec);
                                               dev_env_table / env_role: the development table of t207_d5_agg.py
    read_json(path, what)                      strict JSON (no duplicate keys, no NaN / Infinity)
CLI (python3 -I t207_d5_qparse.py ...):
    spec  RUNDIR --path def|log --qs Q --nul CODES|- --hostin 0|1 [--cover C] [--keep K]
                          derive the expected content from the run's own records, check the dump directory (before the
                          processing), write RUNDIR/qref_spec.json
    check RUNDIR [--allow-v1] [--dumps]      validate RUNDIR/qref.txt against the spec re-derived from the records
                          (--dumps: and the dump directory after the processing)
    thin  RUNDIR [--allow-v1]   delete the dumps that are not every (cover*keep)-th, only those listed in RUNDIR/qref_manifest.txt
                          whose sha256 still matches and that are not symbolic links; requires qref.txt to load strictly and
                          qref_proc.json with rc 0 and check_rc 0 for the same arguments; nothing is deleted if one check fails
Keep rule: with keep > 0, gcd(cover*keep, P) = 1 is required (P = C/gcd(C, stride), C = NSTM*N_FL), so that the kept dumps visit every
block/flavour position of the call cycle equally often.  The kept dumps allow pointwise re-checking only; the recomputable result is the
full qref.txt plus the manifest.  'check' prints the per-nst / br / nfe counts of the kept dumps.

The expected content is derived from  parameters (Nbin, Nsweep, Nwrap, Beta, Dtau), info (Sweeps, Bins, Nstm, N_FL, the
GPU debug env line, the GPU tau pass line), stamp.txt and provenance.json (env of the run, return code, hashes of the
records), run.log (the TAUQREF line of every dump), and the arguments of the processing (path, QS, null codes, HOSTIN,
cover, keep) stored in qref_spec.json:
  - passes = Nbin*Nsweep, NSTM = ceil(Ltrot/Nwrap); dumps from pass >= FROM, one per STRIDE-th call of qref_post, calls
    ordered pass, nst, nfe; numbered 1.. in that order (qref/q<number>.bin); B, Bi exist for nst < NSTM;
  - qref.txt has exactly the dumps with (number-1) mod cover == 0 (cover 1: all; cover K: a directory that was thinned to
    every K-th dump before this processing), in order, each as one 'QREF file ...' line with the members of the expected
    ordered list and the expected count of numbers (8; 5 for qself, qbr), each number in the ES12.4E3 format and finite,
    errors >= 0, propagated numbers exactly -1 where the dump has no B, n0 identical to cpu, br == hbr == the br of run.log,
    quad self-checks <= QSELF_MAX; and a final line 'QREFEND dumps D ok D err 0'; any other line is an error.
"""
import argparse, hashlib, json, math, os, re, struct, sys
from math import gcd

QSELF_MAX = 1e-25          # input check of the proposed criterion (207 dev5): quad vs itself (row-permuted QR input, other branch)
NUMRE = re.compile(r'^-?[0-9]\.[0-9]{4}E[+-][0-9]{3}$')
REF_DEF = 'cpu,nvh,n4,n100004,n1,n100001,n2,n1608,n4808,n5,n100005'
REF_LOG = 'cpu,nvh,n4,n100004,n1,n100001,n2,n5,n100005'
DIAG_NAMES = ('qself', 'qbr', 'din', 'cpuh', 'cpuhx')     # members that are not error statistics of a G
SPEC_VERSION = 3            # 3 (fix round 4): the run part carries binary_sha256 (the ALF binary of provenance.json = stamp.txt)
DUMP_MAGIC = 20791         # hdr(1) of a dump (Prog/gpu_tau_mod.F90 qref_post; t207_d5_qref.F90 refuses anything else)
# Allow-list of the GPU debug environment (ALF_GPU_*) of a run that enters the accuracy statistics: a table {role: exact env}.  A run has a
# role only if its env equals one row exactly (no extra, no missing variable, same values); a run matching no row, or more than one, is
# refused.  The gate (t207_d5_gate.py) builds the table from its spec file (env_table: 'main' exact, every control role = main minus
# omit_env plus extra_env).  The development aggregation (t207_d5_agg.py) uses dev_env_table(): main = DEV_MAIN_ENV plus the run's own
# QREF / FROM numbers, control = the same plus TAU_NS=0 (default path) or TAUSWAP=8 (LOG path), with or without QVAR (the dev5 control runs
# c_*_dev1 have none).  Any other ALF_GPU_* variable (for instance ALF_GPU_TAUSWAP=<mask>, which makes the 'device' result a host result)
# means the run does not test the device: refused.
DEV_MAIN_ENV = {'ALF_GPU_DEBUG': '1', 'ALF_GPU_TAUDIAG_QHOST': '1', 'ALF_GPU_TAUDIAG_QVAR': '8'}
DEV_COUNT_KEYS = ('ALF_GPU_TAUDIAG_QREF', 'ALF_GPU_TAUDIAG_FROM')       # non-negative integers (QREF >= 1)
DEV_CTL_EXTRA = {'def': ('ALF_GPU_TAU_NS', '0'), 'log': ('ALF_GPU_TAUSWAP', '8')}
ENVNAME = re.compile(r'^ALF_GPU_[A-Z0-9_]+$')
RECORDS = ('parameters', 'info', 'run.log', 'stamp.txt', 'provenance.json')


class QError(Exception):
    pass


class Point(object):
    __slots__ = ('hd', 'm', 'names', 'hasb')

    def __init__(self, hd, m, names):
        self.hd, self.m, self.names, self.hasb = hd, m, names, hd['hasb'] == 1


class Run(object):
    def __init__(self, spec, points, dumps, qref_sha256=None, lines=None):
        self.spec, self.points, self.dumps = spec, points, dumps
        self.qref_sha256, self.lines = qref_sha256, lines          # sha256 of the bytes parsed; the 'QREF file' lines of those bytes

    @property
    def n_with_b(self):
        return sum(1 for p in self.points if p.hasb)


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 20), b''):
            h.update(blk)
    return h.hexdigest()


def read_text(rundir, name):
    p = os.path.join(rundir, name)
    if not os.path.isfile(p):
        raise QError('record %s missing in %s' % (name, rundir))
    with open(p, 'rb') as f:
        raw = f.read()
    try:
        return raw.decode('utf-8')
    except UnicodeDecodeError:
        raise QError('record %s is not valid UTF-8' % name)


def _no_dup_keys(pairs):
    d = {}
    for k, v in pairs:
        if k in d:
            raise QError('duplicate key %r' % (k,))
        d[k] = v
    return d


def _no_const(c):
    raise QError('non-finite constant %s' % c)


def read_json(p, what):
    """a JSON file, strictly: must exist, parse, have no duplicate keys and no NaN / Infinity constants"""
    if not os.path.isfile(p):
        raise QError('%s missing: %s' % (what, p))
    with open(p, 'rb') as f:
        return loads_strict(f.read(), what, p)


def loads_strict(raw, what, p):
    """the JSON value of the bytes raw (read and hashed by the caller), strictly as read_json"""
    try:
        return json.loads(raw.decode('utf-8'), object_pairs_hook=_no_dup_keys, parse_constant=_no_const)
    except QError as e:
        raise QError('%s %s: %s' % (what, p, e))
    except (ValueError, UnicodeDecodeError) as e:
        raise QError('%s %s is not valid JSON: %s' % (what, p, e))


def one(rx, txt, what, flags=re.M):
    m = rx.findall(txt) if hasattr(rx, 'findall') else re.findall(rx, txt, flags)
    if len(m) != 1:
        raise QError('%s: expected exactly one match, found %d' % (what, len(m)))
    return m[0]


def fnum(s):
    return float(s.replace('d', 'e').replace('D', 'e'))


ENVKEY = re.compile(r'ALF_GPU_[A-Z0-9_]+=')


def parse_env_info(line):
    """'ALF_GPU_A=1, ALF_GPU_B=8,16, ...' (values may contain commas): split at the keys"""
    pos = [(m.start(), m.end()) for m in ENVKEY.finditer(line)]
    if not pos or pos[0][0] != 0:
        raise QError('GPU debug env line not parsable: %r' % line)
    d = {}
    for i, (s, e) in enumerate(pos):
        nxt = pos[i + 1][0] if i + 1 < len(pos) else len(line)
        val = line[e:nxt]
        if i + 1 < len(pos):
            if not val.endswith(', '):
                raise QError('GPU debug env line not parsable: %r' % line)
            val = val[:-2]
        k = line[s:e - 1]
        if k in d:
            raise QError('duplicate %s in the GPU debug env line' % k)
        d[k] = val.strip()
    return d


def parse_env_tokens(tokens, what):
    d = {}
    for t in tokens:
        if not ENVKEY.match(t):
            raise QError('%s: not an ALF_GPU_ assignment: %r' % (what, t))
        k, v = t.split('=', 1)
        if k in d:
            raise QError('%s: duplicate %s' % (what, k))
        d[k] = v
    return d


def isint(x):
    return isinstance(x, int) and not isinstance(x, bool)


def check_env_dict(d, what):
    """a dict of ALF_GPU_* names to non-empty string values"""
    if not isinstance(d, dict):
        raise QError('%s is not a dictionary' % what)
    for k, v in d.items():
        if not (isinstance(k, str) and ENVNAME.match(k)):
            raise QError('%s: %r is not an ALF_GPU_* variable name' % (what, k))
        if not isinstance(v, str) or v == '' or v != v.strip():
            raise QError('%s: value of %s is %r (a non-empty string without surrounding blanks expected)' % (what, k, v))
    return d


def env_table(env_roles):
    """the allow-list table {role: exact env} of a gate spec's 'env_roles': {'main': {VAR: value, ...}, ROLE: {'extra_env': {...},
    'omit_env': [VAR, ...]}, ...}.  A control role is main minus omit_env plus extra_env; extra_env names no variable of main (no silent
    override), omit_env names only variables of main; no two roles have the same env."""
    if not isinstance(env_roles, dict) or 'main' not in env_roles:
        raise QError('env_roles: a dictionary with the role main is required')
    main = dict(check_env_dict(env_roles['main'], 'env_roles.main'))
    for k in DEV_COUNT_KEYS:
        if k not in main or not re.match(r'^[0-9]+$', main[k]) or (k.endswith('QREF') and int(main[k]) < 1):
            raise QError('env_roles.main: %s missing or not a valid number (%r)' % (k, main.get(k)))
    table = {'main': main}
    for role, r in env_roles.items():
        if role == 'main':
            continue
        if not (isinstance(role, str) and re.match(r'^[A-Za-z0-9_.-]+$', role)):
            raise QError('env_roles: bad role name %r' % (role,))
        if not isinstance(r, dict) or set(r) != {'extra_env', 'omit_env'}:
            raise QError('env_roles.%s: exactly the keys extra_env and omit_env are required' % role)
        extra = check_env_dict(r['extra_env'], 'env_roles.%s.extra_env' % role)
        omit = r['omit_env']
        if not isinstance(omit, list) or not all(isinstance(k, str) and k in main for k in omit) or len(set(omit)) != len(omit):
            raise QError('env_roles.%s.omit_env %r: a list of distinct variables of main expected' % (role, omit))
        if not extra and not omit:
            raise QError('env_roles.%s: extra_env and omit_env both empty (the env of main)' % role)
        for k in extra:
            if k in main:
                raise QError('env_roles.%s.extra_env: %s is a variable of main (no override: omit it and add it)' % (role, k))
        e = {k: v for k, v in main.items() if k not in omit}
        e.update(extra)
        table[role] = e
    roles = sorted(table)
    for i, a in enumerate(roles):
        for b in roles[i + 1:]:
            if table[a] == table[b]:
                raise QError('env_roles: roles %s and %s have the same env' % (a, b))
    return table


def role_of_env(env, table):
    """the role whose env equals the run's env exactly, else QError"""
    if not isinstance(env, dict):
        raise QError('env is not a dictionary')
    hit = [role for role, want in table.items() if env == want]
    if len(hit) != 1:
        raise QError('env allow-list: the run has ALF_GPU_* %s, which is %s; the allowed environments are: %s' % (
            ' '.join('%s=%s' % kv for kv in sorted(env.items())), 'not an allowed one' if not hit else 'more than one role %s' % hit,
            '; '.join('%s: %s' % (r, ' '.join('%s=%s' % kv for kv in sorted(t.items()))) for r, t in sorted(table.items()))))
    return hit[0]


def dev_env_table(env, variant):
    """the development table of t207_d5_agg.py: QREF / FROM are taken from the run's own env (non-negative integers, QREF >= 1)"""
    if not isinstance(env, dict):
        raise QError('env is not a dictionary')
    if variant not in DEV_CTL_EXTRA:
        raise QError('env allow-list: unknown path %r' % (variant,))
    main = dict(DEV_MAIN_ENV)
    for k in DEV_COUNT_KEYS:
        if k not in env or not re.match(r'^[0-9]+$', str(env[k])) or (k.endswith('QREF') and int(env[k]) < 1):
            raise QError('env allow-list: %s missing or not a valid number: %r' % (k, env.get(k)))
        main[k] = env[k]
    k, v = DEV_CTL_EXTRA[variant]
    return env_table({'main': main, 'ctl': {'extra_env': {k: v}, 'omit_env': ['ALF_GPU_TAUDIAG_QVAR']},
                      'ctl_qvar': {'extra_env': {k: v}, 'omit_env': []}})


def env_role(env, variant):
    """development roles of t207_d5_agg.py: 'main' or 'ctl' (rows ctl / ctl_qvar of dev_env_table), else QError"""
    r = role_of_env(env, dev_env_table(env, variant))
    return 'ctl' if r.startswith('ctl') else r


def parse_int_list(s, what):
    if s == '-' or s == '':
        return []
    out = []
    for t in s.split(','):
        if not re.match(r'^[0-9]+$', t):
            raise QError('%s: not a comma list of integers: %r' % (what, s))
        out.append(int(t))
    return out


def derive_run(rundir):
    """Everything the expected content depends on that comes from the run's own records (no processing arguments)."""
    for r in RECORDS:
        read_text(rundir, r)
    prov = read_json(os.path.join(rundir, 'provenance.json'), 'record provenance.json')
    if not isinstance(prov, dict):
        raise QError('provenance.json is not a JSON object')
    if not isint(prov.get('rc')) or prov.get('rc') != 0:
        raise QError('provenance.json: rc = %r (the run did not end normally)' % prov.get('rc'))
    for key, name in (('parameters_sha256', 'parameters'), ('run_log_sha256', 'run.log'), ('info_sha256', 'info')):
        if prov.get(key) != sha256_file(os.path.join(rundir, name)):
            raise QError('provenance.json %s does not match the record %s' % (key, name))
    sm = re.findall(r'^seeds ([0-9a-f]{16}) parameters ([0-9a-f]{16}) *$', read_text(rundir, 'stamp.txt'), re.M)
    if len(sm) != 1:
        raise QError('stamp.txt: seeds/parameters line: expected exactly one, found %d' % len(sm))
    ps = prov.get('seeds_sha256')
    if not (isinstance(ps, str) and re.match(r'^[0-9a-f]{64}$', ps) and ps.startswith(sm[0][0]) and prov['parameters_sha256'].startswith(sm[0][1])):
        raise QError('provenance.json seeds_sha256 / parameters_sha256 do not match the seeds line of stamp.txt')
    bs = prov.get('binary_sha256')
    bm = re.findall(r'^binary ([0-9a-f]{64})(?: .*)?$', read_text(rundir, 'stamp.txt'), re.M)
    if not (isinstance(bs, str) and re.match(r'^[0-9a-f]{64}$', bs)) or len(bm) != 1 or bm[0] != bs:
        raise QError('provenance.json binary_sha256 %r and the binary line of stamp.txt %s do not agree (the ALF binary of the run)' % (bs, bm))
    par = read_text(rundir, 'parameters')
    info = read_text(rundir, 'info')
    log = read_text(rundir, 'run.log')
    stamp = read_text(rundir, 'stamp.txt')
    nbin = int(one(r'^Nbin = ([0-9]+)\b', par, 'parameters Nbin'))
    nsweep = int(one(r'^Nsweep = ([0-9]+)\b', par, 'parameters Nsweep'))
    nwrap = int(one(r'^Nwrap = ([0-9]+)\b', par, 'parameters Nwrap'))
    beta = fnum(one(r'^Beta = ([-+0-9.dDeE]+)', par, 'parameters Beta'))
    dtau = fnum(one(r'^Dtau = ([-+0-9.dDeE]+)', par, 'parameters Dtau'))
    if nbin < 1 or nsweep < 1 or nwrap < 1 or not dtau > 0 or not beta > 0:
        raise QError('parameters: implausible Nbin/Nsweep/Nwrap/Beta/Dtau')
    ltrot = int(round(beta / dtau))
    if abs(beta / dtau - ltrot) > 1e-6 or ltrot < 1:
        raise QError('parameters: Beta/Dtau = %r is not an integer' % (beta / dtau))
    nstm = (ltrot + nwrap - 1) // nwrap
    passes = nbin * nsweep
    # info: the same numbers as the program printed
    i_sw = int(one(r'^ *Sweeps +: +([0-9]+) *$', info, 'info Sweeps'))
    i_bi = int(one(r'^ *Bins +: +([0-9]+) *$', info, 'info Bins'))
    i_nstm = int(one(r'^ *Nstm +: +([0-9]+) *$', info, 'info Nstm'))
    i_nfl = int(one(r'^ *N_FL +: +([0-9]+) *$', info, 'info N_FL'))
    lt = re.findall(r'^ *dtau,Ltrot_eff: +(\S+) +([0-9]+) *$', info, re.M)
    if len(lt) != 1:
        raise QError('info: dtau,Ltrot_eff line: expected exactly one, found %d' % len(lt))
    if i_sw != nsweep or i_bi != nbin:
        raise QError('info Sweeps/Bins %d/%d differ from parameters Nsweep/Nbin %d/%d' % (i_sw, i_bi, nsweep, nbin))
    if int(lt[0][1]) != ltrot or abs(fnum(lt[0][0]) - dtau) > 1e-12:
        raise QError('info dtau,Ltrot_eff %s %s differ from parameters (Dtau %r, Ltrot %d)' % (lt[0][0], lt[0][1], dtau, ltrot))
    if i_nstm != nstm:
        raise QError('info Nstm %d differs from ceil(Ltrot/Nwrap) = %d' % (i_nstm, nstm))
    tp = re.findall(r'^ *GPU tau pass \(207\): ACTIVE .*? N = ([0-9]+), Nstm = ([0-9]+), flavours ([0-9]+), CGR2_2 QR: (.*?) *$', info, re.M)
    if len(tp) != 1:
        raise QError('info: GPU tau pass (207) line: expected exactly one, found %d' % len(tp))
    n, t_nstm, nfl, qrv = int(tp[0][0]), int(tp[0][1]), int(tp[0][2]), tp[0][3]
    if t_nstm != nstm or nfl != i_nfl:
        raise QError('info GPU tau pass line (Nstm %d, flavours %d) differs from info Nstm %d / N_FL %d' % (t_nstm, nfl, nstm, i_nfl))
    if qrv == 'pre-pivoted (default)':
        variant = 'def'
    elif qrv == 'host ZGEQP3 (STABLOG)':
        variant = 'log'
    else:
        raise QError('info: unknown CGR2_2 QR variant %r' % qrv)
    # the env of the run: info line, stamp.txt, provenance.json must agree
    de = re.findall(r'^ *GPU debug env: (ALF_GPU_.*)$', info, re.M)
    if len(de) != 1:
        raise QError('info: GPU debug env line: expected exactly one, found %d' % len(de))
    env_info = parse_env_info(de[0].rstrip())
    se = re.findall(r'^env (.*)$', stamp, re.M)
    if len(se) != 1:
        raise QError('stamp.txt: env line: expected exactly one, found %d' % len(se))
    env_stamp = parse_env_tokens(se[0].split(), 'stamp.txt env')
    if env_info != env_stamp:
        raise QError('env of info %r differs from env of stamp.txt %r' % (env_info, env_stamp))
    if parse_env_tokens(prov.get('env', []), 'provenance env') != env_info:
        raise QError('env of provenance.json differs from the env of info')
    qref = env_info.get('ALF_GPU_TAUDIAG_QREF')
    if qref is None or not re.match(r'^[0-9]+$', qref) or int(qref) < 1:
        raise QError('the run has no ALF_GPU_TAUDIAG_QREF=<stride> (stride >= 1) in its env: not a QREF run')
    stride = int(qref)
    frm = env_info.get('ALF_GPU_TAUDIAG_FROM', '0')
    if not re.match(r'^[0-9]+$', frm):
        raise QError('ALF_GPU_TAUDIAG_FROM=%r' % frm)
    frm = int(frm)
    qhost = env_info.get('ALF_GPU_TAUDIAG_QHOST', '0')
    if qhost not in ('0', '1'):
        raise QError('ALF_GPU_TAUDIAG_QHOST=%r' % qhost)
    qvar = parse_int_list(env_info.get('ALF_GPU_TAUDIAG_QVAR', ''), 'ALF_GPU_TAUDIAG_QVAR')
    tdg = env_info.get('ALF_GPU_TAUDIAG', '0')
    hasdg = tdg in ('1', '2')
    # run.log: the DEBUG line and the TAUQREF line of every dump
    dl = re.findall(r'^ *GPU tau pass: DEBUG ALF_GPU_TAUDIAG_QREF=([0-9]+) active \(dumps to qref/, from pass ([0-9]+), '
                    r'host CGR2_2 outputs ([TF]), device variant masks:((?: [0-9]+)*)\)\s*$', log, re.M)
    if len(dl) != 1:
        raise QError('run.log: GPU tau pass DEBUG line: expected exactly one, found %d' % len(dl))
    if (int(dl[0][0]), int(dl[0][1]), dl[0][2] == 'T', [int(x) for x in dl[0][3].split()]) != (stride, frm, qhost == '1', qvar):
        raise QError('run.log DEBUG line %r differs from the env (stride %d, from %d, qhost %s, qvar %s)' % (dl[0], stride, frm, qhost, qvar))
    # the simulated call sequence of qref_post
    dumps = []
    cnt = 0
    for p in range(max(frm, 1), passes + 1):
        for nst in range(1, nstm + 1):
            for nfe in range(1, nfl + 1):
                cnt += 1
                if (cnt - 1) % stride == 0:
                    dumps.append((len(dumps) + 1, p, nst, nfe))
    if not dumps:
        raise QError('no dump expected (FROM %d, passes %d)' % (frm, passes))
    tq = re.findall(r'^ *TAUQREF .*$', log, re.M)
    if len(tq) != len(dumps):
        raise QError('run.log has %d TAUQREF lines, the run parameters imply %d dumps' % (len(tq), len(dumps)))
    br = []
    for (k, p, nst, nfe), line in zip(dumps, tq):
        m = re.match(r'^ TAUQREF file ([0-9]+) pass ([0-9]+) nst ([0-9]+) nfe ([0-9]+) br ([0-9]+) (qref/q[0-9]{6}\.bin)$', line)
        if not m:
            raise QError('run.log TAUQREF line not parsable: %r' % line)
        v = tuple(int(x) for x in m.groups()[:5])
        if v[:4] != (k, p, nst, nfe) or m.group(6) != 'qref/q%06d.bin' % k:
            raise QError('run.log TAUQREF line %r differs from the simulated call sequence (file %d pass %d nst %d nfe %d)' % (line, k, p, nst, nfe))
        if v[4] not in (1, 2):
            raise QError('run.log TAUQREF line %r: br not 1 or 2' % line)
        br.append(v[4])
    return dict(passes=passes, ltrot=ltrot, nwrap=nwrap, nstm=nstm, nfl=nfl, n=n, variant=variant, stride=stride, frm=frm,
                qhost=qhost == '1', qvar=qvar, taudiag=tdg, hasdg=hasdg, env=env_info, ndumps=len(dumps), br=br,
                seeds_sha256=ps, parameters_sha256=prov['parameters_sha256'], binary_sha256=bs,
                records_sha256={r: sha256_file(os.path.join(rundir, r)) for r in RECORDS})


def keep_rule(nstm, nfl, stride, m):
    """(P, gcd(m, P)): number of distinct call positions visited by the dumps, and the gcd with the thinning step m = cover*keep"""
    c = nstm * nfl
    pp = c // gcd(c, stride)
    return pp, gcd(m, pp)


def make_spec(rundir, path, qs, nul, hostin, cover=1, keep=0):
    if path not in ('def', 'log'):
        raise QError('path must be def or log')
    run = derive_run(rundir)
    if run['variant'] != path:
        raise QError('the run is a %s run (info), the processing says %s' % (run['variant'], path))
    if not (isint(qs) and isint(hostin) and isint(cover) and isint(keep) and isinstance(nul, list) and all(isint(c) and c >= 0 for c in nul)):
        raise QError('bad processing arguments: qs, hostin, cover, keep must be integers and nul a list of integers >= 0 (got qs %r hostin %r '
                     'cover %r keep %r nul %r)' % (qs, hostin, cover, keep, nul))
    if qs < 0 or hostin not in (0, 1) or cover < 1 or keep < 0:
        raise QError('bad processing arguments: qs %r hostin %r cover %r keep %r' % (qs, hostin, cover, keep))
    if len(nul) != len(set(nul)):
        raise QError('duplicate null code')
    if keep > 0:
        # Keep rule (review fix 5).  Call number c (0-based, counted from the first dump pass) is in the block of the pass at position c mod C,
        # C = NSTM*N_FL (blocks x flavours); dump number k is call (k-1)*stride, so the dumps visit P = C/gcd(C, stride) positions cyclically.
        # Keeping every m-th dump (m = cover*keep) visits all P positions equally often iff gcd(m, P) = 1.  Otherwise the kept dumps are a
        # biased subset of the blocks/flavours (e.g. every 40th of stride 2, C = 8: only nst = 1) on which the accuracy bound cannot separate
        # right from wrong.  Refused.  The kept dumps allow pointwise re-checking only; what is recomputable from the evidence is the full
        # qref.txt plus the manifest.
        pp, g = keep_rule(run['nstm'], run['nfl'], run['stride'], cover * keep)
        if g != 1:
            raise QError('keep rule: gcd(cover*keep = %d, P = %d) = %d != 1 (C = NSTM*N_FL = %d, stride %d): the kept dumps would cover only a '
                         'subset of the block/flavour positions' % (cover * keep, pp, g, run['nstm'] * run['nfl'], run['stride']))
    return dict(version=SPEC_VERSION, run=run, args=dict(path=path, qs=qs, nul=nul, hostin=hostin, cover=cover, keep=keep))


def cover_dumps(spec):
    """the expected dumps covered by qref.txt: (number, pass, nst, nfe, br)"""
    run = spec['run']
    out = []
    k = 0
    cnt = 0
    for p in range(max(run['frm'], 1), run['passes'] + 1):
        for nst in range(1, run['nstm'] + 1):
            for nfe in range(1, run['nfl'] + 1):
                cnt += 1
                if (cnt - 1) % run['stride'] == 0:
                    k += 1
                    if (k - 1) % spec['args']['cover'] == 0:
                        out.append((k, p, nst, nfe, run['br'][k - 1]))
    return out


def kept_dumps(spec):
    """the covered dumps that stay after the thinning (all covered ones if keep is 0)"""
    m = spec['args']['cover'] * max(spec['args']['keep'], 1)
    return [d for d in cover_dumps(spec) if (d[0] - 1) % m == 0]


def kept_summary(spec):
    kd = kept_dumps(spec)
    run = spec['run']
    def cnt(i, rng):
        return ' '.join('%d:%d' % (v, sum(1 for d in kd if d[i] == v)) for v in rng)
    return ('QPARSE kept dumps (every %d-th of the run, stride %d): %d of %d; per nst %s; per br %s; per nfe %s'
            % (spec['args']['cover'] * max(spec['args']['keep'], 1), run['stride'], len(kd), run['ndumps'],
               cnt(2, range(1, run['nstm'] + 1)), cnt(4, (1, 2)), cnt(3, range(1, run['nfl'] + 1))))


def expected_names(spec, k, hasb, host='nvh'):
    """the ordered member names of dump number k derived from the run's records and the processing arguments; host: the label the offline
    program prints for the in-run host CGR2_2 of ALF_GPU_TAUDIAG_QHOST=1 (t207_d5_qref.F90 prints 'nvh'; a gate spec declares it)"""
    a, run = spec['args'], spec['run']
    names = ['dev'] + ['var%d' % m for m in run['qvar']]
    if run['qhost']:
        names.append(host)
    names.append('cpu')
    names += ['n%d' % c for c in a['nul']]
    names += ['u1', 'u2']
    if a['qs'] > 0 and (k - 1) % a['qs'] == 0:
        names += ['qself', 'qbr']
    if run['hasdg'] and a['hostin'] == 1:
        names += ['din', 'cpuh', 'cpuhx']
    return names


def load_spec(rundir):
    stored = read_json(os.path.join(rundir, 'qref_spec.json'), 'qref_spec.json (written by the processing before the program runs)')
    a = stored.get('args') if isinstance(stored, dict) else None
    if not isinstance(stored, dict) or stored.get('version') != SPEC_VERSION or not isinstance(a, dict) \
            or set(a) != {'path', 'qs', 'nul', 'hostin', 'cover', 'keep'}:
        raise QError('qref_spec.json: bad version or args')
    fresh = make_spec(rundir, a['path'], a['qs'], a['nul'], a['hostin'], a['cover'], a['keep'])
    if json.loads(json.dumps(fresh)) != stored:
        raise QError('qref_spec.json differs from the spec derived from the records now (records changed after the processing began?)')
    return stored


def fnum_strict(tok, what):
    if not NUMRE.match(tok):
        raise QError('%s: %r is not a number in the format ES12.4E3' % (what, tok))
    v = float(tok)
    if not math.isfinite(v):
        raise QError('%s: %r is not finite' % (what, tok))
    return v


def parse_header(line, what):
    t = line.split(' | ')[0].split()
    keys = ['file', 'pass', 'nst', 'nfe', 'br', 'hbr', 'n', 'ver', 'hasb', 'dg', 'tq']
    if len(t) != 1 + 2 * len(keys) or t[0] != 'QREF' or t[1::2] != keys:
        raise QError('%s: header not parsable: %r' % (what, line.split(' | ')[0]))
    hd = {}
    for k, v in zip(keys[:-1], t[2::2]):
        if not re.match(r'^[0-9]+$', v):
            raise QError('%s: header field %s = %r is not an integer' % (what, k, v))
        hd[k] = int(v)
    if not re.match(r'^[0-9]+\.[0-9]{2}$', t[-1]):
        raise QError('%s: header field tq = %r' % (what, t[-1]))
    hd['tq'] = float(t[-1])
    return hd


def check_ef(e, n, what):
    """Invariant of the program's output, true by construction for every member (errs and qdiff of t207_d5_qref.F90 compute e1..e4 =
    max|d_b| / max|Gq_b| per block and eF = sqrt(sum|d|^2 / sum|Gq|^2) over the 4 n^2 entries): with D = max_b max|d_b| and M = the largest
    |Gq| entry, sum|d|^2 <= 4n^2 D^2, sum|Gq|^2 >= M^2, and D/M <= e of the block holding D, so eF <= 2n max(e1..e4); sum|d|^2 >= max|d_b*|^2
    and sum|Gq|^2 <= 4n^2 M^2 for the block b* holding M, so eF >= e_b*/(2n) >= min(e1..e4)/(2n).  Slack 1e-4 relative for the 5-digit
    printing.  It catches an error number moved between the e slots only when it breaks these bounds (weak: factor 2n either side).
    A SANITY CHECK ONLY: a forged line satisfies it easily.  Authenticity of the numbers comes from re-running the program on the kept dumps
    (t207_d5_gate.py) and, beyond that, from the harness's stamps (207-5)."""
    lo, hi = min(e[0:4]) / (2 * n), 2 * n * max(e[0:4])
    if not (lo * (1 - 1e-4) <= e[4] <= hi * (1 + 1e-4)):
        raise QError('%s: eF = %.4e outside [min(e1..e4)/(2n), 2n max(e1..e4)] = [%.4e, %.4e] (impossible by construction: numbers '
                     'shifted or not from the program; a sanity check only)' % (what, e[4], lo, hi))


def load(qref_txt, allow_v1=False, members=None, host='nvh', raw=None):
    """members (optional, from a gate spec): the exact ordered list of the run's statistical members (every member except the diagnostic
    ones DIAG_NAMES, whose presence follows from qs and hostin); the list derived from the run's records, with the in-run host member
    labelled host, must equal it, else QError.  raw (optional): the bytes of qref_txt as the caller read them: these bytes are parsed, the
    file is not read again.  The Run carries the sha256 of the parsed bytes (qref_sha256) and their 'QREF file' lines (lines)."""
    rundir = os.path.dirname(os.path.abspath(qref_txt))
    spec = load_spec(rundir)
    run = spec['run']
    exp = cover_dumps(spec)
    if members is not None:
        derived = [x for x in expected_names(spec, 1, 1, host) if x not in DIAG_NAMES]
        if derived != members:
            raise QError('member list %s (derived from the run records, host member labelled %r) differs from the member list of the spec %s'
                         % (derived, host, members))
    if raw is None:
        with open(qref_txt, 'rb') as f:
            raw = f.read()
    qsha = hashlib.sha256(raw).hexdigest()
    try:
        txt = raw.decode('ascii')
    except UnicodeDecodeError:
        raise QError('%s: non-ASCII bytes' % qref_txt)
    if not txt.endswith('\n'):
        raise QError('%s: does not end with a newline (truncated?)' % qref_txt)
    lines = txt.split('\n')[:-1]
    if not lines:
        raise QError('%s is empty' % qref_txt)
    end = lines[-1]
    me = re.match(r'^QREFEND dumps ([0-9]+) ok ([0-9]+) err ([0-9]+)$', end)
    if not me:
        raise QError('last line is not QREFEND (program did not finish): %r' % end[:80])
    body = lines[:-1]
    for ln in body:
        if ln.startswith('QREFERR'):
            raise QError('the program reported an error: ' + ln[:400])
    if (int(me.group(1)), int(me.group(2)), int(me.group(3))) != (len(exp), len(exp), 0):
        raise QError('QREFEND %s differs from dumps %d ok %d err 0' % (me.group(0), len(exp), len(exp)))
    if len(body) != len(exp):
        raise QError('%d QREF lines, %d expected' % (len(body), len(exp)))
    points = []
    for (k, p, nst, nfe, br), line in zip(exp, body):
        what = 'dump %d' % k
        if not line.startswith('QREF file '):
            raise QError('%s: unexpected line %r' % (what, line[:80]))
        hd = parse_header(line, what)
        if hd['ver'] != 2 and not (hd['ver'] == 1 and allow_v1):
            raise QError('%s: dump version %d (version 1 dumps only with allow_v1)' % (what, hd['ver']))
        hasb = 1 if (nst < run['nstm'] and hd['ver'] == 2) else 0     # a version 1 dump has no B, Bi
        want = dict(file=k, nfe=nfe, br=br, hbr=br, n=run['n'], hasb=hasb, dg=1 if run['hasdg'] else 0)
        want['pass'] = p
        want['nst'] = nst
        for key, v in want.items():
            if hd[key] != v:
                raise QError('%s: header %s = %d, expected %d' % (what, key, hd[key], v))
        if not hd['tq'] >= 0:
            raise QError('%s: tq' % what)
        names = expected_names(spec, k, hasb, host)
        segs = line.split(' | ')[1:]
        if len(segs) != len(names):
            raise QError('%s: %d members, expected %d (%s)' % (what, len(segs), len(names), ' '.join(names)))
        m = {}
        for name, seg in zip(names, segs):
            t = seg.split()
            if not t or t[0] != name:
                raise QError('%s: member %r where %r was expected (order: %s)' % (what, t[0] if t else seg[:20], name, ' '.join(names)))
            nn = 5 if name in ('qself', 'qbr') else 8
            if len(t) != 1 + nn:
                raise QError('%s: member %s has %d numbers, expected %d' % (what, name, len(t) - 1, nn))
            vals = [fnum_strict(x, '%s member %s' % (what, name)) for x in t[1:]]
            if any(v < 0 for v in vals[:5]):
                raise QError('%s: member %s: negative error number' % (what, name))
            if nn == 8:
                if hasb:
                    if any(v < 0 for v in vals[5:]):
                        raise QError('%s: member %s: negative or missing propagated number although B exists' % (what, name))
                elif t[6:9] != ['-1.0000E+000'] * 3:
                    raise QError('%s: member %s: propagated numbers must be exactly -1 (no B), got %s' % (what, name, ' '.join(t[6:9])))
            elif max(vals) > QSELF_MAX:
                raise QError('%s: quad self-check %s = %.3e > %.0e' % (what, name, max(vals), QSELF_MAX))
            check_ef(vals, run['n'], '%s member %s' % (what, name))
            m[name] = vals
        if 'n0' in m:
            if segs[names.index('n0')].split()[1:] != segs[names.index('cpu')].split()[1:]:
                raise QError('%s: null code 0 differs from the oracle (cpu)' % what)
        points.append(Point(hd, m, names))
    return Run(spec, points, exp, qref_sha256=qsha, lines=list(body))


# ---------------------------------------------------------------- dump directory

def dump_name(k):
    return 'qref/q%06d.bin' % k


def read_manifest(rundir):
    return read_manifest_sha(rundir)[0]


def read_manifest_sha(rundir):
    """(manifest {name: sha256}, sha256 of the manifest bytes that were parsed)"""
    p = os.path.join(rundir, 'qref_manifest.txt')
    if not os.path.isfile(p):
        raise QError('qref_manifest.txt missing')
    d = {}
    with open(p, 'rb') as f:
        rawm = f.read()
    try:
        txtm = rawm.decode('ascii')
    except UnicodeDecodeError:
        raise QError('qref_manifest.txt: non-ASCII bytes')
    for ln in txtm.split('\n')[:-1]:
        m = re.match(r'^([0-9a-f]{64})  (qref/q[0-9]{6}\.bin)$', ln)
        if not m or m.group(2) in d:
            raise QError('qref_manifest.txt: bad or duplicate line %r' % ln)
        d[m.group(2)] = m.group(1)
    return d, hashlib.sha256(rawm).hexdigest()


def check_dumps(rundir, spec, stride, manifest=None):
    """the files in qref/ are exactly the covered dumps with (number-1) mod stride == 0; if a manifest is given their sha256 matches"""
    want = set(dump_name(k) for (k, p, nst, nfe, br) in cover_dumps(spec) if (k - 1) % stride == 0)
    qd = os.path.join(rundir, 'qref')
    if not os.path.isdir(qd):
        raise QError('directory qref/ missing')
    have = set('qref/' + f for f in os.listdir(qd))
    if have != want:
        raise QError('qref/ holds %d files, expected %d; missing %s, unexpected %s' % (
            len(have), len(want), sorted(want - have)[:3], sorted(have - want)[:3]))
    if manifest is not None:
        for f in sorted(want):
            if f not in manifest:
                raise QError('%s not in qref_manifest.txt' % f)
            if sha256_file(os.path.join(rundir, f)) != manifest[f]:
                raise QError('%s: sha256 differs from qref_manifest.txt' % f)
    return len(want)


def verify_dumps(rundir, spec, parsed=None):
    """the dump evidence after the processing: qref_manifest.txt lists exactly the covered dumps, and qref/ holds exactly the kept ones
    (every cover*max(keep,1)-th) with the sha256 of the manifest.  Returns (manifest {name: sha256}, sha256 of the manifest bytes parsed).
    parsed: (manifest, sha256) of read_manifest_sha already called by the caller (the manifest is then not read again)."""
    manifest, msha = parsed if parsed is not None else read_manifest_sha(rundir)
    if set(manifest) != set(dump_name(k) for (k, p, nst, nfe, br) in cover_dumps(spec)):
        raise QError('qref_manifest.txt does not list exactly the covered dumps')
    check_dumps(rundir, spec, spec['args']['cover'] * max(spec['args']['keep'], 1), manifest)
    return manifest, msha


def check_dump_headers(rundir, spec, allow_v1=False):
    """Every kept dump (qref/qNNNNNN.bin) is a dump of this run at its place in the call sequence: the header as t207_d5_qref.F90 reads it
    (13 little-endian int32: magic 20791, version, file, pass, nst, nfe, br, N, path flag, number of QVAR masks, QHOST, TAUDIAG, hasb; version 1
    has 12, no hasb), the QVAR masks that follow it, and the exact file size the header implies (the program's formula).  Every field is compared
    with what the run's own records say (run.log TAUQREF line of the dump, info N / path, env QVAR / QHOST / TAUDIAG, hasb = nst < NSTM).  A text
    file or other data posing as a dump is refused here (bad magic or size).  Returns the number of dumps checked."""
    run = spec['run']
    islog = run['variant'] == 'log'
    n = 0
    for (k, p, nst, nfe, br) in kept_dumps(spec):
        what = dump_name(k)
        f = os.path.join(rundir, what)
        with open(f, 'rb') as fh:
            head = fh.read(52 + 4 * 64)
        size = os.path.getsize(f)
        if len(head) < 48:
            raise QError('%s: %d bytes, shorter than a dump header (not a dump)' % (what, len(head)))
        h = list(struct.unpack('<12i', head[:48]))
        if h[0] != DUMP_MAGIC or h[1] not in (1, 2):
            raise QError('%s: bad dump header (magic %d, version %d; a dump has magic %d, version 2): not a dump' % (what, h[0], h[1], DUMP_MAGIC))
        if h[1] == 1 and not allow_v1:
            raise QError('%s: dump version 1 (no B, Bi; only with allow_v1, never by the gate)' % what)
        if h[1] == 2:
            if len(head) < 52:
                raise QError('%s: short header' % what)
            h.append(struct.unpack('<i', head[48:52])[0])
            hlen = 52
        else:
            h.append(0)
            hlen = 48
        nv = h[9]
        hasb = 1 if (h[1] == 2 and nst < run['nstm']) else 0
        want = [('file', 2, k), ('pass', 3, p), ('nst', 4, nst), ('nfe', 5, nfe), ('br', 6, br), ('N', 7, run['n']), ('path flag', 8, 1 if islog else 0),
                ('number of QVAR masks', 9, len(run['qvar'])), ('QHOST', 10, 1 if run['qhost'] else 0), ('TAUDIAG', 11, 1 if run['hasdg'] else 0),
                ('hasb', 12, hasb)]
        for name, i, v in want:
            if h[i] != v:
                raise QError('%s: header %s = %d, the run records say %d' % (what, name, h[i], v))
        if len(head) < hlen + 4 * nv:
            raise QError('%s: short read (masks)' % what)
        masks = list(struct.unpack('<%di' % nv, head[hlen:hlen + 4 * nv])) if nv else []
        if masks != run['qvar']:
            raise QError('%s: QVAR masks %s in the dump, the run env says %s' % (what, masks, run['qvar']))
        nn = run['n']
        sb = 32 * nn * nn + (8 if islog else 16) * nn
        esz = 4 * (12 + (1 if h[1] == 2 else 0)) + 4 * nv + 2 * sb + (1 + nv + h[10]) * 64 * nn * nn + h[11] * sb + hasb * 32 * nn * nn
        if size != esz:
            raise QError('%s: size %d bytes, the header implies %d (truncated, padded or not a dump)' % (what, size, esz))
        n += 1
    return n


HEX64 = re.compile(r'^[0-9a-f]{64}$')


def proc_record(rundir, spec, prog=None, src=None, qref_sha256=None, manifest_sha256=None):
    """RUNDIR/qref_proc.json, the record t207_d5_proc.sh writes after the program and the check (one implementation for thin, agg, gate):
    rc and check_rc are the integer 0, the arguments are those of qref_spec.json, 'prog' and 'prog_src' are sha256 values (equal to prog /
    src if given), 'qref_sha256' is the sha256 of qref.txt (qref_sha256 given: the hash of the bytes the caller parsed; else of the file as it
    is now), and 'manifest_sha256' (fix round 4) is a sha256, equal to manifest_sha256 if given (the hash of the manifest bytes the caller
    parsed).  Returns the record."""
    pj = os.path.join(rundir, 'qref_proc.json')
    if not os.path.isfile(pj):
        raise QError('qref_proc.json missing in %s (every run needs the record of its processing, t207_d5_proc.sh)' % rundir)
    pr = read_json(pj, 'qref_proc.json')
    if not isinstance(pr, dict):
        raise QError('qref_proc.json is not a JSON object')
    if not (isint(pr.get('rc')) and isint(pr.get('check_rc')) and pr['rc'] == 0 and pr['check_rc'] == 0):
        raise QError('qref_proc.json: rc %r / check_rc %r are not both 0 (the processing did not pass)' % (pr.get('rc'), pr.get('check_rc')))
    a, pa = spec['args'], pr.get('args')
    want = dict(path=a['path'], qs=a['qs'], nul=','.join(str(c) for c in a['nul']), hostin=a['hostin'], cover=a['cover'], keep=a['keep'])
    if not isinstance(pa, dict) or set(pa) != set(want) or any(type(pa[k]) is not type(v) or pa[k] != v for k, v in want.items()):
        raise QError('qref_proc.json args %r differ from the arguments in qref_spec.json %r' % (pa, want))
    for key, what, given in (('prog', 'program', prog), ('prog_src', 'program source', src)):
        v = pr.get(key)
        if not (isinstance(v, str) and HEX64.match(v)):
            raise QError('qref_proc.json: no valid %s sha256 (%s = %r)' % (what, key, v))
        if given is not None and v != given:
            raise QError('qref_proc.json: %s sha256 %s (%s) differs from the required %s' % (what, v, key, given))
    qs = pr.get('qref_sha256')
    if not (isinstance(qs, str) and HEX64.match(qs)):
        raise QError('qref_proc.json: no valid qref_sha256 (%r): the processing did not record the hash of qref.txt' % (qs,))
    if qref_sha256 is None:
        qt = os.path.join(rundir, 'qref.txt')
        if not os.path.isfile(qt):
            raise QError('qref.txt missing in %s' % rundir)
        qref_sha256 = sha256_file(qt)
    if qref_sha256 != qs:
        raise QError('qref.txt: sha256 %s differs from qref_sha256 %s in qref_proc.json (changed after the processing)' % (qref_sha256, qs))
    ms = pr.get('manifest_sha256')
    if not (isinstance(ms, str) and HEX64.match(ms)):
        raise QError('qref_proc.json: no valid manifest_sha256 (%r): the processing did not record the hash of qref_manifest.txt' % (ms,))
    if manifest_sha256 is not None and manifest_sha256 != ms:
        raise QError('qref_manifest.txt: sha256 %s differs from manifest_sha256 %s in qref_proc.json (changed after the processing)'
                     % (manifest_sha256, ms))
    return pr


def thin(rundir, allow_v1=False):
    spec = load_spec(rundir)
    keep = spec['args']['keep']
    if keep < 1:
        raise QError('keep 0: nothing to thin')
    # a passing check is required here, not only in the caller (review fix 6): qref.txt loads strictly, and the processing recorded
    # program exit status 0 and check status 0 for exactly these arguments
    run = load(os.path.join(rundir, 'qref.txt'), allow_v1)
    manifest, msha = read_manifest_sha(rundir)
    try:
        proc_record(rundir, spec, qref_sha256=run.qref_sha256, manifest_sha256=msha)
    except QError as e:
        raise QError('%s: nothing deleted' % e)
    stride = spec['args']['cover'] * keep
    cov = cover_dumps(spec)
    if set(manifest) != set(dump_name(k) for (k, p, nst, nfe, br) in cov):
        raise QError('qref_manifest.txt does not list exactly the covered dumps')
    check_dumps(rundir, spec, spec['args']['cover'], manifest)         # all present, hashes as in the manifest
    qd = os.path.realpath(os.path.join(rundir, 'qref'))
    todel = []
    for (k, p, nst, nfe, br) in cov:
        if (k - 1) % stride == 0:
            continue
        f = os.path.join(rundir, dump_name(k))
        if os.path.islink(f) or os.path.dirname(os.path.realpath(f)) != qd:
            raise QError('%s is a symbolic link or outside qref/: nothing deleted' % f)
        todel.append((dump_name(k), f))
    for name, f in todel:
        os.remove(f)
    with open(os.path.join(rundir, 'qref_thinned.txt'), 'w') as out:
        for name, f in todel:
            out.write(name + '\n')
    return len(todel), len(cov) - len(todel)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('spec')
    s.add_argument('rundir'); s.add_argument('--path', required=True); s.add_argument('--qs', type=int, required=True)
    s.add_argument('--nul', required=True); s.add_argument('--hostin', type=int, required=True)
    s.add_argument('--cover', type=int, default=1); s.add_argument('--keep', type=int, default=0)
    c = sub.add_parser('check')
    c.add_argument('rundir'); c.add_argument('--allow-v1', action='store_true'); c.add_argument('--dumps', action='store_true')
    t = sub.add_parser('thin')
    t.add_argument('rundir'); t.add_argument('--allow-v1', action='store_true')
    a = ap.parse_args()
    try:
        if a.cmd == 'spec':
            out = os.path.join(a.rundir, 'qref_spec.json')
            if os.path.exists(out):
                raise QError('%s exists' % out)
            sp = make_spec(a.rundir, a.path, a.qs, parse_int_list(a.nul, '--nul'), a.hostin, a.cover, a.keep)
            nd = check_dumps(a.rundir, sp, sp['args']['cover'])
            with open(out, 'w') as f:
                json.dump(sp, f, indent=1)
            print('QPARSE spec ok: %d dumps expected in the run, %d covered, %d files in qref/, NSTM %d, flavours %d, stride %d, from pass %d'
                  % (sp['run']['ndumps'], len(cover_dumps(sp)), nd, sp['run']['nstm'], sp['run']['nfl'], sp['run']['stride'], sp['run']['frm']))
        elif a.cmd == 'check':
            r = load(os.path.join(a.rundir, 'qref.txt'), a.allow_v1)
            msg = ''
            if a.dumps:
                sp = r.spec
                nd = check_dumps(a.rundir, sp, sp['args']['cover'] * max(sp['args']['keep'], 1), read_manifest(a.rundir))
                msg = ', %d dumps in qref/ match the manifest' % nd
            print('QPARSE check ok: %d points (%d with B)%s' % (len(r.points), r.n_with_b, msg))
            if r.spec['args']['keep'] > 0:
                print(kept_summary(r.spec))
        else:
            nd, nk = thin(a.rundir, a.allow_v1)
            print('QPARSE thin ok: %d dumps deleted, %d kept' % (nd, nk))
            print(kept_summary(load_spec(a.rundir)))
    except QError as e:
        sys.stderr.write('QPARSE ERROR: %s\n' % e)
        sys.exit(1)


if __name__ == '__main__':
    main()
