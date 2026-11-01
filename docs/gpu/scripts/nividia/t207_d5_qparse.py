#!/usr/bin/env python3
"""Ticket 207-2: strict reader of the output of the offline quad-reference program (t207_d5_qref.F90) and of the records of
the run that produced the dumps.  Nothing is skipped silently: every deviation from the expected content is a QError.

Library (callers load it by file path and run python3 -I):
    load(qref_txt, allow_v1=False) -> Run      the validated points of RUNDIR/qref.txt (RUNDIR = its directory)
CLI (python3 -I t207_d5_qparse.py ...):
    spec  RUNDIR --path def|log --qs Q --nul CODES|- --hostin 0|1 [--cover C] [--keep K]
                          derive the expected content from the run's own records, check the dump directory (before the
                          processing), write RUNDIR/qref_spec.json
    check RUNDIR [--allow-v1] [--dumps]      validate RUNDIR/qref.txt against the spec re-derived from the records
                          (--dumps: and the dump directory after the processing)
    thin  RUNDIR          delete the dumps that are not every (cover*keep)-th, only those listed in RUNDIR/qref_manifest.txt
                          whose sha256 still matches and that are not symbolic links; nothing is deleted if one check fails

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
import argparse, hashlib, json, math, os, re, sys

QSELF_MAX = 1e-25          # input check of the proposed criterion (207 dev5): quad vs itself (row-permuted QR input, other branch)
NUMRE = re.compile(r'^-?[0-9]\.[0-9]{4}E[+-][0-9]{3}$')
REF_DEF = 'cpu,nvh,n4,n100004,n1,n100001,n2,n1608,n4808,n5,n100005'
REF_LOG = 'cpu,nvh,n4,n100004,n1,n100001,n2,n5,n100005'
DIAG_NAMES = ('qself', 'qbr', 'din', 'cpuh', 'cpuhx')     # members that are not error statistics of a G
SPEC_VERSION = 1
RECORDS = ('parameters', 'info', 'run.log', 'stamp.txt', 'provenance.json')


class QError(Exception):
    pass


class Point(object):
    __slots__ = ('hd', 'm', 'names', 'hasb')

    def __init__(self, hd, m, names):
        self.hd, self.m, self.names, self.hasb = hd, m, names, hd['hasb'] == 1


class Run(object):
    def __init__(self, spec, points, dumps):
        self.spec, self.points, self.dumps = spec, points, dumps

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
    pp = os.path.join(rundir, 'provenance.json')
    if not os.path.isfile(pp):
        raise QError('record provenance.json missing in %s' % rundir)
    with open(pp) as f:
        prov = json.load(f)
    if prov.get('rc') != 0:
        raise QError('provenance.json: rc = %r (the run did not end normally)' % prov.get('rc'))
    for key, name in (('parameters_sha256', 'parameters'), ('run_log_sha256', 'run.log'), ('info_sha256', 'info')):
        if prov.get(key) != sha256_file(os.path.join(rundir, name)):
            raise QError('provenance.json %s does not match the record %s' % (key, name))
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
                records_sha256={r: sha256_file(os.path.join(rundir, r)) for r in RECORDS})


def make_spec(rundir, path, qs, nul, hostin, cover=1, keep=0):
    if path not in ('def', 'log'):
        raise QError('path must be def or log')
    run = derive_run(rundir)
    if run['variant'] != path:
        raise QError('the run is a %s run (info), the processing says %s' % (run['variant'], path))
    if qs < 0 or hostin not in (0, 1) or cover < 1 or keep < 0:
        raise QError('bad processing arguments')
    if len(nul) != len(set(nul)):
        raise QError('duplicate null code')
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


def expected_names(spec, k, hasb):
    a, run = spec['args'], spec['run']
    names = ['dev'] + ['var%d' % m for m in run['qvar']]
    if run['qhost']:
        names.append('nvh')
    names.append('cpu')
    names += ['n%d' % c for c in a['nul']]
    names += ['u1', 'u2']
    if a['qs'] > 0 and (k - 1) % a['qs'] == 0:
        names += ['qself', 'qbr']
    if run['hasdg'] and a['hostin'] == 1:
        names += ['din', 'cpuh', 'cpuhx']
    return names


def load_spec(rundir):
    p = os.path.join(rundir, 'qref_spec.json')
    if not os.path.isfile(p):
        raise QError('%s missing (the processing writes it before the program runs)' % p)
    with open(p) as f:
        stored = json.load(f)
    a = stored.get('args')
    if stored.get('version') != SPEC_VERSION or not isinstance(a, dict):
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


def load(qref_txt, allow_v1=False):
    rundir = os.path.dirname(os.path.abspath(qref_txt))
    spec = load_spec(rundir)
    run = spec['run']
    exp = cover_dumps(spec)
    with open(qref_txt, 'rb') as f:
        raw = f.read()
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
            raise QError(ln[:400])
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
        names = expected_names(spec, k, hasb)
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
            m[name] = vals
        if 'n0' in m:
            if segs[names.index('n0')].split()[1:] != segs[names.index('cpu')].split()[1:]:
                raise QError('%s: null code 0 differs from the oracle (cpu)' % what)
        points.append(Point(hd, m, names))
    return Run(spec, points, exp)


# ---------------------------------------------------------------- dump directory

def dump_name(k):
    return 'qref/q%06d.bin' % k


def read_manifest(rundir):
    p = os.path.join(rundir, 'qref_manifest.txt')
    if not os.path.isfile(p):
        raise QError('qref_manifest.txt missing')
    d = {}
    with open(p) as f:
        for ln in f.read().split('\n')[:-1]:
            m = re.match(r'^([0-9a-f]{64})  (qref/q[0-9]{6}\.bin)$', ln)
            if not m or m.group(2) in d:
                raise QError('qref_manifest.txt: bad or duplicate line %r' % ln)
            d[m.group(2)] = m.group(1)
    return d


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


def thin(rundir):
    spec = load_spec(rundir)
    keep = spec['args']['keep']
    if keep < 1:
        raise QError('keep 0: nothing to thin')
    stride = spec['args']['cover'] * keep
    manifest = read_manifest(rundir)
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
    t.add_argument('rundir')
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
        else:
            nd, nk = thin(a.rundir)
            print('QPARSE thin ok: %d dumps deleted, %d kept' % (nd, nk))
    except QError as e:
        sys.stderr.write('QPARSE ERROR: %s\n' % e)
        sys.exit(1)


if __name__ == '__main__':
    main()
