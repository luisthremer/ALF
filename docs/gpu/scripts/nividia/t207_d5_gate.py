#!/usr/bin/env python3
"""Ticket 207-2: verdict of the accuracy test against the quad reference from the aggregate of t207_d5_agg.py.

usage: t207_d5_gate.py --bound B --expect-points N --qs Q [--qs Q ...] --prog-sha SHA256 [--require-controls] AGG.json   (agg run, measure pTT)
--bound, --expect-points and --prog-sha may be given once only; --prog-sha is required (program identity is mandatory).
Verdict on the pooled statistic S_ens(pTT) of the runs in AGG.json (one path):
  INPUT_ERROR (exit 2)  the file is not a complete aggregate or not the accepted configuration: wrong measure, a run whose points used
                        differ from the number derived from its parameters or from --expect-points, runs with different member lists,
                        missing or non-finite numbers, pooled counts that are not the sum of the runs, a main run without 'dev'; with
                        --require-controls: no control member.  Per run (inputs[*].args and env, written by t207_d5_agg.py):
                          cover == 1 (every dump processed) and no version 1 dump; qs equal to one of the --qs values and > 0; 0 in the
                          null codes and the null codes contain the reference list; the reference list is exactly the frozen list of
                          the run's path; the env (ALF_GPU_*) is exactly the allowed one for its role (main / ctl, see t207_d5_qparse.py);
                          no run twice (qref.txt sha256, run.log sha256, role+seeds+parameters); every run's program hash is 64 hex digits
                          and equals --prog-sha; hostin in {0, 1}, keep >= 0 and the keep rule gcd(cover*keep, P) = 1 (P from the run's
                          nstm, nfl, stride, written by the aggregation).
                        Counts and means are tied to the per-run data: the member set of every run is exactly the one expected for its role
                        and env (main: dev, var<QVAR>, cpu, nvh, u1, u2, n<code>; ctl: the same, but only 'ctl_dev' under 'runs'); every
                        runs[l][m].S_ens.n equals n_used of the run; every pooled n is the sum of n_used over the runs that have the member
                        and every pooled mean equals the n-weighted mean of the run means (1e-12 relative).
                        Any exception (a structurally broken file) is INPUT_ERROR too: the line 'GATE FAIL' is printed only on a real FAIL
  FAIL (exit 1)         the device S_ens (member 'dev' of the main runs) is > B, or a control present in the file (table CONTROLS: the
                        stand-in 'var8' of the main runs, 'ctl_dev' of production control runs) has its statistic <= B;
                        with --require-controls every control of the table marked required must be present (else INPUT_ERROR)
  PASS (exit 0)         otherwise
The bound is an argument and is not stored in the aggregate.  Every number of the decision is printed.
"""
import importlib.util, json, math, os, sys
sys.dont_write_bytecode = True      # the loader below imports a sibling script by path: no __pycache__ in the scripts directory
_spec = importlib.util.spec_from_file_location('t207_d5_qparse', os.path.join(os.path.dirname(os.path.abspath(__file__)), 't207_d5_qparse.py'))
QP = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(QP)


class InputError(Exception):
    pass


def need(d, k, what):
    if not isinstance(d, dict) or k not in d:
        raise InputError('%s: key %r missing' % (what, k))
    return d[k]


def fin(x, what):
    if not isinstance(x, (int, float)) or isinstance(x, bool) or not math.isfinite(x):
        raise InputError('%s: %r is not a finite number' % (what, x))
    return x


HEX64 = '0123456789abcdef'
# Control kinds of the gate, one row each.  'member': name in the aggregate's 'pooled'; 'role': role of the runs carrying it; 'stat': an
# entry of STATS (the pooled statistic that must exceed the bound); 'required': demanded by --require-controls.  A later control (graded
# fault ALF_GPU_FAULT_TAU=eps:<x>, a pointwise statistic of emax/eF) is a new row, if needed a new STATS entry and a new member name in
# expected_members() / run_members(); nothing else in decide() is specific to var8 / ctl_dev.
CONTROLS = (
    dict(member='var8', role='main', stat='S_ens', required=True, what='stand-in (variant var8 of the main runs)'),
    dict(member='ctl_dev', role='ctl', stat='S_ens', required=True, what='production control run (TAU_NS=0 / TAUSWAP=8)'),
)
STATS = {'S_ens': lambda entry: entry['S_ens']}   # (a malformed entry raises: caught as INPUT_ERROR by main) aggregate entry (S_cpu, S_ens) -> summary dict with mean, n


def expected_members(i):
    """the members the aggregation writes in inputs[l].members of a run: from its null codes and env"""
    m = {'dev', 'cpu', 'nvh', 'u1', 'u2'} | {'n%d' % c for c in i['args']['nul']}
    if 'ALF_GPU_TAUDIAG_QVAR' in i['env']:
        m.add('var' + i['env']['ALF_GPU_TAUDIAG_QVAR'])
    return sorted(m)


def run_members(i):
    """the keys of runs[l]: every member of a main run, only the renamed device of a control run"""
    return expected_members(i) if i['role'] == 'main' else ['ctl_dev']


def isint(x):
    return isinstance(x, int) and not isinstance(x, bool)


def check_run(lab, i, refs, qs_ok, expect, progsha):
    """the configuration of one run (review fixes 1 to 3, 7): refuse unless it is the accepted one"""
    a = need(i, 'args', 'input ' + lab)
    if need(i, 'v1', lab) is not False:
        raise InputError('run %s: contains version 1 dumps' % lab)
    cover = need(a, 'cover', lab + ' args')
    if not (isinstance(cover, int) and not isinstance(cover, bool) and cover == 1):
        raise InputError('run %s: cover %r: not every dump of the run was processed (qref.txt must cover all of them)' % (lab, cover))
    qs = need(a, 'qs', lab + ' args')
    if not (isinstance(qs, int) and not isinstance(qs, bool) and qs > 0 and qs in qs_ok):
        raise InputError('run %s: qs %r is not > 0 and one of the --qs values %s' % (lab, qs, sorted(qs_ok)))
    nul = need(a, 'nul', lab + ' args')
    if not isinstance(nul, list) or not all(isinstance(c, int) and not isinstance(c, bool) for c in nul) or 0 not in nul:
        raise InputError('run %s: null codes %r (must be a list of integers containing 0)' % (lab, nul))
    path = need(i, 'variant', lab)
    if path not in ('def', 'log'):
        raise InputError('run %s: path %r' % (lab, path))
    if need(a, 'path', lab + ' args') != path:
        raise InputError('run %s: args path %r differs from variant %r' % (lab, a['path'], path))
    if refs != (QP.REF_DEF if path == 'def' else QP.REF_LOG).split(','):
        raise InputError('run %s: the reference list %s is not the frozen list of the %s path' % (lab, refs, path))
    for r in refs:
        if r not in ('cpu', 'nvh') and int(r[1:]) not in nul:      # refs equal the frozen list: names are cpu, nvh, n<code>
            raise InputError('run %s: reference member %s: null code %s is not in the null codes %s' % (lab, r, r[1:], nul))
    if need(i, 'n_used', lab) != expect:
        raise InputError('run %s: %r points used, --expect-points is %d' % (lab, i['n_used'], expect))
    try:
        role = QP.env_role(need(i, 'env', lab), path)
    except QP.QError as e:
        raise InputError('run %s: %s' % (lab, e))
    if role != i['role']:
        raise InputError('run %s: role %s, but its env is that of a %s run' % (lab, i['role'], role))
    ps = need(i, 'prog_sha256', lab)
    if not (isinstance(ps, str) and len(ps) == 64 and all(c in HEX64 for c in ps)):
        raise InputError('run %s: program sha256 %r is not 64 hex digits (no program identity recorded; differs from --prog-sha)' % (lab, ps))
    if ps != progsha:
        raise InputError('run %s: program sha256 %r differs from --prog-sha' % (lab, ps))
    hostin, keep = need(a, 'hostin', lab + ' args'), need(a, 'keep', lab + ' args')
    if not (isint(hostin) and hostin in (0, 1)):
        raise InputError('run %s: hostin %r is not 0 or 1' % (lab, hostin))
    if not (isint(keep) and keep >= 0):
        raise InputError('run %s: keep %r is not an integer >= 0' % (lab, keep))
    if keep > 0:
        geo = [need(i, k, lab) for k in ('nstm', 'nfl', 'stride')]
        if not all(isint(x) and x >= 1 for x in geo):
            raise InputError('run %s: nstm, nfl, stride %r are not positive integers' % (lab, geo))
        pp, g = QP.keep_rule(geo[0], geo[1], geo[2], cover * keep)
        if g != 1:
            raise InputError('run %s: keep rule: gcd(cover*keep = %d, P = %d) = %d != 1: the kept dumps are a biased subset of the block/flavour '
                             'positions' % (lab, cover * keep, pp, g))


def decide(agg, bound, require_controls, expect, qs_ok, progsha):
    out = []
    if need(agg, 'measure', 'aggregate') != 'pTT':
        raise InputError('measure is %r, the gate uses pTT' % agg['measure'])
    refs = need(agg, 'refs', 'aggregate')
    inputs = need(agg, 'inputs', 'aggregate')
    runs = need(agg, 'runs', 'aggregate')
    pooled = need(agg, 'pooled', 'aggregate')
    if not isinstance(refs, list) or 'cpu' not in refs or len(refs) < 2 or len(set(refs)) != len(refs):
        raise InputError('bad reference list %r' % (refs,))
    if not inputs or set(inputs) != set(runs):
        raise InputError('runs and inputs differ')
    main = [l for l, i in inputs.items() if need(i, 'role', 'input ' + l) == 'main']
    if not main:
        raise InputError('no main run')
    members0 = None
    variants = set()
    seen = {}
    for lab, i in inputs.items():
        nu, ne, npt = need(i, 'n_used', lab), need(i, 'n_expected', lab), need(i, 'n_points', lab)
        if not (isinstance(nu, int) and nu == ne and nu > 0 and npt >= nu):
            raise InputError('run %s: %r points used, %r derived from its parameters (of %r points)' % (lab, nu, ne, npt))
        check_run(lab, i, refs, qs_ok, expect, progsha)
        for what, key in (('qref.txt sha256', need(i, 'qref_sha256', lab)), ('run.log sha256', need(i, 'run_log_sha256', lab)),
                          ('run record (role, seeds, parameters)', (i['role'], need(i, 'seeds_sha256', lab), need(i, 'parameters_sha256', lab)))):
            if (what, str(key)) in seen:
                raise InputError('run %s is the same run as %s (same %s)' % (lab, seen[(what, str(key))], what))
            seen[(what, str(key))] = lab
        mem = need(i, 'members', lab)
        for r in refs:
            if r not in mem:
                raise InputError('run %s: reference member %s missing' % (lab, r))
        if i['role'] == 'main':
            if 'dev' not in mem:
                raise InputError('main run %s has no member dev' % lab)
            if members0 is None:
                members0 = mem
            elif mem != members0:
                raise InputError('main runs %s and %s have different member lists' % (main[0], lab))
        elif i['role'] != 'ctl':
            raise InputError('run %s: role %r' % (lab, i['role']))
        if mem != expected_members(i):
            raise InputError('run %s (%s): member list %s differs from the expected %s' % (lab, i['role'], mem, expected_members(i)))
        if sorted(runs[lab]) != run_members(i):
            raise InputError('run %s (%s): the per-run statistics are for %s, expected %s' % (lab, i['role'], sorted(runs[lab]), run_members(i)))
        for m in run_members(i):
            nm = STATS['S_ens'](runs[lab][m])['n']
            if nm - nu != 0 or not isint(nm):
                raise InputError('run %s member %s: S_ens n %r differs from the %d points used of the run' % (lab, m, nm, nu))
        variants.add(need(i, 'variant', lab))
        out.append('run %-8s %-4s points %d used %d expected %d' % (lab, i['role'], npt, nu, ne))
    if len(variants) != 1:
        raise InputError('runs of different paths: %s' % sorted(variants))
    names = set()
    for l in inputs:
        names |= set(run_members(inputs[l]))
    if not isinstance(pooled, dict) or set(pooled) != names:
        raise InputError('pooled members %s differ from the members of the runs %s' % (sorted(pooled) if isinstance(pooled, dict) else pooled, sorted(names)))

    def pooled_s(name, stat='S_ens'):
        """pooled statistic of a member, tied to the per-run data: n = sum of the points used of the runs that have the member, mean = their
        n-weighted mean (1e-12 relative); returns (mean, n)"""
        e = STATS[stat](need(pooled, name, 'pooled'))
        mean, n = fin(need(e, 'mean', name), name + ' mean'), need(e, 'n', name)
        holders = [l for l in runs if name in runs[l]]
        nruns = sum(STATS[stat](runs[l][name])['n'] for l in holders)
        nexp = sum(inputs[l]['n_used'] for l in holders)
        if n != nexp or nruns != nexp or n < 1:
            raise InputError('pooled %s: n %r, points used by the runs %r (sum of the run counts %r)' % (name, n, nexp, nruns))
        wsum = 0.0
        for l in holders:
            s = STATS[stat](runs[l][name])
            wsum += s['n'] * fin(need(s, 'mean', '%s %s' % (l, name)), '%s %s mean' % (l, name))
        wm = wsum / nexp
        if not abs(mean - wm) <= 1e-12 * max(1.0, abs(mean)):
            raise InputError('pooled %s: mean %.17g differs from the n-weighted mean of the runs %.17g' % (name, mean, wm))
        return mean, n

    dev, ndev = pooled_s('dev')
    for name in sorted(pooled):                    # every pooled member is checked against the per-run data, not only the decisive ones
        pooled_s(name)
    ctls = []
    for c in CONTROLS:
        if c['member'] in pooled:
            ctls.append(c)
        elif require_controls and c['required']:
            raise InputError('--require-controls: no control %s (%s) in the aggregate' % (c['member'], c['what']))
    verdict = 'PASS'
    reasons = []
    out.append('device S_ens(pTT) pooled over %d points: %+.4f  (bound %+.4f)' % (ndev, dev, bound))
    if not dev <= bound:
        verdict = 'FAIL'; reasons.append('device S_ens %+.4f > bound %+.4f' % (dev, bound))
    for c in ctls:
        v, n = pooled_s(c['member'], c['stat'])
        ok = v > bound
        out.append('control %-7s %s(pTT) pooled over %d points: %+.4f  (must be > %+.4f): %s' % (c['member'], c['stat'], n, v, bound, 'killed' if ok else 'NOT KILLED'))
        if not ok:
            verdict = 'FAIL'; reasons.append('control %s %s %+.4f <= bound %+.4f' % (c['member'], c['stat'], v, bound))
    return verdict, reasons, out


def main():
    a = sys.argv[1:]
    req = '--require-controls' in a
    a = [x for x in a if x != '--require-controls']
    bound = expect = progsha = None
    qs_ok = set()
    files = []
    given = set()
    try:
        i = 0
        while i < len(a):
            if a[i] in ('--bound', '--expect-points', '--qs', '--prog-sha') and i + 1 < len(a):
                v = a[i + 1]
                if a[i] in given and a[i] != '--qs':
                    raise InputError('%s given twice' % a[i])
                given.add(a[i])
                if a[i] == '--bound':
                    bound = float(v)
                    if not math.isfinite(bound):
                        raise ValueError('--bound is not finite')
                elif a[i] == '--expect-points':
                    expect = int(v)
                    if expect < 1:
                        raise ValueError('--expect-points < 1')
                elif a[i] == '--qs':
                    qs_ok.add(int(v))
                else:
                    progsha = v
                    if len(v) != 64 or any(c not in '0123456789abcdef' for c in v):
                        raise ValueError('--prog-sha is not a sha256')
                i += 2
            else:
                files.append(a[i]); i += 1
        if bound is None or expect is None or not qs_ok or progsha is None or len(files) != 1:
            raise InputError('usage: t207_d5_gate.py --bound B --expect-points N --qs Q [--qs Q ...] --prog-sha SHA256 [--require-controls] AGG.json '
                             '(all of --bound, --expect-points, --qs, --prog-sha and the file are required; see the header of the script)')
        with open(files[0]) as f:
            agg = json.load(f)
        verdict, reasons, lines = decide(agg, bound, req, expect, qs_ok, progsha)
    except InputError as e:
        print('GATE INPUT_ERROR: %s' % e)
        sys.exit(2)
    except Exception as e:                      # every other exception (broken file, wrong types, a bug) is not a verdict (review fix 4)
        print('GATE INPUT_ERROR: %s: %s' % (type(e).__name__, e))
        sys.exit(2)
    for l in lines:
        print('  ' + l)
    print('GATE %s%s' % (verdict, (': ' + '; '.join(reasons)) if reasons else ''))
    sys.exit(0 if verdict == 'PASS' else 1)


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except BaseException as e:                  # nothing escapes as a traceback with another exit status
        print('GATE INPUT_ERROR: %s: %s' % (type(e).__name__, e))
        sys.exit(2)
