#!/usr/bin/env python3
"""Ticket 207-2: verdict of the accuracy test against the quad reference from the aggregate of t207_d5_agg.py.

usage: t207_d5_gate.py --bound B --expect-points N --qs Q [--qs Q ...] [--prog-sha SHA256] [--require-controls] AGG.json   (agg run, measure pTT)
Verdict on the pooled statistic S_ens(pTT) of the runs in AGG.json (one path):
  INPUT_ERROR (exit 2)  the file is not a complete aggregate or not the accepted configuration: wrong measure, a run whose points used
                        differ from the number derived from its parameters or from --expect-points, runs with different member lists,
                        missing or non-finite numbers, pooled counts that are not the sum of the runs, a main run without 'dev'; with
                        --require-controls: no control member.  Per run (inputs[*].args and env, written by t207_d5_agg.py):
                          cover == 1 (every dump processed) and no version 1 dump; qs equal to one of the --qs values and > 0; 0 in the
                          null codes and the null codes contain the reference list; the reference list is exactly the frozen list of
                          the run's path; the env (ALF_GPU_*) is exactly the allowed one for its role (main / ctl, see t207_d5_qparse.py);
                          no run twice (qref.txt sha256, run.log sha256, role+seeds+parameters); with --prog-sha every run's program hash
                          equals it.
                        Any exception (a structurally broken file) is INPUT_ERROR too: the line 'GATE FAIL' is printed only on a real FAIL
  FAIL (exit 1)         the device S_ens (member 'dev' of the main runs) is > B, or a control present in the file (the stand-in
                        'var8' of the main runs, 'ctl_dev' of production control runs) has S_ens <= B
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
    if progsha is not None and need(i, 'prog_sha256', lab) != progsha:
        raise InputError('run %s: program sha256 %r differs from --prog-sha' % (lab, i['prog_sha256']))


def decide(agg, bound, require_controls, expect, qs_ok, progsha=None):
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
        variants.add(need(i, 'variant', lab))
        out.append('run %-8s %-4s points %d used %d expected %d' % (lab, i['role'], npt, nu, ne))
    if len(variants) != 1:
        raise InputError('runs of different paths: %s' % sorted(variants))

    def pooled_s(name):
        e = need(need(pooled, name, 'pooled'), 'S_ens', 'pooled ' + name)
        mean, n = fin(need(e, 'mean', name), name + ' mean'), need(e, 'n', name)
        nexp = sum(runs[l][name]['S_ens']['n'] for l in runs if name in runs[l])
        if n != nexp or n < 1:
            raise InputError('pooled %s: n %r, sum over the runs %r' % (name, n, nexp))
        for l in runs:
            if name in runs[l]:
                fin(runs[l][name]['S_ens']['mean'], '%s %s mean' % (l, name))
        return mean, n

    dev, ndev = pooled_s('dev')
    if ndev != sum(inputs[l]['n_used'] for l in main):
        raise InputError('pooled dev n %d differs from the points of the main runs' % ndev)
    ctl = [n for n in ('var8', 'ctl_dev') if n in pooled]
    if require_controls and not ctl:
        raise InputError('--require-controls: no control member (var8, ctl_dev) in the aggregate')
    verdict = 'PASS'
    reasons = []
    out.append('device S_ens(pTT) pooled over %d points: %+.4f  (bound %+.4f)' % (ndev, dev, bound))
    if not dev <= bound:
        verdict = 'FAIL'; reasons.append('device S_ens %+.4f > bound %+.4f' % (dev, bound))
    for c in ctl:
        v, n = pooled_s(c)
        ok = v > bound
        out.append('control %-7s S_ens(pTT) pooled over %d points: %+.4f  (must be > %+.4f): %s' % (c, n, v, bound, 'killed' if ok else 'NOT KILLED'))
        if not ok:
            verdict = 'FAIL'; reasons.append('control %s S_ens %+.4f <= bound %+.4f' % (c, v, bound))
    return verdict, reasons, out


def main():
    a = sys.argv[1:]
    req = '--require-controls' in a
    a = [x for x in a if x != '--require-controls']
    bound = expect = progsha = None
    qs_ok = set()
    files = []
    try:
        i = 0
        while i < len(a):
            if a[i] in ('--bound', '--expect-points', '--qs', '--prog-sha') and i + 1 < len(a):
                v = a[i + 1]
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
        if bound is None or expect is None or not qs_ok or len(files) != 1:
            raise InputError('usage: t207_d5_gate.py --bound B --expect-points N --qs Q [--qs Q ...] [--prog-sha SHA256] [--require-controls] AGG.json '
                             '(all of --bound, --expect-points, --qs and the file are required; see the header of the script)')
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
