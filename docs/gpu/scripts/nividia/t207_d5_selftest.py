#!/usr/bin/env python3
"""Ticket 207-2: mutation self-tests of the QREF analysis pipeline (nividia, CPU only).   usage: t207_d5_selftest.py WORKDIR [text|dump|all]

Needs WORKDIR/runs/{m_def_dev1..4,d_def_dev1,m_log_dev1}/hub12x12_cb_p produced by t207_d5_repro.sh (the new pipeline on the kept dev5
dumps) and WORKDIR/bin/qref_pp.out, qref_log.out.  Everything it writes goes to WORKDIR/st/ (created; refuses to start if it exists).
  baseline   the unmutated inputs must be accepted (agg exit 0) and the gate must say PASS (device S_ens(pTT) pooled over the four
             default-path seed sets <= bound, both controls > bound)
  text       mutations of qref.txt, of the run's records and of the spec: each must end in the error named (the strict parser, agg,
             stat) or, for the swap of device and var8, in GATE FAIL
  gate       mutations of the aggregate JSON and of the bound: INPUT_ERROR / FAIL as named
  dump       mutations of real dumps (NaN, Infinity, zeroed state, truncation, extra bytes, bad header) run through t207_d5_proc.sh
             (program + strict check): the program must print a QREFERR line for the dump, exit with status 3, the check must fail and no
             dump may be deleted; plus the directory checks (extra / missing dump, stale manifest, symlink and sha256 guards of the thinning)
Output: one table line per test, generated from the actual outcomes; exit status 1 if any test did not give the expected outcome.
The bound used is 0.15 (the proposed criterion of 207 dev5; the self-test is not a calibration and does not change any bound).
"""
import json, os, re, shutil, struct, subprocess, sys

W = os.path.abspath(sys.argv[1])
WHAT = sys.argv[2] if len(sys.argv) > 2 else 'all'
HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
BOUND = 0.15
ST = os.path.join(W, 'st')
RESULTS = []


def sh(argv, env=None, timeout=3600):
    p = subprocess.run(argv, capture_output=True, text=True, env=env, timeout=timeout)
    return p.returncode, (p.stdout + p.stderr)


def py(script, *args):
    return sh([PY, '-I', os.path.join(HERE, script)] + list(args))


def record(name, expect, ok, detail):
    RESULTS.append((name, expect, ok, detail))
    print('%-34s %-34s %-5s %s' % (name, expect, 'ok' if ok else 'WRONG', detail.replace('\n', ' ')[:150]), flush=True)


def base(tag):
    return os.path.join(W, 'runs', tag, 'hub12x12_cb_p')


REC = ('parameters', 'info', 'run.log', 'stamp.txt', 'provenance.json', 'qref_spec.json')


def mkrun(name, src_tag, qtxt=None, rec=None, spec=None, fixprov=False):
    """a run directory st/name/hub12x12_cb_p: the records of src_tag, qref.txt = qtxt (default: the original), optional edits of records"""
    d = os.path.join(ST, name, 'hub12x12_cb_p')
    os.makedirs(d)
    for f in REC:
        shutil.copy2(os.path.join(base(src_tag), f), os.path.join(d, f))
    if qtxt is None:
        shutil.copy2(os.path.join(base(src_tag), 'qref.txt'), os.path.join(d, 'qref.txt'))
    else:
        with open(os.path.join(d, 'qref.txt'), 'wb') as f:
            f.write(qtxt if isinstance(qtxt, bytes) else qtxt.encode('latin-1'))
    if rec:
        for f, fn in rec.items():
            p = os.path.join(d, f)
            t = open(p).read()
            t2 = fn(t)
            assert t2 != t, 'record mutation %s had no effect' % f
            open(p, 'w').write(t2)
    if fixprov:                       # keep provenance.json consistent with the edited records, so that the inner checks are reached
        pj = os.path.join(d, 'provenance.json')
        pv = json.load(open(pj))
        import hashlib
        for key, f in (('parameters_sha256', 'parameters'), ('run_log_sha256', 'run.log'), ('info_sha256', 'info')):
            pv[key] = hashlib.sha256(open(os.path.join(d, f), 'rb').read()).hexdigest()
        json.dump(pv, open(pj, 'w'), indent=1)
    if spec:
        p = os.path.join(d, 'qref_spec.json')
        s = json.load(open(p))
        spec(s)
        json.dump(s, open(p, 'w'), indent=1)
    return os.path.join(d, 'qref.txt')


def orig(tag):
    return open(os.path.join(base(tag), 'qref.txt')).read()


# ---- text editing helpers
def edit_line(txt, i, fn):
    ls = txt.split('\n')
    ls[i] = fn(ls[i])
    return '\n'.join(ls)


def seg_edit(name, fn):
    def f(line):
        segs = line.split(' | ')
        for j, s in enumerate(segs):
            if j and s.split()[0] == name:
                t = s.split()
                segs[j] = ' '.join(fn(t))
                return ' | '.join(segs)
        raise AssertionError('member %s not in the line' % name)
    return f


def set_tok(name, i, val):
    def fn(t):
        t = list(t)
        t[1 + i] = val
        return t
    return seg_edit(name, fn)


A_SET = ['dev1', 'dev2', 'dev3', 'dev4']


def agg_args(out, mut_path=None, mut_label='dev1', refs='def', measure='pTT', extra=None, sets=None):
    a = [out, measure, '--ref', refs] + (extra or [])
    for s in (sets or A_SET):
        p = mut_path if (s == mut_label and mut_path) else os.path.join(base('m_def_' + s), 'qref.txt')
        a.append('main:%s=%s' % (s, p))
    return a


def expect_error(name, qpath, substr, label='dev1', **kw):
    out = os.path.join(ST, 'agg_%s.json' % name)
    rc, msg = py('t207_d5_agg.py', *agg_args(out, qpath, label, **kw))
    ok = rc == 1 and substr in msg and not os.path.exists(out)
    record(name, 'error: %s' % substr[:24], ok, 'rc %d: %s' % (rc, msg.strip().split('\n')[-1][:110]))


def text_tests():
    t0 = orig('m_def_dev1')
    nl = t0.split('\n')
    assert len(nl) == 14 + 2 and nl[-1] == ''
    L = lambda f: f                                                              # noqa: E731
    tests = [
        ('nan_dev_pTT', lambda t: edit_line(t, 1, set_tok('dev', 5, 'NaN')), 'not a number in the format ES12.4E3'),
        ('inf_dev_pTT', lambda t: edit_line(t, 1, set_tok('dev', 5, 'Infinity')), 'not a number in the format'),
        ('minus_inf', lambda t: edit_line(t, 4, set_tok('cpu', 0, '-Infinity')), 'not a number in the format'),
        ('nan_in_n0', lambda t: edit_line(t, 2, set_tok('n0', 3, 'NaN')), 'not a number in the format'),
        ('eless_number', lambda t: edit_line(t, 1, set_tok('var8', 2, '1.798+308')), 'not a number in the format'),
        ('old_ES10.3_format', lambda t: edit_line(t, 1, set_tok('nvh', 1, '1.234E-13')), 'not a number in the format'),
        ('missing_member_n4', lambda t: edit_line(t, 2, lambda l: ' | '.join(s for s in l.split(' | ') if not s.startswith('n4 '))), 'members, expected'),
        ('missing_member_dev', lambda t: edit_line(t, 2, lambda l: ' | '.join(s for s in l.split(' | ') if not s.startswith('dev '))), 'members, expected'),
        ('missing_point', lambda t: '\n'.join(nl[:5] + nl[6:]), '13 QREF lines, 14 expected'),
        ('missing_point_and_end_fixed', lambda t: '\n'.join(nl[:5] + nl[6:-2] + ['QREFEND dumps 13 ok 13 err 0', '']), 'differs from dumps 14'),
        ('missing_QREFEND', lambda t: '\n'.join(nl[:-2] + ['']), 'last line is not QREFEND'),
        ('truncated_mid_line', lambda t: t[:-60], 'does not end with a newline'),
        ('field_moved_dev_to_var8', lambda t: edit_line(edit_line(t, 3, seg_edit('dev', lambda x: x[:-1])), 3,
                                                       seg_edit('var8', lambda x: [x[0], '1.0000E-010'] + x[1:])), 'member dev has 7 numbers'),
        ('unknown_extra_member', lambda t: edit_line(t, 3, lambda l: l + ' | zz ' + ' '.join(['1.0000E-010'] * 8)), 'members, expected'),
        ('renamed_member', lambda t: edit_line(t, 3, lambda l: l.replace(' | cpu ', ' | cpux ', 1)), "member 'cpux' where 'cpu' was expected"),
        ('member_order_swapped', lambda t: edit_line(t, 3, lambda l: _swap_members(l, 'n4', 'n100004')), 'was expected (order'),
        ('duplicate_point', lambda t: '\n'.join(nl[:3] + [nl[1]] + nl[4:]), 'header file'),
        ('points_swapped', lambda t: '\n'.join([nl[0], nl[2], nl[1]] + nl[3:]), 'header file'),
        ('header_pass_changed', lambda t: edit_line(t, 2, lambda l: l.replace(' pass ', ' pass 9', 1)), 'header pass'),
        ('header_hasb_flipped', lambda t: edit_line(t, 2, lambda l: l.replace(' hasb 1 ', ' hasb 0 ', 1)), 'header hasb'),
        ('header_br_and_hbr_flipped', lambda t: edit_line(t, 2, lambda l: re.sub(r' br ([12]) hbr ([12])', lambda m: ' br %d hbr %d' % (3 - int(m.group(1)), 3 - int(m.group(2))), l, 1)), 'header br'),
        ('header_br_flipped', lambda t: edit_line(t, 2, lambda l: re.sub(r' br ([12]) hbr ([12])', lambda m: ' br %s hbr %s' % (m.group(1), 3 - int(m.group(2))), l, 1)), 'header hbr'),
        ('header_n_changed', lambda t: edit_line(t, 2, lambda l: l.replace(' n 144 ', ' n 145 ', 1)), 'header n'),
        ('n0_differs_from_cpu', lambda t: edit_line(t, 2, set_tok('n0', 7, '9.9999E-010')), 'null code 0 differs'),
        ('QREFERR_line', lambda t: edit_line(t, 2, lambda l: 'QREFERR dump qref/q000041.bin : non-finite number in B'), 'QREFERR dump'),
        ('QREFEND_err_counts', lambda t: t.replace('QREFEND dumps 14 ok 14 err 0', 'QREFEND dumps 14 ok 13 err 1'), 'QREFEND'),
        ('negative_error', lambda t: edit_line(t, 1, set_tok('dev', 0, '-1.0000E-013')), 'negative error number'),
        ('minus1_on_point_with_B', lambda t: edit_line(t, 1, set_tok('dev', 5, '-1.0000E+000')), 'negative or missing propagated number'),
        ('qself_too_large', lambda t: edit_line(t, 1, set_tok('qself', 2, '1.0000E-020')), 'quad self-check'),
        ('blank_line', lambda t: '\n'.join(nl[:3] + [''] + nl[3:]), '15 QREF lines, 14 expected'),
        ('trailing_garbage_line', lambda t: '\n'.join(nl[:-2] + ['QREF junk', nl[-2], '']), '15 QREF lines, 14 expected'),
        ('non_ascii_byte', lambda t: t.replace('QREF file 1 ', 'QREF file 1\xb7 ', 1), 'non-ASCII'),
        ('duplicate_member_cpu', lambda t: edit_line(t, 3, lambda l: _dup_member(l, 'cpu')), 'members, expected'),
    ]
    for name, mut, substr in tests:
        q = mkrun(name, 'm_def_dev1', mut(t0))
        expect_error(name, q, substr)
    # zero error (well-formed): the parser accepts it, agg refuses to take a log
    q = mkrun('zero_error_dev_pTT', 'm_def_dev1', edit_line(t0, 1, set_tok('dev', 5, '0.0000E+000')))
    expect_error('zero_error_dev_pTT', q, 'not positive and finite')
    # well-formed huge value (the old clamp value): accepted as data (a very large error, no longer a placeholder): the pooled S_ens moves
    q = mkrun('overflow_value_dev_pTT', 'm_def_dev1', edit_line(t0, 1, set_tok('dev', 5, '1.7977E+308')))
    expect_error('overflow_value_dev_pTT', q, 'is not finite')
    q = mkrun('huge_value_dev_pTT', 'm_def_dev1', edit_line(t0, 1, set_tok('dev', 5, '1.0000E+100')))
    out = os.path.join(ST, 'agg_huge.json')
    rc, msg = py('t207_d5_agg.py', *agg_args(out, q))
    rc2, g = py('t207_d5_gate.py', '--bound', str(BOUND), out) if rc == 0 else (-1, msg)
    record('huge_value_dev_pTT', 'GATE FAIL', rc2 == 1 and 'GATE FAIL' in g, g.strip().split('\n')[-1][:110])
    # dev and var8 swapped at every point: the gate catches it
    def swap(t):
        ls = t.split('\n')
        for i, l in enumerate(ls):
            if l.startswith('QREF file '):
                segs = l.split(' | ')
                d = [j for j, s in enumerate(segs) if j and s.split()[0] == 'dev'][0]
                v = [j for j, s in enumerate(segs) if j and s.split()[0] == 'var8'][0]
                nd, nv = segs[d].split()[1:], segs[v].split()[1:]
                segs[d] = 'dev ' + ' '.join(nv)
                segs[v] = 'var8 ' + ' '.join(nd)
                ls[i] = ' | '.join(segs)
        return '\n'.join(ls)
    sets = {s: os.path.join(base('m_def_' + s), 'qref.txt') for s in A_SET}
    for lab in A_SET:
        sets[lab] = mkrun('swap_dev_var8_' + lab, 'm_def_' + lab, swap(orig('m_def_' + lab)))
    out = os.path.join(ST, 'agg_swap.json')
    rc, msg = py('t207_d5_agg.py', out, 'pTT', '--ref', 'def', *['main:%s=%s' % (s, sets[s]) for s in A_SET])
    rc2, g = py('t207_d5_gate.py', '--bound', str(BOUND), out) if rc == 0 else (-1, msg)
    record('dev_var8_swapped_all_sets', 'GATE FAIL', rc == 0 and rc2 == 1 and 'GATE FAIL' in g, ' '.join(l.strip() for l in g.strip().split('\n') if 'S_ens' in l or 'GATE' in l)[:150])
    # records and spec
    rtests = [
        ('records_run.log_edited_after_run', dict(**{'run.log': lambda t: t.replace(' TAUQREF file 2 ', ' TAUQREFX file 2 ', 1)}), None, False, 'does not match the record run.log'),
        ('records_info_edited_after_run', dict(info=lambda t: re.sub(r'(Nstm +: +)4', r'\g<1>5', t, 1)), None, False, 'does not match the record info'),
        ('records_parameters_edited_after_run', dict(parameters=lambda t: t.replace('Nsweep = 20 ', 'Nsweep = 21 ', 1)), None, False, 'does not match the record parameters'),
        ('parameters_Nsweep_changed', dict(parameters=lambda t: t.replace('Nsweep = 20 ', 'Nsweep = 21 ', 1)), None, True, 'Sweeps/Bins'),
        ('parameters_Nwrap_changed', dict(parameters=lambda t: t.replace('Nwrap = 10 ', 'Nwrap = 8 ', 1)), None, True, 'Nstm'),
        ('runlog_TAUQREF_line_dropped', dict(**{'run.log': lambda t: t.replace(' TAUQREF file 2 ', ' TAUQREFX file 2 ', 1)}), None, True, 'TAUQREF lines'),
        ('runlog_TAUQREF_br_changed', dict(**{'run.log': lambda t: re.sub(r'( TAUQREF file 21 pass \d+ nst \d+ nfe \d+ br )(\d)', lambda m: m.group(1) + str(3 - int(m.group(2))), t, 1)}), None, True, 'qref_spec.json differs'),
        ('runlog_TAUQREF_pass_changed', dict(**{'run.log': lambda t: t.replace(' TAUQREF file 21 pass ', ' TAUQREF file 21 pass 9', 1)}), None, True, 'simulated call sequence'),
        ('stamp_env_var_removed', dict(**{'stamp.txt': lambda t: t.replace(' ALF_GPU_TAUDIAG_QHOST=1', '', 1)}), None, False, 'differs from env of stamp.txt'),
        ('stamp_env_FROM_changed', dict(**{'stamp.txt': lambda t: t.replace('ALF_GPU_TAUDIAG_FROM=331', 'ALF_GPU_TAUDIAG_FROM=333', 1)}), None, False, 'differs from env of stamp.txt'),
        ('provenance_rc_1', dict(**{'provenance.json': lambda t: t.replace('"rc": 0', '"rc": 1', 1)}), None, False, 'rc = 1'),
        ('info_Nstm_changed', dict(info=lambda t: re.sub(r'(Nstm +: +)4', r'\g<1>5', t, 1)), None, True, 'info Nstm 5'),
        ('spec_null_code_dropped', None, lambda s: s['args']['nul'].remove(4), False, 'members, expected'),
        ('spec_qs_changed', None, lambda s: s['args'].__setitem__('qs', 40), False, 'members, expected'),
        ('spec_cover_changed', None, lambda s: s['args'].__setitem__('cover', 40), False, 'differs from dumps 7'),
        ('spec_run_field_changed', None, lambda s: s['run'].__setitem__('nfl', 3), False, 'qref_spec.json differs'),
        ('spec_missing', None, 'rm', False, 'qref_spec.json'),
    ]
    for name, rec, spec, fixprov, substr in rtests:
        if spec == 'rm':
            q = mkrun(name, 'm_def_dev1')
            os.remove(os.path.join(os.path.dirname(q), 'qref_spec.json'))
        else:
            q = mkrun(name, 'm_def_dev1', rec=rec, spec=spec, fixprov=fixprov)
        expect_error(name, q, substr)
    # dump version 1 (smoke run s_log): no B, Bi; every propagated number must be exactly -1; refused without --allow-v1
    sl = orig('s_log')
    sll = sl.split('\n')
    q = os.path.join(base('s_log'), 'qref.txt')
    rc, msg = py('t207_d5_agg.py', os.path.join(ST, 'x0.json'), 'emax', '--ref', 'log', 'main:s=' + q)
    record('v1_dumps_without_allow_v1', 'error: only with allow_v1', rc == 1 and 'only with allow_v1' in msg, msg.strip().split('\n')[-1][:110])
    rc, msg = py('t207_d5_agg.py', os.path.join(ST, 'x0.json'), 'emax', '--ref', 'log', '--allow-v1', 'main:s=' + q)
    record('v1_dumps_with_allow_v1_emax', 'ok', rc == 0, 'agg accepted the v1 run' if rc == 0 else msg.strip().split('\n')[-1][:110])
    rc, msg = py('t207_d5_agg.py', os.path.join(ST, 'x0b.json'), 'pTT', '--ref', 'log', '--allow-v1', 'main:s=' + q)
    record('v1_dumps_pTT_has_no_points', 'error: 0 points', rc == 1 and '0 points used' in msg, msg.strip().split('\n')[-1][:110])
    q = mkrun('nob_propagated_not_minus1', 's_log', edit_line(sl, 2, set_tok('dev', 5, '1.0000E-010')))
    rc, msg = py('t207_d5_agg.py', os.path.join(ST, 'x1.json'), 'emax', '--ref', 'log', '--allow-v1', 'main:s=' + q)
    record('nob_propagated_not_minus1', 'error: exactly -1', rc == 1 and 'exactly -1' in msg, msg.strip().split('\n')[-1][:110])
    q = mkrun('nob_header_hasb_1', 's_log', edit_line(sl, 2, lambda l: l.replace(' hasb 0 ', ' hasb 1 ', 1)))
    rc, msg = py('t207_d5_agg.py', os.path.join(ST, 'x1b.json'), 'emax', '--ref', 'log', '--allow-v1', 'main:s=' + q)
    record('v1_header_hasb_1', 'error: header hasb', rc == 1 and 'header hasb' in msg, msg.strip().split('\n')[-1][:110])
    dt = orig('d_def_dev1')
    q = mkrun('d_base', 'd_def_dev1')
    rc, msg = py('t207_d5_agg.py', os.path.join(ST, 'agg_d_base.json'), 'pTT', '--ref', 'def', 'main:d=' + q)
    record('d_def_baseline_agg', 'ok', rc == 0, 'agg accepted the run with din/cpuh/cpuhx members' if rc == 0 else msg.strip().split('\n')[-1][:100])
    q = mkrun('nob_cpuh_missing', 'd_def_dev1', edit_line(dt, 1, lambda l: ' | '.join(s for s in l.split(' | ') if not s.startswith('cpuh '))))
    rc, msg = py('t207_d5_agg.py', os.path.join(ST, 'x2.json'), 'pTT', '--ref', 'def', 'main:d=' + q)
    record('din_member_cpuh_missing', 'error: members, expected', rc == 1 and 'members, expected' in msg, msg.strip().split('\n')[-1][:110])
    q = mkrun('din_nan', 'd_def_dev1', edit_line(dt, 1, set_tok('din', 4, 'NaN')))
    rc, msg = py('t207_d5_stat.py', os.path.join(ST, 'x3.json'), 'din=' + q)
    record('din_NaN_via_stat', 'error: not a number', rc == 1 and 'not a number' in msg, msg.strip().split('\n')[-1][:110])
    q = mkrun('stat_nan_dev', 'm_def_dev1', edit_line(t0, 1, set_tok('dev', 5, 'NaN')))
    rc, msg = py('t207_d5_stat.py', os.path.join(ST, 'x4.json'), 'dev1=' + q)
    record('NaN_via_stat', 'error: not a number', rc == 1 and 'not a number' in msg, msg.strip().split('\n')[-1][:110])
    rc, msg = py('t207_d5_stat.py', os.path.join(ST, 'x5.json'), 'dev1=' + os.path.join(base('m_def_dev1'), 'qref.txt'))
    record('stat_baseline', 'ok', rc == 0, 'stat accepted the unmutated run' if rc == 0 else msg[-110:])
    # old-format record (dev5) is refused (the pipeline reads only its own output)
    old = os.path.join(os.environ.get('OLDD5', os.path.expanduser('~/alf-gpu/work/t207/dev5')), 'runs/m_def_dev1/hub12x12_cb_p/qref.txt')
    d = os.path.join(ST, 'old_format', 'hub12x12_cb_p')
    os.makedirs(d)
    for f in REC:
        shutil.copy2(os.path.join(base('m_def_dev1'), f), os.path.join(d, f))
    shutil.copy2(old, os.path.join(d, 'qref.txt'))
    rc, msg = py('t207_d5_stat.py', os.path.join(ST, 'x6.json'), 'old=' + os.path.join(d, 'qref.txt'))
    record('dev5_old_format_refused', 'error: last line is not QREFEND', rc == 1 and 'last line is not QREFEND' in msg, msg.strip().split('\n')[-1][:110])
    # reference list tests
    q = os.path.join(base('m_def_dev1'), 'qref.txt')
    for name, refs, extra, substr, which in [
        ('ref_unknown_name', 'cpu,foo', None, 'not a reference name', None),
        ('ref_dev_listed', 'cpu,dev', None, 'not a reference name', None),
        ('ref_n0_listed', 'cpu,n0', None, 'not a reference name', None),
        ('ref_duplicate', 'cpu,n4,n4', None, 'duplicate name', None),
        ('ref_without_cpu', 'nvh,n4', None, 'must contain cpu', None),
        ('ref_member_not_in_run', 'cpu,nvh,n999', None, 'n999 is missing', None),
        ('ref_def_on_LOG_run', 'def', None, 'n1608 is missing', 'm_log_dev1'),
    ]:
        if which:
            q2 = os.path.join(base(which), 'qref.txt')
            rc, msg = py('t207_d5_agg.py', os.path.join(ST, 'x7.json'), 'pTT', '--ref', refs, 'main:a=' + q2)
        else:
            rc, msg = py('t207_d5_agg.py', os.path.join(ST, 'x7.json'), 'pTT', '--ref', refs, 'main:a=' + q)
        record(name, 'error: ' + substr, rc == 1 and substr in msg, msg.strip().split('\n')[-1][:110])
    rc, msg = py('t207_d5_agg.py', os.path.join(ST, 'x8.json'), 'pTT', 'main:a=' + q)
    record('ref_option_omitted', 'error: usage', rc == 1 and 'usage' in msg, msg.strip().split('\n')[-1][:110])
    rc, msg = py('t207_d5_agg.py', os.path.join(ST, 'x9.json'), 'pTT', '--ref', 'def', 'main:a=' + q, 'main:a=' + q)
    record('duplicate_label', 'error: duplicate label', rc == 1 and 'duplicate label' in msg, msg.strip().split('\n')[-1][:110])
    # LOG path baseline with its own list; and the same data refused with the default list (above)
    out = os.path.join(ST, 'agg_log1.json')
    rc, msg = py('t207_d5_agg.py', out, 'pTT', '--ref', 'log', 'main:dev1=' + os.path.join(base('m_log_dev1'), 'qref.txt'))
    record('log_run_with_log_list', 'ok', rc == 0, 'agg accepted' if rc == 0 else msg.strip().split('\n')[-1][:110])


def _swap_members(line, a, b):
    segs = line.split(' | ')
    i = [j for j, s in enumerate(segs) if j and s.split()[0] == a][0]
    k = [j for j, s in enumerate(segs) if j and s.split()[0] == b][0]
    segs[i], segs[k] = segs[k], segs[i]
    return ' | '.join(segs)


def _dup_member(line, a):
    segs = line.split(' | ')
    i = [j for j, s in enumerate(segs) if j and s.split()[0] == a][0]
    segs.insert(i, segs[i])
    return ' | '.join(segs)


def gate_tests():
    out = os.path.join(ST, 'agg_base.json')
    rc, msg = py('t207_d5_agg.py', *agg_args(out))
    record('baseline_agg_def_dev1-4', 'ok', rc == 0, 'agg accepted all four runs' if rc == 0 else msg.strip().split('\n')[-1][:110])
    if rc != 0:
        return
    rc, g = py('t207_d5_gate.py', '--bound', str(BOUND), '--require-controls', out)
    record('baseline_gate', 'GATE PASS', rc == 0 and 'GATE PASS' in g, ' '.join(l.strip() for l in g.strip().split('\n') if 'S_ens' in l or 'GATE' in l)[:160])
    agg = json.load(open(out))

    def variant(name, fn, expect_rc, substr, bound=BOUND, extra=()):
        a = json.loads(json.dumps(agg))
        r = fn(a)
        p = os.path.join(ST, 'gate_%s.json' % name)
        with open(p, 'w') as f:
            json.dump(a, f)
        rc, g = py('t207_d5_gate.py', '--bound', str(bound), *(list(extra) + [p]))
        record(name, 'rc %d: %s' % (expect_rc, substr), rc == expect_rc and substr in g, g.strip().split('\n')[-1][:120])

    variant('gate_bound_below_device', lambda a: None, 1, 'GATE FAIL', bound=-0.5)
    variant('gate_bound_above_controls', lambda a: None, 1, 'NOT KILLED', bound=0.5)
    variant('gate_dev_mean_NaN', lambda a: a['pooled']['dev']['S_ens'].__setitem__('mean', float('nan')), 2, 'INPUT_ERROR')
    variant('gate_wrong_measure', lambda a: a.__setitem__('measure', 'emax'), 2, 'INPUT_ERROR')
    variant('gate_n_used_short', lambda a: a['inputs']['dev1'].__setitem__('n_used', a['inputs']['dev1']['n_used'] - 1), 2, 'INPUT_ERROR')
    variant('gate_pooled_n_wrong', lambda a: a['pooled']['dev']['S_ens'].__setitem__('n', 3), 2, 'INPUT_ERROR')
    variant('gate_ref_missing_in_run', lambda a: a['inputs']['dev2']['members'].remove('n4'), 2, 'INPUT_ERROR')
    variant('gate_no_dev_in_main', lambda a: a['inputs']['dev2']['members'].remove('dev'), 2, 'INPUT_ERROR')
    variant('gate_no_controls_required', lambda a: (a['pooled'].pop('var8'), None), 2, 'INPUT_ERROR', extra=['--require-controls'])
    variant('gate_runs_inputs_differ', lambda a: a['inputs'].pop('dev3'), 2, 'INPUT_ERROR')
    variant('gate_variant_mixed', lambda a: a['inputs']['dev3'].__setitem__('variant', 'log'), 2, 'INPUT_ERROR')
    rc, g = py('t207_d5_gate.py', '--bound', str(BOUND), os.path.join(ST, 'does_not_exist.json'))
    record('gate_missing_file', 'rc 2: INPUT_ERROR', rc == 2 and 'INPUT_ERROR' in g, g.strip()[-100:])
    with open(os.path.join(ST, 'gate_notjson.json'), 'w') as f:
        f.write('{"measure": "pTT", ')
    rc, g = py('t207_d5_gate.py', '--bound', str(BOUND), os.path.join(ST, 'gate_notjson.json'))
    record('gate_truncated_json', 'rc 2: INPUT_ERROR', rc == 2 and 'INPUT_ERROR' in g, g.strip()[-100:])
    rc, g = py('t207_d5_gate.py', '--bound', 'nan', out)
    record('gate_bound_nan', 'rc 2: INPUT_ERROR', rc == 2 and 'INPUT_ERROR' in g, g.strip()[-100:])


# ---------------------------------------------------------------- dump mutations
def dump_layout(path, log):
    with open(path, 'rb') as f:
        h = struct.unpack('<13i', f.read(52))
    n, nv, nvh, dg, hasb = h[7], h[9], h[10], h[11], h[12]
    base_ = 52 + 4 * nv
    sb = 32 * n * n + (8 if log else 16) * n          # one state: U, D (default) or L (LOG), V
    off = dict(U1=base_, D1=base_ + 16 * n * n, V1=base_ + 16 * n * n + (8 if log else 16) * n, S2=base_ + sb)
    off['U2'] = off['S2']; off['V2'] = off['S2'] + 16 * n * n + (8 if log else 16) * n
    off['dev'] = base_ + 2 * sb
    off['var'] = off['dev'] + 64 * n * n
    off['nvh'] = off['dev'] + 64 * n * n * (1 + nv)
    size = os.path.getsize(path)
    off['B'] = size - 32 * n * n
    off['Bi'] = size - 16 * n * n
    off['end'] = size
    return h, n, off


NAN = struct.pack('<d', float('nan'))
INF = struct.pack('<d', float('inf'))


def put(path, pos, b):
    with open(path, 'r+b') as f:
        f.seek(pos)
        f.write(b)


def dump_tests():
    prog = {'def': os.path.join(W, 'bin', 'qref_pp.out'), 'log': os.path.join(W, 'bin', 'qref_log.out')}
    NU = {'def': '0,4', 'log': '0,4'}
    muts = [
        # name, path, expected substring of the QREFERR reason, function(dump path, layout)
        ('B_NaN', 'def', 'non-finite number in B', lambda p, n, o: put(p, o['B'] + 16 * 5, NAN)),
        ('Bi_Inf', 'def', 'non-finite number in Bi', lambda p, n, o: put(p, o['Bi'] + 8 * 7, INF)),
        ('right_U_NaN', 'def', 'non-finite U in the right state', lambda p, n, o: put(p, o['U1'] + 8, NAN)),
        ('right_D_Inf', 'def', 'non-finite D in the right state', lambda p, n, o: put(p, o['D1'] + 16 * 3, INF)),
        ('left_V_NaN', 'def', 'non-finite V in the left state', lambda p, n, o: put(p, o['V2'] + 8 * 100, NAN)),
        ('dev_output_NaN', 'def', 'device output (dev)', lambda p, n, o: put(p, o['dev'] + 8 * 11, NAN)),
        ('var8_output_Inf', 'def', 'device variant var8', lambda p, n, o: put(p, o['var'] + 8 * 11, INF)),
        ('nvh_output_NaN', 'def', 'host output (nvh)', lambda p, n, o: put(p, o['nvh'] + 8 * 11, NAN)),
        ('left_V_zeroed_singular', 'def', 'non-finite number in the quad reference', lambda p, n, o: put(p, o['V2'], b'\0' * (16 * n * n))),
        ('log_L_NaN', 'log', 'non-finite L in the right state', lambda p, n, o: put(p, o['D1'] + 8 * 2, NAN)),
        ('truncated_100_bytes', 'def', 'header implies', lambda p, n, o: os.truncate(p, o['end'] - 100)),
        ('extra_16_bytes', 'def', 'header implies', lambda p, n, o: open(p, 'ab').write(b'\0' * 16)),
        ('bad_magic', 'def', 'bad dump header', lambda p, n, o: put(p, 0, struct.pack('<i', 1))),
        ('header_n_145', 'def', 'N differs from the first dump', lambda p, n, o: put(p, 28, struct.pack('<i', 145))),
        ('header_path_flag_flipped', 'def', 'does not match this build', lambda p, n, o: put(p, 32, struct.pack('<i', 1))),
        ('header_hasb_cleared', 'def', 'header implies', lambda p, n, o: put(p, 48, struct.pack('<i', 0))),
        ('empty_file', 'def', 'short header', lambda p, n, o: os.truncate(p, 0)),
    ]
    # directories: records of m_<p>_dev2 (cover 40: dumps 1, 41, ..., 241), 7 dumps, copies
    def mkdumpdir(name, p, mutate=None, k=3):
        src = base('m_%s_dev2' % p)
        d = os.path.join(ST, 'dumpruns', 'dump_' + name, 'hub12x12_cb_p')
        os.makedirs(os.path.join(d, 'qref'))
        for f in ('parameters', 'info', 'run.log', 'stamp.txt', 'provenance.json'):
            shutil.copy2(os.path.join(src, f), os.path.join(d, f))
        files = sorted(os.listdir(os.path.join(src, 'qref')))
        assert len(files) == 7, files
        for f in files:
            shutil.copyfile(os.path.realpath(os.path.join(src, 'qref', f)), os.path.join(d, 'qref', f))
        if mutate:
            fp = os.path.join(d, 'qref', files[k - 1])
            h, n, off = dump_layout(fp, p == 'log')
            mutate(fp, n, off)
        return d, files

    def proc(d, p, keep, qs='0', nul='0,4', cover='40', hostin='0', symlinks=False):
        env = dict(os.environ, D5W=W, D5RUNS=os.path.join(ST, 'dumpruns'), DUMP_COVER=cover, TMPDIR=os.path.join(W, 'tmp'))
        return sh(['bash', os.path.join(HERE, 't207_d5_proc.sh'), d, p, qs, nul, hostin, str(keep)], env=env)

    # baseline: processed, checked, thinned (KEEP 2: dumps 1, 81, 161, 241 stay; 41, 121, 201 go)
    d, files = mkdumpdir('baseline', 'def')
    rc, out = proc(d, 'def', 2)
    left = sorted(os.listdir(os.path.join(d, 'qref')))
    want = [files[0], files[2], files[4], files[6]]
    thinned = open(os.path.join(d, 'qref_thinned.txt')).read().split() if os.path.exists(os.path.join(d, 'qref_thinned.txt')) else None
    ok = rc == 0 and left == want and thinned == ['qref/' + files[i] for i in (1, 3, 5)]
    record('dump_baseline_and_thinning', 'rc 0, 3 of 7 dumps thinned', ok, 'rc %d, left %d dumps, thinned %s' % (rc, len(left), thinned))
    chk = open(os.path.join(d, 'qref_check.txt')).read().strip().split('\n')
    record('dump_baseline_check_after_thin', 'check --dumps ok', 'dumps in qref/ match the manifest' in chk[-1], chk[-1][:110])
    dd = os.path.join(ST, 'dumpruns', 'dump_baseline', 'hub12x12_cb_p')
    # re-processing a processed directory is refused
    rc, out = proc(dd, 'def', 2)
    record('dump_reprocess_refused', 'rc 1: exists', rc == 1 and 'exists' in out, out.strip().split('\n')[-1][:100])
    # a run directory outside the designated tree is refused
    d, files = mkdumpdir('outside_runs_tree', 'def')
    rc, out = sh(['bash', os.path.join(HERE, 't207_d5_proc.sh'), d, 'def', '0', '0,4', '0', '0'], env=dict(os.environ, D5W=W, D5RUNS=os.path.join(W, 'runs'), DUMP_COVER='40'))
    record('dump_proc_refuses_dir_outside_runs', 'rc 1: refusing', rc == 1 and 'refusing' in out and not os.path.exists(os.path.join(d, 'qref_spec.json')), out.strip().split('\n')[-1][:100])
    # the dump mutations
    for name, p, substr, fn in muts:
        d, files = mkdumpdir(name, p, fn)
        fp = os.path.join(d, 'qref', files[2])
        rc, out = proc(d, p, 2)
        txt = open(os.path.join(d, 'qref.txt')).read() if os.path.exists(os.path.join(d, 'qref.txt')) else ''
        err = [l for l in txt.split('\n') if l.startswith('QREFERR')]
        endl = [l for l in txt.split('\n') if l.startswith('QREFEND')]
        chk = open(os.path.join(d, 'qref_check.txt')).read() if os.path.exists(os.path.join(d, 'qref_check.txt')) else ''
        prc = json.load(open(os.path.join(d, 'qref_proc.json'))) if os.path.exists(os.path.join(d, 'qref_proc.json')) else {}
        left = os.listdir(os.path.join(d, 'qref'))
        ok = (rc == 3 and len(err) == 1 and files[2] in err[0] and substr in err[0] and prc.get('rc') == 3 and prc.get('check_rc') == 1
              and endl == ['QREFEND dumps 7 ok 6 err 1'] and len(left) == 7 and 'QPARSE ERROR' in chk)
        record('dump_' + name, 'QREFERR, rc 3, nothing deleted', ok, 'rc %d; %s; left %d' % (rc, (err[0] if err else out.strip().split('\n')[-1])[-100:], len(left)))
    # a version 2 dump without B, Bi (hasb cleared, B and Bi cut off): the program must print the propagated numbers as -1 and not fail; the check
    # refuses it because the run's parameters say that this point (nst 1 < NSTM) has B
    d, files = mkdumpdir('B_dropped', 'def', lambda p, n, o: (os.truncate(p, o['end'] - 32 * n * n), put(p, 48, struct.pack('<i', 0))))
    rc, out = proc(d, 'def', 2)
    txt = open(os.path.join(d, 'qref.txt')).read()
    prc = json.load(open(os.path.join(d, 'qref_proc.json')))
    chk = open(os.path.join(d, 'qref_check.txt')).read()
    l3 = [l for l in txt.split('\n') if l.startswith('QREF file %d ' % int(files[2][1:7]))]
    ok = rc == 3 and prc['rc'] == 0 and prc['check_rc'] == 1 and 'header hasb = 0, expected 1' in chk and len(l3) == 1 and l3[0].count('-1.0000E+000') == 3 * 8 and len(os.listdir(os.path.join(d, 'qref'))) == 7
    record('dump_B_dropped_hasb_cleared', 'program rc 0, check error, kept', ok, 'rc %d, program rc %s, %s' % (rc, prc['rc'], chk.strip()[-80:]))
    # directory checks
    d, files = mkdumpdir('extra_file', 'def')
    shutil.copyfile(os.path.join(d, 'qref', files[0]), os.path.join(d, 'qref', 'q000002.bin'))
    rc, out = proc(d, 'def', 2)
    record('dump_extra_file_in_qref', 'rc 2: unexpected', rc == 2 and 'unexpected' in out and not os.path.exists(os.path.join(d, 'qref.txt')), out.strip().split('\n')[-1][:100])
    d, files = mkdumpdir('missing_file', 'def')
    os.remove(os.path.join(d, 'qref', files[4]))
    rc, out = proc(d, 'def', 2)
    record('dump_missing_file_in_qref', 'rc 2: missing', rc == 2 and 'missing' in out, out.strip().split('\n')[-1][:100])
    d, files = mkdumpdir('wrong_cover', 'def')
    rc, out = proc(d, 'def', 2, cover='20')
    record('dump_wrong_cover_argument', 'rc 2', rc == 2 and 'missing' in out, out.strip().split('\n')[-1][:100])
    d, files = mkdumpdir('wrong_path', 'def')
    rc, out = proc(d, 'log', 2)
    record('dump_wrong_path_argument', 'rc 2: run is a def run', rc == 2 and 'is a def run' in out, out.strip().split('\n')[-1][:100])
    d, files = mkdumpdir('stale_manifest', 'def')
    sh(['bash', '-c', 'cd %s && sha256sum qref/q*.bin > qref_manifest.txt && echo 0000000000000000000000000000000000000000000000000000000000000000  qref/q000001.bin >> qref_manifest.txt' % d])
    rc, out = proc(d, 'def', 2)
    record('dump_stale_manifest', 'rc 2: manifest differs', rc == 2 and 'manifest' in out, out.strip().split('\n')[-1][:100])
    # symbolic links are not thinned
    d, files = mkdumpdir('symlinks', 'def')
    for f in files:
        fp = os.path.join(d, 'qref', f)
        tgt = os.path.join(ST, 'dumpruns', 'dump_symlinks', 'orig_' + f)
        shutil.move(fp, tgt)
        os.symlink(tgt, fp)
    rc, out = proc(d, 'def', 2)
    left = os.listdir(os.path.join(d, 'qref'))
    record('dump_symlinks_not_thinned', 'rc 4, nothing deleted', rc == 4 and len(left) == 7 and 'symbolic link' in out + open(os.path.join(d, 'qref_check.txt')).read() + out,
           'rc %d, left %d: %s' % (rc, len(left), out.strip().split('\n')[-1][:90]))
    # thinning guard: a dump changed after the manifest was written is not deleted (spec + manifest by hand, no program run)
    d, files = mkdumpdir('thin_sha_guard', 'def')
    env = dict(os.environ, D5W=W)
    rc, out = py('t207_d5_qparse.py', 'spec', d, '--path', 'def', '--qs', '0', '--nul', '0,4', '--hostin', '0', '--cover', '40', '--keep', '2')
    sh(['bash', '-c', 'cd %s && sha256sum qref/q*.bin > qref_manifest.txt' % d])
    with open(os.path.join(d, 'qref', files[1]), 'ab') as f:
        f.write(b'x')
    rc, out = py('t207_d5_qparse.py', 'thin', d)
    record('thin_guard_sha_mismatch', 'rc 1, nothing deleted', rc == 1 and len(os.listdir(os.path.join(d, 'qref'))) == 7, out.strip().split('\n')[-1][:100])
    # the LOG path baseline (one run through the program with the LOG build)
    d, files = mkdumpdir('log_baseline', 'log')
    rc, out = proc(d, 'log', 0)
    record('dump_log_baseline', 'rc 0', rc == 0, out.strip().split('\n')[-1][:100])


def main():
    if os.path.exists(ST):
        sys.exit('refusing: %s exists' % ST)
    os.makedirs(ST)
    if WHAT in ('text', 'all'):
        gate_tests()
        text_tests()
    if WHAT in ('dump', 'all'):
        dump_tests()
    bad = [r for r in RESULTS if not r[2]]
    print('\nself-test: %d tests, %d with the expected outcome, %d WRONG' % (len(RESULTS), len(RESULTS) - len(bad), len(bad)))
    with open(os.path.join(ST, 'results.txt'), 'w') as f:
        for name, expect, ok, detail in RESULTS:
            f.write('%s\t%s\t%s\t%s\n' % (name, expect, 'ok' if ok else 'WRONG', detail.replace('\n', ' ')))
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
