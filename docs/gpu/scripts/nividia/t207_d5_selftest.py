#!/usr/bin/env python3
"""Ticket 207-2: mutation self-tests of the QREF analysis pipeline (nividia, CPU only).   usage: t207_d5_selftest.py WORKDIR [text|dump|all]

Needs WORKDIR/runs/{m_def_dev1..4,c_def_dev1,d_def_dev1,m_log_dev1,s_def,s_log}/hub12x12_cb_p produced by t207_d5_repro.sh (the new pipeline on the kept dev5
dumps) and WORKDIR/bin/qref_pp.out, qref_log.out.  Everything it writes goes to WORKDIR/st/ (created; refuses to start if it exists).
  baseline   the unmutated inputs must be accepted (agg exit 0) and the gate must say PASS (device S_ens(pTT) pooled over the four
             default-path seed sets <= bound, both controls > bound)
  text       mutations of qref.txt, of the run's records and of the spec: each must end in the error named (the strict parser, agg,
             stat) or, for the swap of device and var8, in GATE FAIL
  gate       mutations of the aggregate JSON, of the gate's arguments and of the bound: INPUT_ERROR / FAIL as named.  The kept dev5 dumps are a
             subset (cover 20 / 40), which the gate refuses (cover != 1, tested unmodified); the gate logic is tested on the same aggregate with
             inputs[*].args.cover set to 1 ('legalised', a doctored aggregate: the data are the kept subset, the tests are about the gate's
             decisions, one planted defect at a time against that PASSing baseline)
  dump       mutations of real dumps (NaN, Infinity, zeroed state, truncation, extra bytes, bad header) run through t207_d5_proc.sh
             (program + strict check): the program must print a QREFERR line for the dump, exit with status 3, the check must fail and no
             dump may be deleted; plus the directory checks (extra / missing dump, stale manifest, symlink and sha256 guards of the thinning).
             The dev5 v2 dumps are kept at an even spacing (cover 20/40), which the keep rule refuses for thinning, so the thinning is tested
             on the cover-1 version 1 smoke runs s_def / s_log (stride 1, KEEP 3); the v2 mutations run with KEEP 0.
Output: one table line per test, generated from the actual outcomes; exit status 1 if any test did not give the expected outcome.
The bound used is 0.15 (the proposed criterion of 207 dev5; the self-test is not a calibration and does not change any bound).
"""
import hashlib, importlib.util, json, os, re, shutil, struct, subprocess, sys
sys.dont_write_bytecode = True

W = os.path.abspath(sys.argv[1])
WHAT = sys.argv[2] if len(sys.argv) > 2 else 'all'
HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
_sp = importlib.util.spec_from_file_location('qp', os.path.join(HERE, 't207_d5_qparse.py'))
QP = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(QP)
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


REC = ('parameters', 'info', 'run.log', 'stamp.txt', 'provenance.json', 'qref_spec.json', 'qref_proc.json')


def mkrun(name, src_tag, qtxt=None, rec=None, spec=None, fixprov=False):
    """a run directory st/name/hub12x12_cb_p: the records of src_tag, qref.txt = qtxt (default: the original), optional edits of records"""
    d = os.path.join(ST, name, 'hub12x12_cb_p')
    os.makedirs(d)
    for f in REC:
        if os.path.exists(os.path.join(base(src_tag), f)):
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
        for key, f in (('parameters_sha256', 'parameters'), ('run_log_sha256', 'run.log'), ('info_sha256', 'info')):
            pv[key] = hashlib.sha256(open(os.path.join(d, f), 'rb').read()).hexdigest()
        st_ = open(os.path.join(d, 'stamp.txt')).read()               # the stamp carries the first 16 digits of the parameters hash
        st2 = re.sub(r'( parameters )[0-9a-f]{16}', lambda m: m.group(1) + pv['parameters_sha256'][:16], st_)
        open(os.path.join(d, 'stamp.txt'), 'w').write(st2)
        json.dump(pv, open(pj, 'w'), indent=1)
        sp = os.path.join(d, 'qref_spec.json')
        st = json.load(open(sp))                         # the stored spec carries the record hashes: keep it consistent as well
        st['run']['records_sha256'] = {r: hashlib.sha256(open(os.path.join(d, r), 'rb').read()).hexdigest() for r in QP.RECORDS}
        json.dump(st, open(sp, 'w'), indent=1)
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
G_SET = ['dev2', 'dev3', 'dev4']          # the gate tests use the three runs with the same number of points (7 with B each, cover 40)


def legal(a):
    """a doctored copy of an aggregate: every run claims cover 1 and a null-error-statistic qs of 20 (the kept dev5 dumps are a subset, which
    the gate refuses unmodified).  All other fields stay as the aggregation wrote them."""
    a = json.loads(json.dumps(a))
    for i in a['inputs'].values():
        i['args']['cover'] = 1
        i['args']['qs'] = 20
    return a


def gate_run(agg, name, expect=7, qs=(20,), bound=None, extra=(), args=None):
    """write the aggregate dict and run the gate on it; returns (rc, output)"""
    p = os.path.join(ST, 'gate_%s.json' % name)
    with open(p, 'w') as f:
        json.dump(agg, f)
    a = ['--bound', str(BOUND if bound is None else bound), '--expect-points', str(expect)] if args is None else list(args)
    if args is None:
        for q in qs:
            a += ['--qs', str(q)]
    return py('t207_d5_gate.py', *(a + list(extra) + [p]))


def agg_args(out, mut_path=None, mut_label='dev1', refs='def', measure='pTT', extra=None, sets=None):
    a = [out, measure, '--ref', refs] + (extra or [])
    for s in (sets or A_SET):
        p = mut_path if (s == mut_label and mut_path) else os.path.join(base('m_def_' + s), 'qref.txt')
        a.append('main:%s=%s' % (s, p))
    return a


def expect_error(name, qpath, substr, label='dev1', absent=None, **kw):
    out = os.path.join(ST, 'agg_%s.json' % name)
    rc, msg = py('t207_d5_agg.py', *agg_args(out, qpath, label, **kw))
    ok = rc == 1 and substr in msg and not os.path.exists(out) and (absent is None or absent not in msg)
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
        ('QREFERR_line', lambda t: edit_line(t, 2, lambda l: 'QREFERR dump qref/q000041.bin : non-finite number in B'),
         'the program reported an error: QREFERR dump qref/q000041.bin : non-finite number in B'),
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
        expect_error(name, q, substr, absent='unexpected line' if name == 'QREFERR_line' else None)
    # zero error (well-formed): the parser accepts it, agg refuses to take a log
    q = mkrun('zero_error_dev_pTT', 'm_def_dev1', edit_line(t0, 1, set_tok('dev', 5, '0.0000E+000')))
    expect_error('zero_error_dev_pTT', q, 'not positive and finite')
    # well-formed huge value (the old clamp value): accepted as data (a very large error, no longer a placeholder): the pooled S_ens moves
    q = mkrun('overflow_value_dev_pTT', 'm_def_dev1', edit_line(t0, 1, set_tok('dev', 5, '1.7977E+308')))
    expect_error('overflow_value_dev_pTT', q, 'is not finite')
    q = mkrun('huge_value_dev_pTT', 'm_def_dev2', edit_line(orig('m_def_dev2'), 1, set_tok('dev', 5, '1.0000E+100')))
    out = os.path.join(ST, 'agg_huge.json')
    rc, msg = py('t207_d5_agg.py', *agg_args(out, q, 'dev2', sets=G_SET))
    rc2, g = gate_run(legal(json.load(open(out))), 'huge') if rc == 0 else (-1, msg)
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
    sets = {}
    for lab in G_SET:
        sets[lab] = mkrun('swap_dev_var8_' + lab, 'm_def_' + lab, swap(orig('m_def_' + lab)))
    out = os.path.join(ST, 'agg_swap.json')
    rc, msg = py('t207_d5_agg.py', out, 'pTT', '--ref', 'def', *['main:%s=%s' % (s, sets[s]) for s in G_SET])
    rc2, g = gate_run(legal(json.load(open(out))), 'swap') if rc == 0 else (-1, msg)
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
    rc, msg = py('t207_d5_stat.py', os.path.join(ST, 'agg_d_base.json'), 'din=' + q)
    record('d_def_baseline_stat', 'ok', rc == 0, 'stat accepted the run with din/cpuh/cpuhx members' if rc == 0 else msg.strip().split('\n')[-1][:100])
    rc, msg = py('t207_d5_agg.py', os.path.join(ST, 'agg_d_env.json'), 'pTT', '--ref', 'def', 'main:d=' + q)
    record('d_def_run_as_main_env_refused', 'error: env allow-list', rc == 1 and 'env allow-list' in msg and 'ALF_GPU_TAUDIAG=2' in msg, msg.strip().split('\n')[-1][:110])
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


def mkenv(name, src, fn):
    """a copy of run src whose env (ALF_GPU_*) is changed by fn(dict) consistently in info, stamp.txt, provenance.json and the stored spec"""
    q = mkrun(name, src)
    d = os.path.dirname(q)
    sp = json.load(open(os.path.join(d, 'qref_spec.json')))
    new = fn(dict(sp['run']['env']))
    t = open(os.path.join(d, 'info')).read()
    t2 = re.sub(r'^( *GPU debug env: ).*$', lambda m: m.group(1) + ', '.join('%s=%s' % kv for kv in sorted(new.items())), t, flags=re.M)
    assert t2 != t
    open(os.path.join(d, 'info'), 'w').write(t2)
    t = open(os.path.join(d, 'stamp.txt')).read()
    t2 = re.sub(r'^env .*$', lambda m: 'env ' + ' '.join('%s=%s' % kv for kv in new.items()), t, flags=re.M)
    assert t2 != t
    open(os.path.join(d, 'stamp.txt'), 'w').write(t2)
    pv = json.load(open(os.path.join(d, 'provenance.json')))
    pv['env'] = ['%s=%s' % kv for kv in new.items()]
    pv['info_sha256'] = hashlib.sha256(open(os.path.join(d, 'info'), 'rb').read()).hexdigest()
    json.dump(pv, open(os.path.join(d, 'provenance.json'), 'w'), indent=1)
    sp['run']['env'] = new
    sp['run']['records_sha256'] = {r: hashlib.sha256(open(os.path.join(d, r), 'rb').read()).hexdigest() for r in QP.RECORDS}
    json.dump(sp, open(os.path.join(d, 'qref_spec.json'), 'w'), indent=1)
    return q


def agg_direct(name, refs, runs, extra=()):
    out = os.path.join(ST, 'agg_%s.json' % name)
    rc, msg = py('t207_d5_agg.py', out, 'pTT', '--ref', refs, *(list(extra) + ['%s:%s=%s' % r for r in runs]))
    return rc, msg, out


def last(msg):
    return msg.strip().split('\n')[-1][:110]


def agg_tests():
    """role / env allow-list, duplicates, program hash, keep rule (review fixes 2, 3, 5, 7)"""
    Q = lambda tag: os.path.join(base(tag), 'qref.txt')                                         # noqa: E731
    # the baseline with a control run: the control run has the same seeds as dev1 (another role: allowed)
    rc, msg, out = agg_direct('ctl_base', 'def', [('main', s, Q('m_def_' + s)) for s in A_SET] + [('ctl', 'c1', Q('c_def_dev1'))])
    record('agg_with_ctl_run_baseline', 'ok', rc == 0, 'agg accepted 4 main runs and the control run (same seeds as dev1)' if rc == 0 else last(msg))
    # env allow-list and role
    q = mkenv('env_TAUSWAP_on_def_main', 'm_def_dev2', lambda e: dict(e, ALF_GPU_TAUSWAP='8'))
    rc, msg, out = agg_direct('x_env1', 'def', [('main', 'a', q)])
    record('env_TAUSWAP_on_default_path_main', 'error: env allow-list', rc == 1 and 'env allow-list' in msg and 'TAUSWAP=8' in msg, last(msg))
    q = mkenv('env_TAUSWAP_mask_log_main', 'm_log_dev1', lambda e: dict(e, ALF_GPU_TAUSWAP='8'))
    rc, msg, out = agg_direct('x_env2', 'log', [('main', 'a', q)])
    record('env_TAUSWAP_on_LOG_run_given_as_main', 'error: role main given', rc == 1 and 'role main given' in msg, last(msg))
    rc, msg, out = agg_direct('x_env2b', 'log', [('ctl', 'a', q)])
    record('env_TAUSWAP_8_on_LOG_run_as_ctl_ok', 'ok', rc == 0, 'agg accepted the LOG run with TAUSWAP=8 as ctl' if rc == 0 else last(msg))
    q = mkenv('env_unknown_var_main', 'm_def_dev2', lambda e: dict(e, ALF_GPU_FOO='1'))
    rc, msg, out = agg_direct('x_env3', 'def', [('main', 'a', q)])
    record('env_unknown_ALF_GPU_variable_main', 'error: env allow-list', rc == 1 and 'env allow-list' in msg and 'ALF_GPU_FOO' in msg, last(msg))
    q = mkenv('env_TAU_NS_1_ctl', 'c_def_dev1', lambda e: dict(e, ALF_GPU_TAU_NS='1'))
    rc, msg, out = agg_direct('x_env4', 'def', [('ctl', 'a', q)])
    record('env_TAU_NS_1_ctl', 'error: env allow-list', rc == 1 and 'env allow-list' in msg, last(msg))
    q = mkenv('env_TAU_NS_removed_ctl', 'c_def_dev1', lambda e: {k: v for k, v in e.items() if k != 'ALF_GPU_TAU_NS'})
    rc, msg, out = agg_direct('x_env5', 'def', [('ctl', 'a', q)])
    record('env_ctl_without_TAU_NS_has_no_QVAR', 'error: env allow-list', rc == 1 and 'env allow-list' in msg, last(msg))
    rc, msg, out = agg_direct('x_env6', 'def', [('main', 'a', Q('c_def_dev1'))])
    record('role_main_for_control_run', 'error: role main given', rc == 1 and 'role main given' in msg, last(msg))
    rc, msg, out = agg_direct('x_env7', 'def', [('ctl', 'a', Q('m_def_dev2'))])
    record('role_ctl_for_main_run', 'error: role ctl given', rc == 1 and 'role ctl given' in msg, last(msg))
    q = mkenv('env_QHOST_0', 'm_def_dev2', lambda e: dict(e, ALF_GPU_TAUDIAG_QHOST='0'))
    rc, msg, out = agg_direct('x_env8', 'def', [('main', 'a', q)])
    record('env_QHOST_0_main', 'error (any)', rc == 1, last(msg))
    # duplicates
    rc, msg, out = agg_direct('x_dup1', 'def', [('main', 'a', Q('m_def_dev2')), ('main', 'b', Q('m_def_dev2'))])
    record('dup_same_run_two_labels', 'error: same run', rc == 1 and 'is the same run as a' in msg and 'qref.txt sha256' in msg, last(msg))
    q2 = mkrun('dup_copied_dir', 'm_def_dev2')
    rc, msg, out = agg_direct('x_dup2', 'def', [('main', 'a', Q('m_def_dev2')), ('main', 'b', q2)])
    record('dup_copied_run_directory', 'error: same run', rc == 1 and 'qref.txt sha256' in msg, last(msg))
    t2 = orig('m_def_dev2')
    q2 = mkrun('dup_edited_qref', 'm_def_dev2', edit_line(t2, 1, set_tok('u1', 0, '1.2345E-012')))
    rc, msg, out = agg_direct('x_dup3', 'def', [('main', 'a', Q('m_def_dev2')), ('main', 'b', q2)])
    record('dup_other_qref_same_records', 'error: same run.log', rc == 1 and 'run.log sha256' in msg, last(msg))
    q2 = mkrun('dup_same_seeds', 'm_def_dev2', edit_line(t2, 1, set_tok('u1', 0, '1.2345E-012')), rec={'run.log': lambda t: t + ' harmless extra line\n'}, fixprov=True)
    rc, msg, out = agg_direct('x_dup4', 'def', [('main', 'a', Q('m_def_dev2')), ('main', 'b', q2)])
    record('dup_same_seeds_and_parameters', 'error: same run record', rc == 1 and 'role, seeds, parameters' in msg, last(msg))
    rc, msg, out = agg_direct('x_dup5', 'def', [('main', 'a', Q('m_def_dev2')), ('main', 'b', Q('m_def_dev3'))])
    record('dup_two_different_runs_ok', 'ok', rc == 0, 'agg accepted dev2 and dev3' if rc == 0 else last(msg))
    # program hash
    prog = json.load(open(os.path.join(base('m_def_dev2'), 'qref_proc.json')))['prog']
    rc, msg, out = agg_direct('x_prog1', 'def', [('main', s, Q('m_def_' + s)) for s in G_SET], extra=['--prog-sha', prog])
    ag = json.load(open(out)) if rc == 0 else {}
    record('prog_sha_matching', 'ok, recorded in the aggregate', rc == 0 and all(i['prog_sha256'] == prog for i in ag['inputs'].values()),
           'agg accepted --prog-sha %s' % prog[:16] if rc == 0 else last(msg))
    rc, msg, out = agg_direct('x_prog2', 'def', [('main', 'a', Q('m_def_dev2'))], extra=['--prog-sha', '0' * 64])
    record('prog_sha_mismatch', 'error: differs from --prog-sha', rc == 1 and 'differs from --prog-sha' in msg, last(msg))
    q2 = mkrun('prog_no_proc_json', 'm_def_dev2')
    os.remove(os.path.join(os.path.dirname(q2), 'qref_proc.json'))
    rc, msg, out = agg_direct('x_prog3', 'def', [('main', 'a', q2)], extra=['--prog-sha', prog])
    record('prog_sha_asked_proc_json_missing', 'error: is missing', rc == 1 and 'is missing' in msg, last(msg))
    q2 = mkrun('proc_json_rc1', 'm_def_dev2')
    pj = os.path.join(os.path.dirname(q2), 'qref_proc.json')
    d = json.load(open(pj)); d['check_rc'] = 1; json.dump(d, open(pj, 'w'))
    rc, msg, out = agg_direct('x_prog4', 'def', [('main', 'a', q2)])
    record('proc_json_check_rc_1', 'error: not a passing processing', rc == 1 and 'not a passing processing' in msg, last(msg))
    q2 = mkrun('proc_json_args_differ', 'm_def_dev2')
    pj = os.path.join(os.path.dirname(q2), 'qref_proc.json')
    d = json.load(open(pj)); d['args']['qs'] = 40; json.dump(d, open(pj, 'w'))
    rc, msg, out = agg_direct('x_prog5', 'def', [('main', 'a', q2)])
    record('proc_json_args_differ_from_spec', 'error: not a passing processing', rc == 1 and 'not a passing processing' in msg, last(msg))
    # keep rule: (1) brute-force equivalence with the simulated call sequence over many (NSTM, N_FL, stride, step), (2) the old rule planted, (3) the CLI
    bad = n = 0
    for nstm in (1, 2, 3, 4, 5):
        for nfl in (1, 2):
            c = nstm * nfl
            for stride in range(1, 12):
                allpos = {(j * stride) % c for j in range(2 * c)}                        # positions of all dumps in the pass cycle
                for m in range(1, 91):
                    pp, g = QP.keep_rule(nstm, nfl, stride, m)
                    pos = [((m * j * stride) % c) for j in range(4 * pp * c)]            # kept dump j is dump m*j, call m*j*stride
                    cnt = {x: pos.count(x) for x in set(pos)}
                    balanced = set(cnt) == allpos and len(set(cnt.values())) == 1
                    n += 1
                    if (g == 1) != balanced:
                        bad += 1
    record('keep_rule_equals_simulation', '%d cases, 0 mismatches' % n, bad == 0,
           '%d cases (NSTM 1-5, N_FL 1-2, stride 1-11, step 1-90), %d mismatches between gcd rule and simulated even visit of all positions' % (n, bad))
    pp, g = QP.keep_rule(4, 2, 2, 40)
    visited = sorted({(40 * j * 2) % 8 for j in range(100)})
    record('keep_rule_old_every40th_biased', 'refused (gcd 4), visits 1 of 4 positions', g == 4 and visited == [0] and pp == 4,
           'P %d, gcd %d, positions visited by every 40th dump of stride 2 (C 8): %s' % (pp, g, visited))
    pp, g = QP.keep_rule(4, 2, 2, 39)
    visited = sorted({(39 * j * 2) % 8 for j in range(100)})
    record('keep_rule_39_accepted', 'accepted, all 4 positions', g == 1 and visited == [0, 2, 4, 6], 'P %d, gcd %d, positions visited by every 39th dump of stride 2: %s' % (pp, g, visited))
    rd = os.path.join(ST, 'recdir', 'hub12x12_cb_p')
    os.makedirs(rd)
    for f in QP.RECORDS:
        shutil.copy2(os.path.join(base('m_def_dev2'), f), os.path.join(rd, f))
    for (cover, keep, ok) in ((1, 40, False), (1, 39, True), (40, 1, False), (1, 8, False), (3, 13, True), (1, 0, True)):
        try:
            sp = QP.make_spec(rd, 'def', 20, [0, 4], 0, cover, keep)
            got, msg = True, 'accepted'
        except QP.QError as e:
            got, msg = False, str(e)
        record('keep_rule_make_spec_cover%d_keep%d' % (cover, keep), 'accepted' if ok else 'error: keep rule', got == ok and (ok or 'keep rule' in msg), msg[:110])
        if ok and keep == 39:
            summ = QP.kept_summary(sp)
            record('keep_summary_39', '8 of 280, nst 2 each', '8 of 280; per nst 1:2 2:2 3:2 4:2' in summ, summ[:150])
    rc, out = py('t207_d5_qparse.py', 'spec', rd, '--path', 'def', '--qs', '20', '--nul', '0,4', '--hostin', '0', '--cover', '1', '--keep', '40')
    record('keep_rule_cli_spec_refused', 'rc 1: keep rule', rc == 1 and 'keep rule' in out and not os.path.exists(os.path.join(rd, 'qref_spec.json')), last(out))


def gate_tests():
    out = os.path.join(ST, 'agg_base.json')
    rc, msg = py('t207_d5_agg.py', *agg_args(out))
    record('baseline_agg_def_dev1-4', 'ok', rc == 0, 'agg accepted all four runs' if rc == 0 else last(msg))
    if rc != 0:
        return
    # the gate on the kept dev5 subset itself: refused (cover != 1); the gate on dev2-4 + control, legalised: PASS
    rc, msg = py('t207_d5_agg.py', os.path.join(ST, 'agg_gate.json'), 'pTT', '--ref', 'def',
                 *(['main:%s=%s' % (s, os.path.join(base('m_def_' + s), 'qref.txt')) for s in G_SET] + ['ctl:c1=' + os.path.join(base('c_def_dev1'), 'qref.txt')]))
    record('baseline_agg_gate_set', 'ok', rc == 0, 'agg accepted dev2-4 and the control run' if rc == 0 else last(msg))
    agg0 = json.load(open(os.path.join(ST, 'agg_gate.json')))
    rc, g = gate_run(agg0, 'unmodified_kept_subset')
    record('gate_kept_dev5_subset_unmodified', 'rc 2: INPUT_ERROR (cover)', rc == 2 and 'INPUT_ERROR' in g and 'cover' in g and 'GATE FAIL' not in g, g.strip()[-120:])
    agg = legal(agg0)
    rc, g = gate_run(agg, 'baseline', extra=['--require-controls'])
    record('baseline_gate_legalised', 'GATE PASS', rc == 0 and 'GATE PASS' in g, ' '.join(l.strip() for l in g.strip().split('\n') if 'S_ens' in l or 'GATE' in l)[:160])
    prog = agg['inputs']['dev2']['prog_sha256']

    def variant(name, fn, expect_rc, substr, bound=None, extra=(), expect=7, qs=(20,), args=None, a0=None):
        a = json.loads(json.dumps(agg if a0 is None else a0))
        fn(a)
        rc, g = gate_run(a, name, expect=expect, qs=qs, bound=bound, extra=extra, args=args)
        # the line 'GATE FAIL' appears exactly for a real FAIL (rc 1), never for an input error or a crash
        ok = rc == expect_rc and substr in g and (('GATE FAIL' in g) == (expect_rc == 1))
        record(name, 'rc %d: %s' % (expect_rc, substr), ok, g.strip().split('\n')[-1][:120])

    nop = lambda a: None                                                                         # noqa: E731
    variant('gate_bound_below_device', nop, 1, 'GATE FAIL', bound=-0.5)
    variant('gate_bound_above_controls', nop, 1, 'NOT KILLED', bound=0.5)
    variant('gate_dev_mean_NaN', lambda a: a['pooled']['dev']['S_ens'].__setitem__('mean', float('nan')), 2, 'INPUT_ERROR')
    variant('gate_wrong_measure', lambda a: a.__setitem__('measure', 'emax'), 2, 'INPUT_ERROR')
    variant('gate_n_used_short', lambda a: a['inputs']['dev2'].__setitem__('n_used', a['inputs']['dev2']['n_used'] - 1), 2, 'INPUT_ERROR')
    variant('gate_pooled_n_wrong', lambda a: a['pooled']['dev']['S_ens'].__setitem__('n', 3), 2, 'INPUT_ERROR')
    variant('gate_ref_missing_in_run', lambda a: a['inputs']['dev2']['members'].remove('n4'), 2, 'INPUT_ERROR')
    variant('gate_no_dev_in_main', lambda a: a['inputs']['dev2']['members'].remove('dev'), 2, 'INPUT_ERROR')
    variant('gate_no_controls_required', lambda a: (a['pooled'].pop('var8'), a['pooled'].pop('ctl_dev'), None), 2, 'INPUT_ERROR', extra=['--require-controls'])
    variant('gate_runs_inputs_differ', lambda a: a['inputs'].pop('dev3'), 2, 'INPUT_ERROR')
    variant('gate_variant_mixed', lambda a: a['inputs']['dev3'].__setitem__('variant', 'log'), 2, 'INPUT_ERROR')
    # review fix 1: configuration of every run
    variant('gate_cover_2', lambda a: a['inputs']['dev3']['args'].__setitem__('cover', 2), 2, 'cover')
    variant('gate_v1_dumps', lambda a: a['inputs']['dev3'].__setitem__('v1', True), 2, 'version 1')
    variant('gate_v1_flag_missing', lambda a: a['inputs']['dev3'].pop('v1'), 2, "'v1' missing")
    variant('gate_qs_zero', lambda a: a['inputs']['dev3']['args'].__setitem__('qs', 0), 2, 'qs', qs=(0, 20))
    variant('gate_qs_not_passed', lambda a: a['inputs']['dev3']['args'].__setitem__('qs', 40), 2, 'qs 40')
    variant('gate_qs_two_values_passed_ok', lambda a: a['inputs']['dev3']['args'].__setitem__('qs', 40), 0, 'GATE PASS', qs=(20, 40))
    variant('gate_nul_without_0', lambda a: a['inputs']['dev3']['args']['nul'].remove(0), 2, 'null codes')
    variant('gate_nul_misses_ref_code', lambda a: a['inputs']['dev3']['args']['nul'].remove(1608), 2, 'null code 1608')
    variant('gate_refs_member_dropped', lambda a: a['refs'].remove('n1608'), 2, 'frozen list')
    variant('gate_refs_reordered', lambda a: a['refs'].reverse(), 2, 'frozen list')
    variant('gate_refs_other_path_list', lambda a: a.__setitem__('refs', QP.REF_LOG.split(',')), 2, 'frozen list')
    variant('gate_n_used_vs_expect_points', nop, 2, 'points used', expect=8)
    variant('gate_n_used_tampered_with_expected', lambda a: [i.update(n_used=8, n_expected=8) for i in a['inputs'].values()], 2, 'INPUT_ERROR')
    # review fix 2: env and role
    variant('gate_env_TAUSWAP_main', lambda a: a['inputs']['dev3']['env'].__setitem__('ALF_GPU_TAUSWAP', '8'), 2, 'env allow-list')
    variant('gate_env_unknown_var_main', lambda a: a['inputs']['dev3']['env'].__setitem__('ALF_GPU_FOO', '1'), 2, 'env allow-list')
    variant('gate_env_QVAR_removed_main', lambda a: a['inputs']['dev3']['env'].pop('ALF_GPU_TAUDIAG_QVAR'), 2, 'env allow-list')
    variant('gate_env_QVAR_other_main', lambda a: a['inputs']['dev3']['env'].__setitem__('ALF_GPU_TAUDIAG_QVAR', '8,16'), 2, 'env allow-list')
    variant('gate_env_DEBUG_removed_main', lambda a: a['inputs']['dev3']['env'].pop('ALF_GPU_DEBUG'), 2, 'env allow-list')
    variant('gate_env_QHOST_0_main', lambda a: a['inputs']['dev3']['env'].__setitem__('ALF_GPU_TAUDIAG_QHOST', '0'), 2, 'env allow-list')
    variant('gate_env_QREF_not_a_number', lambda a: a['inputs']['dev3']['env'].__setitem__('ALF_GPU_TAUDIAG_QREF', 'x'), 2, 'env allow-list')
    variant('gate_env_ctl_TAU_NS_1', lambda a: a['inputs']['c1']['env'].__setitem__('ALF_GPU_TAU_NS', '1'), 2, 'env allow-list')
    variant('gate_env_ctl_TAUSWAP_on_default_path',
            lambda a: (a['inputs']['c1']['env'].pop('ALF_GPU_TAU_NS'), a['inputs']['c1']['env'].__setitem__('ALF_GPU_TAUSWAP', '8')), 2, 'env allow-list')
    variant('gate_role_main_but_ctl_env', lambda a: a['inputs']['c1'].__setitem__('role', 'main'), 2, 'role main')
    variant('gate_role_ctl_but_main_env', lambda a: a['inputs']['dev3'].__setitem__('role', 'ctl'), 2, 'role ctl')
    variant('gate_env_missing', lambda a: a['inputs']['dev3'].pop('env'), 2, "'env' missing")
    # review fix 3: the same run twice
    variant('gate_dup_qref_sha256', lambda a: a['inputs']['dev3'].__setitem__('qref_sha256', a['inputs']['dev2']['qref_sha256']), 2, 'same run')
    variant('gate_dup_run_log_sha256', lambda a: a['inputs']['dev3'].__setitem__('run_log_sha256', a['inputs']['dev2']['run_log_sha256']), 2, 'same run')
    variant('gate_dup_seeds_and_parameters', lambda a: a['inputs']['dev3'].__setitem__('seeds_sha256', a['inputs']['dev2']['seeds_sha256']), 2, 'same run')
    variant('gate_seeds_shared_with_other_role_ok', lambda a: a['inputs']['c1'].__setitem__('seeds_sha256', a['inputs']['dev2']['seeds_sha256']), 0, 'GATE PASS')
    # review fix 7: program hash
    variant('gate_prog_sha_matching', nop, 0, 'GATE PASS', extra=['--prog-sha', prog])
    variant('gate_prog_sha_mismatch', nop, 2, 'differs from --prog-sha', extra=['--prog-sha', '0' * 64])
    variant('gate_prog_sha_one_run_differs', lambda a: a['inputs']['dev3'].__setitem__('prog_sha256', '1' * 64), 2, 'differs from --prog-sha', extra=['--prog-sha', prog])
    variant('gate_prog_sha_unrecorded', lambda a: a['inputs']['dev3'].__setitem__('prog_sha256', None), 2, 'differs from --prog-sha', extra=['--prog-sha', prog])
    variant('gate_prog_sha_not_a_hash', nop, 2, 'not a sha256', extra=['--prog-sha', 'abc'])
    # review fix 4: a crash is not a FAIL
    variant('gate_crash_args_not_a_dict', lambda a: a['inputs']['dev3'].__setitem__('args', 'x'), 2, 'INPUT_ERROR')
    variant('gate_crash_inputs_a_string', lambda a: a.__setitem__('inputs', 'x'), 2, 'INPUT_ERROR')
    variant('gate_crash_runs_entry_a_list', lambda a: a['runs'].__setitem__('dev3', []), 2, 'INPUT_ERROR')
    variant('gate_crash_pooled_entry_string', lambda a: a['pooled'].__setitem__('dev', 'x'), 2, 'INPUT_ERROR')
    # real crashes (not caught by the gate's own checks: a TypeError / AttributeError / KeyError inside the decision): exit 2 and the exception named
    variant('gate_REAL_crash_n_points_a_string', lambda a: a['inputs']['dev3'].__setitem__('n_points', 'x'), 2, 'TypeError')
    variant('gate_pooled_n_a_list', lambda a: a['pooled']['dev']['S_ens'].__setitem__('n', [1]), 2, 'INPUT_ERROR')
    variant('gate_REAL_crash_runs_entry_n_a_list', lambda a: a['runs']['dev3']['dev']['S_ens'].__setitem__('n', [7]), 2, 'TypeError')
    variant('gate_REAL_crash_bound_stats_None', lambda a: a['runs']['dev3'].__setitem__('dev', None), 2, 'TypeError')
    variant('gate_crash_env_a_list', lambda a: a['inputs']['dev3'].__setitem__('env', []), 2, 'INPUT_ERROR')
    variant('gate_crash_top_level_list', nop, 2, 'INPUT_ERROR', a0=[1, 2])
    variant('gate_crash_refs_with_non_string', lambda a: a.__setitem__('refs', ['cpu', 5]), 2, 'INPUT_ERROR')
    rc, g = py('t207_d5_gate.py', '--bound', str(BOUND), '--expect-points', '7', '--qs', '20', os.path.join(ST, 'does_not_exist.json'))
    record('gate_missing_file', 'rc 2: INPUT_ERROR', rc == 2 and 'INPUT_ERROR' in g and 'GATE FAIL' not in g, g.strip()[-100:])
    with open(os.path.join(ST, 'gate_notjson.json'), 'w') as f:
        f.write('{"measure": "pTT", ')
    rc, g = py('t207_d5_gate.py', '--bound', str(BOUND), '--expect-points', '7', '--qs', '20', os.path.join(ST, 'gate_notjson.json'))
    record('gate_truncated_json', 'rc 2: INPUT_ERROR', rc == 2 and 'INPUT_ERROR' in g and 'GATE FAIL' not in g, g.strip()[-100:])
    pth = os.path.join(ST, 'gate_baseline.json')
    for name, a in (('gate_bound_nan', ['--bound', 'nan', '--expect-points', '7', '--qs', '20']), ('gate_no_expect_points', ['--bound', '0.15', '--qs', '20']),
                    ('gate_no_qs', ['--bound', '0.15', '--expect-points', '7']), ('gate_no_bound', ['--expect-points', '7', '--qs', '20']),
                    ('gate_expect_points_zero', ['--bound', '0.15', '--expect-points', '0', '--qs', '20']),
                    ('gate_qs_not_an_integer', ['--bound', '0.15', '--expect-points', '7', '--qs', 'x']),
                    ('gate_unknown_option', ['--bound', '0.15', '--expect-points', '7', '--qs', '20', '--allow-cover', '40'])):
        rc, g = py('t207_d5_gate.py', *(a + [pth]))
        record(name, 'rc 2: INPUT_ERROR', rc == 2 and 'INPUT_ERROR' in g and 'GATE FAIL' not in g, g.strip().split('\n')[-1][-100:])


# ---------------------------------------------------------------- dump mutations
def dump_layout(path, log):
    with open(path, 'rb') as f:
        h = struct.unpack('<13i', f.read(52))
    n, nv, nvh, dg = h[7], h[9], h[10], h[11]
    base_ = (52 if h[1] == 2 else 48) + 4 * nv         # version 1: 12 header integers, no hasb
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
    NU = '0,4'
    BIG = struct.pack('<d', 1.5e308)          # finite; sums / products of a few of these overflow in the program
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
        # finite device entries that overflow in the program's own arithmetic (error numbers computed from them are not finite)
        ('dev_G_block3_overflows_in_errors', 'def', 'non-finite error number for dev', lambda p, n, o: put(p, o['dev'] + 2 * 16 * n * n, BIG * (2 * n * n))),
        ('dev_G_block1_overflows_in_errors', 'def', 'non-finite error number for dev', lambda p, n, o: put(p, o['dev'], BIG * (2 * n * n))),
        # e finite, only the propagated number B dG Bi overflows (the bracket test below shows that e stays finite at this size): the maxabs_z / fin_r(ep) branch
        ('dev_G_block3_1e305_only_propagated_overflows', 'def', 'non-finite error number for dev', lambda p, n, o: put(p, o['dev'] + 2 * 16 * n * n, struct.pack('<d', 1e305) * (2 * n * n))),
        ('dev_G_block2_overflows_in_propagation', 'def', 'non-finite error number for dev', lambda p, n, o: put(p, o['dev'] + 16 * n * n, BIG * (2 * n * n))),
    ]

    # directories: copies of the dumps of a repro run (v2 dev5 dumps: cover 40, 7 dumps; v1 smoke runs s_def / s_log: cover 1, 16 dumps)
    def mkdumpdir(name, p, mutate=None, k=3, src_tag=None, ndump=7):
        src = base(src_tag or ('m_%s_dev2' % p))
        d = os.path.join(ST, 'dumpruns', 'dump_' + name, 'hub12x12_cb_p')
        os.makedirs(os.path.join(d, 'qref'))
        for f in ('parameters', 'info', 'run.log', 'stamp.txt', 'provenance.json'):
            shutil.copy2(os.path.join(src, f), os.path.join(d, f))
        files = sorted(os.listdir(os.path.join(src, 'qref')))
        assert len(files) == ndump, files
        for f in files:
            shutil.copyfile(os.path.realpath(os.path.join(src, 'qref', f)), os.path.join(d, 'qref', f))
        if mutate:
            fp = os.path.join(d, 'qref', files[k - 1])
            h, n, off = dump_layout(fp, p == 'log')
            mutate(fp, n, off)
        return d, files

    def proc(d, p, keep, qs='0', nul=NU, cover='40', hostin='0', v1=False):
        env = dict(os.environ, D5W=W, D5RUNS=os.path.join(ST, 'dumpruns'), DUMP_COVER=cover, TMPDIR=os.path.join(W, 'tmp'))
        if v1:
            env['ALLOW_V1'] = '1'
        return sh(['bash', os.path.join(HERE, 't207_d5_proc.sh'), d, p, qs, nul, hostin, str(keep)], env=env)

    def sproc(d, p, keep):
        """the cover-1 version 1 smoke runs (stride 1, C = 8, P = 8): qs 1, the null codes of the dev5 pipeline"""
        return proc(d, p, keep, qs='1', nul={'def': '0,4,100004,1,100001,2,1608,4808,5,100005', 'log': '0,4,100004,1,100001,2,5,100005'}[p], cover='1', v1=True)

    def files_left(d):
        return sorted(os.listdir(os.path.join(d, 'qref')))

    # ---- thinning on the version 1 smoke run s_def: KEEP 3 (gcd(3, 8) = 1): dumps 1, 4, 7, 10, 13, 16 stay, 10 go
    d, files = mkdumpdir('thin_baseline', 'def', src_tag='s_def', ndump=16)
    rc, out = sproc(d, 'def', 3)
    left = files_left(d)
    thinned = open(os.path.join(d, 'qref_thinned.txt')).read().split() if os.path.exists(os.path.join(d, 'qref_thinned.txt')) else None
    keepi = [0, 3, 6, 9, 12, 15]
    ok = rc == 0 and left == [files[i] for i in keepi] and thinned == ['qref/' + files[i] for i in range(16) if i not in keepi]
    record('dump_baseline_and_thinning', 'rc 0, 10 of 16 dumps thinned', ok, 'rc %d, left %d dumps, thinned %s' % (rc, len(left), None if thinned is None else len(thinned)))
    chk = open(os.path.join(d, 'qref_check.txt')).read().strip().split('\n')
    record('dump_baseline_check_after_thin', 'check --dumps ok', 'dumps in qref/ match the manifest' in chk[-2], chk[-2][:110])
    record('dump_check_prints_kept_counts', 'per-nst/br/nfe counts', any(l.startswith('QPARSE kept dumps (every 3-th of the run, stride 1): 6 of 16; per nst 1:2 2:1 3:1 4:2; per br') and 'per nfe 1:3 2:3' in l for l in chk), [l for l in chk if 'kept dumps' in l][0][:150])
    dd = d
    rc, out = sproc(dd, 'def', 3)
    record('dump_reprocess_refused', 'rc 1: exists', rc == 1 and 'exists' in out, out.strip().split('\n')[-1][:100])
    # the keep rule at the processing: KEEP 8 (gcd(8, 8) = 8) is refused before anything runs
    d, files = mkdumpdir('keep_rule_refused', 'def', src_tag='s_def', ndump=16)
    rc, out = sproc(d, 'def', 8)
    record('dump_keep_rule_refused_in_proc', 'rc 2, program not run', rc == 2 and 'keep rule' in out and not os.path.exists(os.path.join(d, 'qref.txt')) and len(files_left(d)) == 16, out.strip().split('\n')[-1][:110])
    d, files = mkdumpdir('keep_rule_even_cover', 'def')
    rc, out = proc(d, 'def', 1)                      # cover 40, keep 1: every 40th, the biased subset of the first review
    record('dump_keep_rule_refused_even_cover', 'rc 2, program not run', rc == 2 and 'keep rule' in out and not os.path.exists(os.path.join(d, 'qref.txt')), out.strip().split('\n')[-1][:110])
    # a failed check: thin must refuse by itself, also when called directly (review fix 6); NaN in the right state of dump 3 of s_def
    d, files = mkdumpdir('thin_failed_check', 'def', lambda p, n, o: put(p, o['U1'] + 8, NAN), src_tag='s_def', ndump=16)
    rc, out = sproc(d, 'def', 3)
    left = files_left(d)
    record('dump_failed_check_nothing_deleted_by_proc', 'rc 3, 16 dumps left', rc == 3 and len(left) == 16, 'rc %d, left %d: %s' % (rc, len(left), out.strip().split('\n')[-1][:80]))
    rc, out = py('t207_d5_qparse.py', 'thin', d, '--allow-v1')
    record('thin_refuses_on_failed_check_directly', 'rc 1, nothing deleted', rc == 1 and len(files_left(d)) == 16 and not os.path.exists(os.path.join(d, 'qref_thinned.txt')), last(out))
    # thin's own preconditions: a complete passing directory by hand (qref.txt, qref_proc.json of the good s_def processing), one defect at a time
    gd = os.path.join(ST, 'dumpruns', 'dump_thin_baseline', 'hub12x12_cb_p')

    def thin_dir(name):
        d, files = mkdumpdir(name, 'def', src_tag='s_def', ndump=16)
        rc, out = py('t207_d5_qparse.py', 'spec', d, '--path', 'def', '--qs', '1', '--nul', '0,4,100004,1,100001,2,1608,4808,5,100005', '--hostin', '0', '--cover', '1', '--keep', '3')
        assert rc == 0, out
        sh(['bash', '-c', 'cd %s && sha256sum qref/q*.bin > qref_manifest.txt' % d])
        shutil.copy2(os.path.join(gd, 'qref.txt'), os.path.join(d, 'qref.txt'))
        return d, files

    d, files = thin_dir('thin_no_proc_json')
    rc, out = py('t207_d5_qparse.py', 'thin', d, '--allow-v1')
    record('thin_refuses_without_qref_proc_json', 'rc 1, nothing deleted', rc == 1 and 'qref_proc.json missing' in out and len(files_left(d)) == 16, last(out))
    pj = json.load(open(os.path.join(gd, 'qref_proc.json')))
    for name, edit, substr in (('check_rc_1', lambda x: x.__setitem__('check_rc', 1), 'not both 0'), ('rc_3', lambda x: x.__setitem__('rc', 3), 'not both 0'),
                               ('keep_differs', lambda x: x['args'].__setitem__('keep', 5), 'differ')):
        d, files = thin_dir('thin_proc_json_' + name)
        x = json.loads(json.dumps(pj)); edit(x)
        json.dump(x, open(os.path.join(d, 'qref_proc.json'), 'w'))
        rc, out = py('t207_d5_qparse.py', 'thin', d, '--allow-v1')
        record('thin_refuses_proc_json_' + name, 'rc 1: %s, nothing deleted' % substr, rc == 1 and substr in out and len(files_left(d)) == 16, last(out))
    d, files = thin_dir('thin_qref_txt_bad')
    shutil.copy2(os.path.join(gd, 'qref_proc.json'), os.path.join(d, 'qref_proc.json'))
    t = open(os.path.join(d, 'qref.txt')).read()
    open(os.path.join(d, 'qref.txt'), 'w').write(t.replace('QREF file 2 ', 'QREF file 3 ', 1))
    rc, out = py('t207_d5_qparse.py', 'thin', d, '--allow-v1')
    record('thin_refuses_bad_qref_txt', 'rc 1, nothing deleted', rc == 1 and 'QPARSE ERROR' in out and len(files_left(d)) == 16, last(out))
    d, files = thin_dir('thin_sha_guard')
    shutil.copy2(os.path.join(gd, 'qref_proc.json'), os.path.join(d, 'qref_proc.json'))
    with open(os.path.join(d, 'qref', files[1]), 'ab') as f:
        f.write(b'x')
    rc, out = py('t207_d5_qparse.py', 'thin', d, '--allow-v1')
    record('thin_guard_sha_mismatch', 'rc 1, nothing deleted', rc == 1 and 'sha256 differs' in out and len(files_left(d)) == 16, last(out))
    d, files = thin_dir('thin_good_by_hand')
    shutil.copy2(os.path.join(gd, 'qref_proc.json'), os.path.join(d, 'qref_proc.json'))
    rc, out = py('t207_d5_qparse.py', 'thin', d, '--allow-v1')
    record('thin_by_hand_good_directory_thins', 'rc 0, 6 dumps kept', rc == 0 and len(files_left(d)) == 6 and 'kept dumps' in out, last(out))
    d, files = thin_dir('thin_v1_without_flag')
    shutil.copy2(os.path.join(gd, 'qref_proc.json'), os.path.join(d, 'qref_proc.json'))
    rc, out = py('t207_d5_qparse.py', 'thin', d)
    record('thin_v1_dumps_need_allow_v1', 'rc 1, nothing deleted', rc == 1 and 'allow_v1' in out and len(files_left(d)) == 16, last(out))
    # symbolic links are not thinned (s_def: all links)
    d, files = mkdumpdir('symlinks', 'def', src_tag='s_def', ndump=16)
    for f in files:
        fp = os.path.join(d, 'qref', f)
        tgt = os.path.join(ST, 'dumpruns', 'dump_symlinks', 'orig_' + f)
        shutil.move(fp, tgt)
        os.symlink(tgt, fp)
    rc, out = sproc(d, 'def', 3)
    left = files_left(d)
    record('dump_symlinks_not_thinned', 'rc 4, nothing deleted', rc == 4 and len(left) == 16 and 'symbolic link' in out + open(os.path.join(d, 'qref_check.txt')).read(),
           'rc %d, left %d: %s' % (rc, len(left), out.strip().split('\n')[-1][:90]))
    # a run directory outside the designated tree is refused
    d, files = mkdumpdir('outside_runs_tree', 'def')
    rc, out = sh(['bash', os.path.join(HERE, 't207_d5_proc.sh'), d, 'def', '0', '0,4', '0', '0'], env=dict(os.environ, D5W=W, D5RUNS=os.path.join(W, 'runs'), DUMP_COVER='40'))
    record('dump_proc_refuses_dir_outside_runs', 'rc 1: refusing', rc == 1 and 'refusing' in out and not os.path.exists(os.path.join(d, 'qref_spec.json')), out.strip().split('\n')[-1][:100])
    # program hash (review fix 7): a binary that differs from its .sha256 file is refused before anything is written
    wd = os.path.join(ST, 'badbin')
    os.makedirs(os.path.join(wd, 'bin'))
    os.makedirs(os.path.join(wd, 'tmp'))
    shutil.copy2(os.path.join(W, 'bin', 'qref_pp.out'), os.path.join(wd, 'bin', 'qref_pp.out'))
    shutil.copy2(os.path.join(W, 'bin', 'qref_pp.out.sha256'), os.path.join(wd, 'bin', 'qref_pp.out.sha256'))
    with open(os.path.join(wd, 'bin', 'qref_pp.out'), 'ab') as f:
        f.write(b'\0')
    d, files = mkdumpdir('badbin', 'def')
    env = dict(os.environ, D5W=wd, D5RUNS=os.path.join(ST, 'dumpruns'), DUMP_COVER='40')
    rc, out = sh(['bash', os.path.join(HERE, 't207_d5_proc.sh'), d, 'def', '0', NU, '0', '0'], env=env)
    record('proc_refuses_binary_not_matching_sha256_file', 'rc 1: refusing', rc == 1 and 'sha256 of' in out and not os.path.exists(os.path.join(d, 'qref_spec.json')), out.strip().split('\n')[-1][:110])
    # the dump mutations (v2 dev5 dumps, cover 40: KEEP 0, the biased thinning is refused)
    for name, p, substr, fn in muts:
        d, files = mkdumpdir(name, p, fn)
        rc, out = proc(d, p, 0)
        txt = open(os.path.join(d, 'qref.txt')).read() if os.path.exists(os.path.join(d, 'qref.txt')) else ''
        err = [l for l in txt.split('\n') if l.startswith('QREFERR')]
        endl = [l for l in txt.split('\n') if l.startswith('QREFEND')]
        chk = open(os.path.join(d, 'qref_check.txt')).read() if os.path.exists(os.path.join(d, 'qref_check.txt')) else ''
        prc = json.load(open(os.path.join(d, 'qref_proc.json'))) if os.path.exists(os.path.join(d, 'qref_proc.json')) else {}
        left = os.listdir(os.path.join(d, 'qref'))
        ok = (rc == 3 and len(err) == 1 and files[2] in err[0] and substr in err[0] and prc.get('rc') == 3 and prc.get('check_rc') == 1
              and endl == ['QREFEND dumps 7 ok 6 err 1'] and len(left) == 7 and 'QPARSE ERROR' in chk and 'the program reported an error' in chk)
        record('dump_' + name, 'QREFERR, rc 3, nothing deleted', ok, 'rc %d; %s; left %d' % (rc, (err[0] if err else out.strip().split('\n')[-1])[-100:], len(left)))
    # bracket for the previous test: the same block at 1e300 is finite everywhere (accepted as data, a huge error); its error numbers e, times 1e5,
    # are still finite, so at 1e305 only the propagated number can be what overflows
    d, files = mkdumpdir('dev_G_block3_1e300_finite', 'def', lambda p, n, o: put(p, o['dev'] + 2 * 16 * n * n, struct.pack('<d', 1e300) * (2 * n * n)))
    rc, out = proc(d, 'def', 0)
    txt = open(os.path.join(d, 'qref.txt')).read() if os.path.exists(os.path.join(d, 'qref.txt')) else ''
    k3 = int(files[2][1:7])
    l3 = [l for l in txt.split('\n') if l.startswith('QREF file %d ' % k3)]
    ev = None
    if l3:
        dv = [x for x in l3[0].split(' | ') if x.startswith('dev ')][0].split()[1:]
        ev = [float(x) for x in dv]
    ok = rc == 0 and ev is not None and max(ev[:5]) * 1e5 < 1.7e308 and ev[5] > 1e298 and 'QREFERR' not in txt
    record('dump_dev_G_block3_1e300_finite_bracket', 'rc 0, e * 1e5 finite', ok, 'rc %d, dev numbers %s' % (rc, None if ev is None else ' '.join('%.3e' % x for x in ev)))
    # a version 2 dump without B, Bi (hasb cleared, B and Bi cut off): the program must print the propagated numbers as -1 and not fail; the check
    # refuses it because the run's parameters say that this point (nst 1 < NSTM) has B
    d, files = mkdumpdir('B_dropped', 'def', lambda p, n, o: (os.truncate(p, o['end'] - 32 * n * n), put(p, 48, struct.pack('<i', 0))))
    rc, out = proc(d, 'def', 0)
    txt = open(os.path.join(d, 'qref.txt')).read()
    prc = json.load(open(os.path.join(d, 'qref_proc.json')))
    chk = open(os.path.join(d, 'qref_check.txt')).read()
    l3 = [l for l in txt.split('\n') if l.startswith('QREF file %d ' % int(files[2][1:7]))]
    ok = rc == 3 and prc['rc'] == 0 and prc['check_rc'] == 1 and 'header hasb = 0, expected 1' in chk and len(l3) == 1 and l3[0].count('-1.0000E+000') == 3 * 8 and len(os.listdir(os.path.join(d, 'qref'))) == 7
    record('dump_B_dropped_hasb_cleared', 'program rc 0, check error, kept', ok, 'rc %d, program rc %s, %s' % (rc, prc['rc'], chk.strip()[-80:]))
    # directory checks
    d, files = mkdumpdir('extra_file', 'def')
    shutil.copyfile(os.path.join(d, 'qref', files[0]), os.path.join(d, 'qref', 'q000002.bin'))
    rc, out = proc(d, 'def', 0)
    record('dump_extra_file_in_qref', 'rc 2: unexpected', rc == 2 and 'unexpected' in out and not os.path.exists(os.path.join(d, 'qref.txt')), out.strip().split('\n')[-1][:100])
    d, files = mkdumpdir('missing_file', 'def')
    os.remove(os.path.join(d, 'qref', files[4]))
    rc, out = proc(d, 'def', 0)
    record('dump_missing_file_in_qref', 'rc 2: missing', rc == 2 and 'missing' in out, out.strip().split('\n')[-1][:100])
    d, files = mkdumpdir('wrong_cover', 'def')
    rc, out = proc(d, 'def', 0, cover='20')
    record('dump_wrong_cover_argument', 'rc 2', rc == 2 and 'missing' in out, out.strip().split('\n')[-1][:100])
    d, files = mkdumpdir('wrong_path', 'def')
    rc, out = proc(d, 'log', 0)
    record('dump_wrong_path_argument', 'rc 2: run is a def run', rc == 2 and 'is a def run' in out, out.strip().split('\n')[-1][:100])
    d, files = mkdumpdir('stale_manifest', 'def')
    sh(['bash', '-c', 'cd %s && sha256sum qref/q*.bin > qref_manifest.txt && echo 0000000000000000000000000000000000000000000000000000000000000000  qref/q000001.bin >> qref_manifest.txt' % d])
    rc, out = proc(d, 'def', 0)
    record('dump_stale_manifest', 'rc 2: manifest differs', rc == 2 and 'manifest' in out, out.strip().split('\n')[-1][:100])
    # the LOG path baseline (one run through the program with the LOG build); its proc record carries the program hash
    d, files = mkdumpdir('log_baseline', 'log')
    rc, out = proc(d, 'log', 0)
    pr = json.load(open(os.path.join(d, 'qref_proc.json')))
    want = open(os.path.join(W, 'bin', 'qref_log.out.sha256')).read().split()[0]
    record('dump_log_baseline', 'rc 0, prog hash recorded', rc == 0 and pr['prog'] == want, 'rc %d, prog %s' % (rc, pr['prog'][:16]))


def main():
    if os.path.exists(ST):
        sys.exit('refusing: %s exists' % ST)
    os.makedirs(ST)
    if WHAT in ('text', 'all'):
        gate_tests()
        text_tests()
        agg_tests()
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
