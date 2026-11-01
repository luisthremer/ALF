#!/usr/bin/env python3
"""Ticket 207-2: mutation self-tests of the QREF analysis pipeline (nividia, CPU only).   usage: t207_d5_selftest.py WORKDIR [text|gate|reproc|dump|all]

Needs WORKDIR/runs/{m_def_dev1..4,c_def_dev1,d_def_dev1,m_log_dev1,s_def,s_log}/hub12x12_cb_p produced by t207_d5_repro.sh (the new pipeline on the kept dev5
dumps) and WORKDIR/bin/qref_pp.out, qref_log.out.  Everything it writes goes to WORKDIR/st/ (created; refuses to start if it exists).
  text       mutations of qref.txt, of the run's records and of the spec: each must end in the error named (the strict parser, agg, stat)
  gate       the gate (design 3) on copies of the run directories m_def_dev2..4 + c_def_dev1 with a gate spec: the baseline, every spec
             defect, every evidence defect (one per copy; records edited consistently where the test is about a check behind them) and every
             planted defect of the numbers: INPUT_ERROR / FAIL / DEV-ONLY as named; every gate run must leave a report with the same verdict
             and exit status.  The kept dev5 dumps are a subset (cover 40), which the gate refuses with a cover-1 spec (tested); the gate logic
             is tested through a spec marked test_only, whose all-criteria-met verdict is TEST-ONLY (exit 4), never PASS.  A genuine cover-1
             stride-3 PASS needs 207-3's runs.  The baseline's numbers are compared with agg's CLI on the same files (1e-12).
  dump       mutations of real dumps (NaN, Infinity, zeroed state, truncation, extra bytes, bad header) run through t207_d5_proc.sh
             (program + strict check): the program must print a QREFERR line for the dump, exit with status 3, the check must fail and no
             dump may be deleted; plus the directory checks (extra / missing dump, stale manifest, symlink and sha256 guards of the thinning).
             The dev5 v2 dumps are kept at an even spacing (cover 20/40), which the keep rule refuses for thinning, so the thinning is tested
             on the cover-1 version 1 smoke runs s_def / s_log (stride 1, KEEP 3); the v2 mutations run with KEEP 0.
Output: one table line per test, generated from the actual outcomes, with the exit status of the last command the test ran (also in
st/results.txt, column last_rc); exit status 1 if any test did not give the expected outcome.
The bound used is 0.15 (the proposed criterion of 207 dev5; the self-test is not a calibration and does not change any bound).
Fix round 4 (review 207-2d): the gate takes --spec-sha and --prog and re-runs the program on the kept dumps.  Most gate tests call the
gate's main() in-process (the same function the CLI runs; the gate keeps re-processing results of identical program / dump bytes / arguments
within one process, so the spec variations do not re-run the program); the baselines and the reviewer's probes run the CLI.  Tests whose
planted defect is an edit of qref.txt numbers (the decision logic) run in-process with the gate's authenticate step replaced by the self-test
(marked 'reprocessing replaced' in the detail; the report records it); each of them has a companion with the real step, which must
end in INPUT_ERROR because the edited numbers are not what the program computes.
  reproc     the gate's re-processing check on every kept dev5 dump of every repro run (every number exactly equal)
Fix round 5 (review 207-2e): the spec pins keep, hostin and min_reprocessed (M1); the re-processed numbers must be exactly equal (N1: a
one-digit difference is now INPUT_ERROR, test F1g changed); end-to-end device FAIL through the real processing (N2, dump mode, test
N2_*); control member that is a reference (N3); minimal program env, LD_PRELOAD refused (N4); test S4 renamed (N7).  Manager items: X1 the host member label from the spec (no hard-coded label in gate / qparse / agg; tests X1_*),
X2 t207_d5_proc.sh sources ~/alf-gpu/env/gnu.sh only if it exists.
207-5 (harness, round-5 review follow-ups): Y1 ctl_dev is a reserved host label; Y2 a kept dump changed while the program re-processes it is
refused (re-hash at read time; Y2b the same copy unchanged keeps the baseline verdict); Z0 t207_d5_pairs.py with agg.ref_list's host argument
(text group); Z1-Z3 t207_d5_proc.sh: interpreter D5_PYTHON (default ~/alf-gpu/venv/bin/python3) refused if not executable, used and recorded in
qref_proc.json with the machine record (dump group).  Localized fault (manager, 2026-11-01): ALF_GPU_FAULT_TAU_AT=<nst>:<nfe> only with
ALF_GPU_FAULT_TAU (E14-E17); changed pre-existing assertions: E13_1/E13_2 ('eps:2.5e-10@2', 'eps:1e-9@3:2' valid) moved to E12_9/E12_10 (now
refused: the '@' form is dropped), so E13 has two cases.
"""
import contextlib, hashlib, importlib.util, io, json, math, os, re, shutil, struct, subprocess, sys
sys.dont_write_bytecode = True

W = os.path.abspath(sys.argv[1])
WHAT = sys.argv[2] if len(sys.argv) > 2 else 'all'
HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
_sp = importlib.util.spec_from_file_location('qp', os.path.join(HERE, 't207_d5_qparse.py'))
QP = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(QP)
_gs = importlib.util.spec_from_file_location('t207_d5_gate_inproc', os.path.join(HERE, 't207_d5_gate.py'))
GATE = importlib.util.module_from_spec(_gs)
_gs.loader.exec_module(GATE)
BOUND = 0.15
ST = os.path.join(W, 'st')
RESULTS = []


LAST = {'rc': None}


def sh(argv, env=None, timeout=3600):
    p = subprocess.run(argv, capture_output=True, text=True, env=env, timeout=timeout)
    LAST['rc'] = p.returncode
    return p.returncode, (p.stdout + p.stderr)


def py(script, *args):
    return sh([PY, '-I', os.path.join(HERE, script)] + list(args))


def record(name, expect, ok, detail):
    """one result line; the exit status of the last command the test ran is recorded with it (results.txt column 'last_rc')"""
    RESULTS.append((name, expect, ok, detail, LAST['rc']))
    print('%-34s %-34s %-5s rc %-4s %s' % (name, expect, 'ok' if ok else 'WRONG', LAST['rc'], detail.replace('\n', ' ')[:150]), flush=True)
    LAST['rc'] = None


def base(tag):
    return os.path.join(W, 'runs', tag, 'hub12x12_cb_p')


REC = ('parameters', 'info', 'run.log', 'stamp.txt', 'provenance.json', 'qref_spec.json', 'qref_proc.json')


def sha_file(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def fix_qsha(d):
    """set qref_sha256 of RUNDIR/qref_proc.json to the hash of its qref.txt as it is now (a consistent edit, used to reach the checks behind it)"""
    pj = os.path.join(d, 'qref_proc.json')
    if os.path.exists(pj):
        pv = json.load(open(pj))
        pv['qref_sha256'] = sha_file(os.path.join(d, 'qref.txt'))
        json.dump(pv, open(pj, 'w'), indent=1)


def fix_msha(d):
    """set manifest_sha256 of RUNDIR/qref_proc.json to the hash of its qref_manifest.txt as it is now (consistent edit, fix round 4)"""
    pj = os.path.join(d, 'qref_proc.json')
    pv = json.load(open(pj))
    pv['manifest_sha256'] = sha_file(os.path.join(d, 'qref_manifest.txt'))
    json.dump(pv, open(pj, 'w'), indent=1)


def replace_dump(d, name, data):
    """qref/NAME (often a symbolic link to a dev5 dump) replaced by a regular file holding data; manifest line and manifest_sha256 follow"""
    p = os.path.join(d, 'qref', name)
    os.remove(p)
    with open(p, 'wb') as f:
        f.write(data)
    mp = os.path.join(d, 'qref_manifest.txt')
    rewrite(mp, lambda t: re.sub(r'^[0-9a-f]{64}(  qref/%s)$' % re.escape(name), lambda m: sha_file(p) + m.group(1), t, flags=re.M))
    fix_msha(d)


def dump_bytes(d, name):
    with open(os.path.realpath(os.path.join(d, 'qref', name)), 'rb') as f:
        return f.read()


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
        fix_qsha(d)                   # the processing record names the mutated qref.txt (consistent edit): the inner checks are reached
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


def prog_of_runs():
    return json.load(open(os.path.join(base('m_def_dev2'), 'qref_proc.json')))['prog']


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


def pairs_tests():
    """207-5: t207_d5_pairs.py runs with agg.ref_list's signature (host label argument, default nvh) on the repro runs m_def_dev1..4"""
    p = subprocess.run([PY, '-I', os.path.join(HERE, 't207_d5_pairs.py'), '0.15', 'def'], capture_output=True, text=True, cwd=W, timeout=3600)
    LAST['rc'] = p.returncode
    rows = [l for l in p.stdout.split('\n') if re.match(r'^  dev[1-4]\+dev[1-4] ', l)]
    record('Z0_pairs_default_host_nvh', 'rc 0, 6 pairs', p.returncode == 0 and len(rows) == 6 and 'path def, bound 0.15' in p.stdout,
           (rows[0] if rows else (p.stdout + p.stderr).strip()[-150:]))
    p = subprocess.run([PY, '-I', os.path.join(HERE, 't207_d5_pairs.py'), '0.15', 'def', 'cpu'], capture_output=True, text=True, cwd=W, timeout=3600)
    LAST['rc'] = p.returncode
    record('Z0b_pairs_reserved_host_label', 'rc != 0: host member label', p.returncode != 0 and 'host member label' in p.stdout + p.stderr,
           (p.stdout + p.stderr).strip()[-150:])


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
        # the by-construction bounds min(e1..e4)/(2n) <= eF <= 2n max(e1..e4) (t207_d5_qparse.check_ef)
        ('eF_above_2n_max_block', lambda t: edit_line(t, 1, set_tok('dev', 4, '1.0000E-005')), 'outside [min(e1..e4)/(2n), 2n max(e1..e4)]'),
        ('eF_below_min_block_over_2n', lambda t: edit_line(t, 1, set_tok('dev', 4, '1.0000E-020')), 'outside [min(e1..e4)/(2n), 2n max(e1..e4)]'),    ]
    for name, mut, substr in tests:
        q = mkrun(name, 'm_def_dev1', mut(t0))
        expect_error(name, q, substr, absent='unexpected line' if name == 'QREFERR_line' else None)
    # zero error (well-formed): the parser accepts it, agg refuses to take a log
    q = mkrun('zero_error_dev_pTT', 'm_def_dev1', edit_line(t0, 1, set_tok('dev', 5, '0.0000E+000')))
    expect_error('zero_error_dev_pTT', q, 'not positive and finite')
    # well-formed huge value (the old clamp value): accepted as data (a very large error, no longer a placeholder): the pooled S_ens moves
    q = mkrun('overflow_value_dev_pTT', 'm_def_dev1', edit_line(t0, 1, set_tok('dev', 5, '1.7977E+308')))
    expect_error('overflow_value_dev_pTT', q, 'is not finite')
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
    q2 = mkrun('prog_no_proc_json_no_prog_sha_arg', 'm_def_dev2')           # fix round 2: the processing record is required whatever the arguments
    os.remove(os.path.join(os.path.dirname(q2), 'qref_proc.json'))
    rc, msg, out = agg_direct('x_prog3b', 'def', [('main', 'a', q2)])
    record('agg_proc_json_missing_without_prog_sha', 'error: is missing', rc == 1 and 'is missing' in msg and not os.path.exists(out), last(msg))
    q2 = mkrun('proc_json_prog_invalid', 'm_def_dev2')
    pj = os.path.join(os.path.dirname(q2), 'qref_proc.json')
    d = json.load(open(pj)); d['prog'] = 'xyz'; json.dump(d, open(pj, 'w'))
    rc, msg, out = agg_direct('x_prog3c', 'def', [('main', 'a', q2)])
    record('agg_proc_json_prog_not_a_hash', 'error: no valid program sha256', rc == 1 and 'no valid program sha256' in msg, last(msg))
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


# ---------------------------------------------------------------- gate (design 3: recomputed from the run directories, no aggregate input)
GSET = ['m_def_dev2', 'm_def_dev3', 'm_def_dev4', 'c_def_dev1_q20']      # three main seed sets with 7 points with B each, one TAU_NS=0 control run processed with qs 20 like the main runs (repro tag c_def_dev1_q20)
NUL_DEF = [0, 4, 100004, 1, 100001, 2, 1608, 4808, 5, 100005]
NUL_LOG = [0, 4, 100004, 1, 100001, 2, 5, 100005]


def mlist(nul, var=True, host='nvh'):
    """the ordered statistical member list of a run (main: with var8; control runs of dev5: without QVAR)"""
    return ['dev'] + (['var8'] if var else []) + [host, 'cpu'] + ['n%d' % c for c in nul] + ['u1', 'u2']
GW = {0: 'GATE PASS', 1: 'GATE FAIL', 2: 'GATE INPUT_ERROR', 3: 'GATE DEV-ONLY', 4: 'GATE TEST-ONLY'}


def gdir():
    return os.path.join(ST, 'g')


def gspec(**kw):
    """the gate spec for the kept dev5 subset (every 40th dump): marked test_only, since the subset has cover 40.  Bound 0.15 = the proposed
    criterion of 207 dev5 (the self-test is not a calibration).  A genuine cover-1 stride-3 run set does not exist yet (207-3)."""
    s = dict(spec_version=1, path='def', host_member='nvh', members={'main': mlist(NUL_DEF), 'ctl_tau_ns0': mlist(NUL_DEF, var=False)},
             refs=QP.ref_names('def', 'nvh'), qs=20, nul=list(NUL_DEF), cover=40, keep=0, hostin=0, min_reprocessed=7,
             test_only='207-2 self-test on the kept dev5 subset (cover 40, 7 points with B per run): tests the gate logic only; '
                       'a genuine cover-1 stride-3 run set is pending (207-3)',
             expect_points={'pTT': 7}, prog_sha=prog_of_runs(), src_sha=sha_file(os.path.join(HERE, 't207_d5_qref.F90')),
             binary_sha256=prov_of('m_def_dev2', 'binary_sha256'), parameters_sha256=prov_of('m_def_dev2', 'parameters_sha256'),
             seeds={'main': [seeds_of(t) for t in GSET[:3]], 'ctl_tau_ns0': [seeds_of('c_def_dev1_q20')]},
             env_roles={'main': {'ALF_GPU_DEBUG': '1', 'ALF_GPU_TAUDIAG_QHOST': '1', 'ALF_GPU_TAUDIAG_QVAR': '8',
                                 'ALF_GPU_TAUDIAG_QREF': '2', 'ALF_GPU_TAUDIAG_FROM': '331'},
                        'ctl_tau_ns0': {'extra_env': {'ALF_GPU_TAU_NS': '0'}, 'omit_env': ['ALF_GPU_TAUDIAG_QVAR']}},
             n_runs={'main': 3, 'ctl_tau_ns0': 1},
             gates=[dict(name='acc_pTT', measure='pTT', statistic='S_ens', pool='pooled', bound=BOUND)],
             controls=[dict(name='var8', kind='variant', member='var8', trips='acc_pTT'),
                       dict(name='tau_ns0', kind='run', role='ctl_tau_ns0', trips='acc_pTT')])
    s.update(kw)
    return s


def smut(fn, **kw):
    s = gspec(**kw)
    fn(s)
    return s


def gcopy(name, tags=None, edits=None, src=None):
    """copies of run directories (records, qref.txt, processing records, qref/ with its symbolic links) under st/g/NAME/LABEL; tags: the
    labels (default GSET), src: {label: source tag} for labels that are not tags; edits: {label: fn(dir)}; returns the run directories"""
    out = []
    tags = tags or GSET
    for lab in (edits or {}):
        assert lab in tags, lab
    for lab in tags:
        d = os.path.join(gdir(), name, lab, 'hub12x12_cb_p')
        shutil.copytree(base((src or {}).get(lab, lab)), d, symlinks=True)
        if edits and lab in edits:
            edits[lab](d)
        out.append(d)
    return out


def fix_records(d):
    """after an edit of records: provenance.json hashes, the parameters prefix of stamp.txt and the stored spec's record hashes follow the
    edit (a consistent forgery, so that the checks behind these are reached)"""
    pj = os.path.join(d, 'provenance.json')
    pv = json.load(open(pj))
    for key, f in (('parameters_sha256', 'parameters'), ('run_log_sha256', 'run.log'), ('info_sha256', 'info')):
        pv[key] = sha_file(os.path.join(d, f))
    json.dump(pv, open(pj, 'w'), indent=1)
    st_ = open(os.path.join(d, 'stamp.txt')).read()
    open(os.path.join(d, 'stamp.txt'), 'w').write(re.sub(r'( parameters )[0-9a-f]{16}', lambda m: m.group(1) + pv['parameters_sha256'][:16], st_))
    sp = os.path.join(d, 'qref_spec.json')
    s = json.load(open(sp))
    s['run']['records_sha256'] = {r: sha_file(os.path.join(d, r)) for r in QP.RECORDS}
    s['run']['parameters_sha256'] = pv['parameters_sha256']
    s['run']['seeds_sha256'] = pv['seeds_sha256']
    s['run']['binary_sha256'] = pv['binary_sha256']
    json.dump(s, open(sp, 'w'), indent=1)


def e_qref(fn):
    """edit of qref.txt (text -> text); qref_sha256 of the processing record follows (consistent)"""
    def e(d):
        p = os.path.join(d, 'qref.txt')
        t = open(p).read()
        t2 = fn(t)
        assert t2 != t
        open(p, 'w').write(t2)
        fix_qsha(d)
    return e


def e_qref_raw(fn):
    """edit of qref.txt without updating the processing record"""
    def e(d):
        p = os.path.join(d, 'qref.txt')
        t = open(p).read()
        t2 = fn(t)
        assert t2 != t
        open(p, 'w').write(t2)
    return e


def e_proc(fn):
    def e(d):
        p = os.path.join(d, 'qref_proc.json')
        x = json.load(open(p))
        fn(x)
        json.dump(x, open(p, 'w'), indent=1)
    return e


def e_args(**kv):
    """the processing arguments changed consistently in qref_spec.json and qref_proc.json"""
    def e(d):
        for f in ('qref_spec.json', 'qref_proc.json'):
            p = os.path.join(d, f)
            x = json.load(open(p))
            for k, v in kv.items():
                x['args'][k] = (','.join(str(c) for c in v) if (f == 'qref_proc.json' and k == 'nul') else v)
            json.dump(x, open(p, 'w'), indent=1)
    return e


def e_env(fn):
    """the run's ALF_GPU_* env changed by fn(dict) consistently in info, stamp.txt, provenance.json and the stored spec"""
    def e(d):
        sp = json.load(open(os.path.join(d, 'qref_spec.json')))
        new = fn(dict(sp['run']['env']))
        t = open(os.path.join(d, 'info')).read()
        t2 = re.sub(r'^( *GPU debug env: ).*$', lambda m: m.group(1) + ', '.join('%s=%s' % kv for kv in sorted(new.items())), t, flags=re.M)
        assert t2 != t
        open(os.path.join(d, 'info'), 'w').write(t2)
        t = open(os.path.join(d, 'stamp.txt')).read()
        open(os.path.join(d, 'stamp.txt'), 'w').write(re.sub(r'^env .*$', lambda m: 'env ' + ' '.join('%s=%s' % kv for kv in new.items()), t, flags=re.M))
        pv = json.load(open(os.path.join(d, 'provenance.json')))
        pv['env'] = ['%s=%s' % kv for kv in new.items()]
        json.dump(pv, open(os.path.join(d, 'provenance.json'), 'w'), indent=1)
        sp['run']['env'] = new
        json.dump(sp, open(os.path.join(d, 'qref_spec.json'), 'w'), indent=1)
        fix_records(d)
    return e


def e_seeds(full):
    """the run's seeds hash set to FULL (64 hex) consistently in stamp.txt, provenance.json and the stored spec"""
    def e(d):
        t = open(os.path.join(d, 'stamp.txt')).read()
        open(os.path.join(d, 'stamp.txt'), 'w').write(re.sub(r'^seeds [0-9a-f]{16}', 'seeds ' + full[:16], t, flags=re.M))
        pv = json.load(open(os.path.join(d, 'provenance.json')))
        pv['seeds_sha256'] = full
        json.dump(pv, open(os.path.join(d, 'provenance.json'), 'w'), indent=1)
        fix_records(d)
    return e


def rewrite(p, fn):
    """file P replaced by fn(its text); the text is read before the file is opened for writing"""
    t = open(p).read()
    t2 = fn(t)
    assert t2 != t, 'rewrite of %s had no effect' % p
    with open(p, 'w') as f:
        f.write(t2)


def e_file(fn):
    def e(d):
        fn(d)
    return e


def chain(*es):
    def e(d):
        for x in es:
            x(d)
    return e


def fmt(x):
    """a number in the program's ES12.4E3 format"""
    m, ex = ('%.4E' % x).split('E')
    return '%sE%s%03d' % (m, ex[0], abs(int(ex)))


def lines_all(fn, only=None):
    """fn applied to every QREF data line (or to the data line number only, 1-based)"""
    def g(t):
        ls = t.split('\n')
        k = 0
        for i, l in enumerate(ls):
            if l.startswith('QREF file '):
                k += 1
                if only is None or k == only:
                    ls[i] = fn(l)
        return '\n'.join(ls)
    return g


def scale(member, idx, f):
    def g(line):
        segs = line.split(' | ')
        for j, s in enumerate(segs):
            t = s.split()
            if j and t[0] == member:
                t[1 + idx] = fmt(float(t[1 + idx]) * f)
                segs[j] = ' '.join(t)
        return ' | '.join(segs)
    return g


def setnum(member, idx, tok):
    def g(line):
        segs = line.split(' | ')
        for j, s in enumerate(segs):
            t = s.split()
            if j and t[0] == member:
                t[1 + idx] = tok
                segs[j] = ' '.join(t)
        return ' | '.join(segs)
    return g


def copy_member(dst, src):
    def g(line):
        segs = line.split(' | ')
        nums = {s.split()[0]: s.split()[1:] for s in segs[1:]}
        for j, s in enumerate(segs):
            if j and s.split()[0] == dst:
                segs[j] = dst + ' ' + ' '.join(nums[src])
        return ' | '.join(segs)
    return g


def swap_members(a, b):
    def g(line):
        segs = line.split(' | ')
        nums = {s.split()[0]: s.split()[1:] for s in segs[1:]}
        for j, s in enumerate(segs):
            n = s.split()[0]
            if j and n in (a, b):
                segs[j] = n + ' ' + ' '.join(nums[b if n == a else a])
        return ' | '.join(segs)
    return g


def drop_member(m):
    return lambda line: ' | '.join(s for s in line.split(' | ') if not s.startswith(m + ' '))


def qbin(path):
    """the offline program of a path (the binaries the repro runs were processed with)"""
    return os.path.join(W, 'bin', 'qref_pp.out' if path == 'def' else 'qref_log.out')


NOAUTH_TAG = 'SKIPPED: reprocessing replaced by the self-test (planted numbers in qref.txt; decision-logic test)'


def gate_inproc(argv, noauth=False):
    """GATE.main(argv) in this process (the function the CLI runs), stdout captured; noauth: the authenticate step replaced (recorded in the
    report as NOAUTH_TAG), restored afterwards"""
    buf = io.StringIO()
    saved = GATE.authenticate
    if noauth:
        GATE.authenticate = lambda runs, s, prog: [dict(dir=x['dir'], dumps=0, numbers=0, equal=0, source=NOAUTH_TAG) for x in runs]
    try:
        with contextlib.redirect_stdout(buf):
            rc = GATE.main(list(argv))
    finally:
        GATE.authenticate = saved
    LAST['rc'] = rc
    return rc, buf.getvalue()


def gate(name, dirs, spec, expect_rc, substr, check=None, raw=None, cli=False, noauth=False, spec_sha=None, prog=None):
    """run the gate on the run directories with the spec (dict, or raw text); the outcome must be the expected exit status with the
    expected verdict as the last line, 'GATE PASS' printed iff rc 0, 'GATE FAIL' iff rc 1, and a report that records the same verdict and
    exit status.  check(report, output) -> (ok, text): an additional condition.  --spec-sha: the sha256 of the spec file written (unless
    spec_sha is given), --prog: the program of the spec's path (unless prog is given).  cli: run the CLI in a subprocess (else in-process);
    noauth: in-process with the authenticate step replaced (the report must then record NOAUTH_TAG for every run)."""
    sp = os.path.join(gdir(), 'spec_%s.json' % name)
    with open(sp, 'w') as f:
        f.write(raw if raw is not None else json.dumps(spec, indent=1))
    rp = os.path.join(gdir(), 'rep_%s.json' % name)
    path = spec.get('path', 'def') if isinstance(spec, dict) and spec.get('path') in ('def', 'log') else 'def'
    argv = ['--spec', sp, '--spec-sha', spec_sha or sha_file(sp), '--prog', prog or qbin(path), '--report', rp] + list(dirs)
    assert not (cli and noauth)
    rc, out = py('t207_d5_gate.py', *argv) if cli else gate_inproc(argv, noauth)
    lines = out.strip().split('\n')
    rep = json.load(open(rp)) if os.path.exists(rp) else None
    ok = (rc == expect_rc and substr in out and lines[-1].startswith(GW.get(expect_rc, '?')) and ('GATE PASS' in out) == (expect_rc == 0)
          and ('GATE FAIL' in out) == (expect_rc == 1) and rep is not None and rep.get('exit_code') == rc
          and rep.get('verdict') == GW.get(rc, '? ?').split()[1])
    if ok and noauth and 'authenticity' in rep:
        ok = all(a['source'] == NOAUTH_TAG for a in rep['authenticity'])
    extra = ('[%s] ' % ('cli' if cli else 'in-process, reprocessing replaced' if noauth else 'in-process'))
    if ok and check is not None:
        ok, e2 = check(rep, out)
        extra += e2
    record(name, 'rc %d: %s' % (expect_rc, substr[:28]), ok, lines[-1][:150] + ' | ' + extra)
    return rc, out, rep


def prog_of(tag):
    return json.load(open(os.path.join(base(tag), 'qref_proc.json')))['prog']


def seeds_of(tag):
    return json.load(open(os.path.join(base(tag), 'provenance.json')))['seeds_sha256']


def prov_of(tag, key):
    return json.load(open(os.path.join(base(tag), 'provenance.json')))[key]


def gate_tests():
    os.makedirs(gdir())
    B = gcopy('base')
    M3 = B[:3]
    Q = lambda d: os.path.join(d, 'qref.txt')                                                    # noqa: E731
    # ---- baseline: the kept dev5 subset through the test_only spec; every number independently from agg's CLI on the same files
    rc, msg = py('t207_d5_agg.py', *agg_args(os.path.join(ST, 'agg_base.json')))
    record('baseline_agg_def_dev1-4', 'ok', rc == 0, 'agg accepted all four runs' if rc == 0 else last(msg))
    out = os.path.join(ST, 'agg_gate.json')
    rca, msga = py('t207_d5_agg.py', out, 'pTT', '--ref', 'def', *(['main:%s=%s' % (t, Q(d)) for t, d in zip(GSET[:3], M3)] + ['ctl:c1=' + Q(B[3])]))
    record('baseline_agg_gate_set', 'ok', rca == 0, 'agg accepted dev2-4 and the control run' if rca == 0 else last(msga))
    A = json.load(open(out))

    def same_as_agg(rep, o):
        g = rep['gates'][0]
        c = {x['name']: x['checks'][0] for x in rep['controls']}
        pairs = [(g['device_value'], A['pooled']['dev']['S_ens']['mean'])] + \
                [(v, A['runs'][t]['dev']['S_ens']['mean']) for v, t in zip(g['device_run_means'], GSET[:3])] + \
                [(c['var8']['value'], A['pooled']['var8']['S_ens']['mean']), (c['tau_ns0']['value'], A['pooled']['ctl_dev']['S_ens']['mean'])]
        dmax = max(abs(a - b) for a, b in pairs)
        sc = rep['scripts']
        hashes = all(sc[f] == sha_file(os.path.join(HERE, f)) for f in ('t207_d5_gate.py', 't207_d5_agg.py', 't207_d5_qparse.py', 't207_d5_qref.F90'))
        ok = dmax <= 1e-12 and hashes and [r['role'] for r in rep['runs']] == ['main'] * 3 + ['ctl_tau_ns0'] and 'control var8' in o and 'control tau_ns0' in o
        # fix round 4: the spec hash (of the file given), the ALF binary of every run and the re-processing of every kept dump are in the report
        sp_ok = rep['spec_sha256'] == sha_file(os.path.join(gdir(), 'spec_baseline_gate_test_only_kept_subset.json'))
        bin_ok = all(r['binary_sha256'] == prov_of('m_def_dev2', 'binary_sha256') for r in rep['runs'])
        au = rep['authenticity']
        au_ok = len(au) == 4 and all(a['dumps'] == 7 and a['numbers'] > 0 and a['equal'] == a['numbers'] and a['source'] == 'program' for a in au)
        ok = ok and sp_ok and bin_ok and au_ok
        return ok, 'device %+.4f, var8 %+.4f, tau_ns0 %+.4f; max |gate - agg| %.1e; script hashes %s; spec sha %s, binary %s; re-processed %s' % (
            g['device_value'], c['var8']['value'], c['tau_ns0']['value'], dmax, 'ok' if hashes else 'WRONG', 'ok' if sp_ok else 'WRONG',
            'ok' if bin_ok else 'WRONG', ' '.join('%d/%d' % (a['numbers'], a['equal']) for a in au))
    rc, o, rep0 = gate('baseline_gate_test_only_kept_subset', B, gspec(), 4, 'GATE TEST-ONLY', check=same_as_agg, cli=True)
    # the same through pool each_run: the expected outcome is derived from agg's per-run numbers (independent of the gate's flow)
    dev_max = max(A['runs'][t]['dev']['S_ens']['mean'] for t in GSET[:3])
    v8_min = min(A['runs'][t]['var8']['S_ens']['mean'] for t in GSET[:3])
    ctl_min = A['runs']['c1']['ctl_dev']['S_ens']['mean']
    exp = 4 if (dev_max <= BOUND and v8_min > BOUND and ctl_min > BOUND) else 1
    gate('gate_pool_each_run_matches_agg_runs', B, smut(lambda s: s['gates'][0].__setitem__('pool', 'each_run')), exp, GW[exp],
         check=lambda r, o: (abs(r['gates'][0]['device_value'] - dev_max) <= 1e-12, 'device max run %+.4f, var8 min run %+.4f, ctl %+.4f' % (dev_max, v8_min, ctl_min)))
    # the LOG path: m_log_dev2..4 + c_log_dev1 (TAUSWAP=8) through a LOG test_only spec; outcome and numbers from agg's CLI
    LS = ['m_log_dev2', 'm_log_dev3', 'm_log_dev4', 'c_log_dev1_q20']
    L = gcopy('base_log', tags=LS)
    outl = os.path.join(ST, 'agg_gate_log.json')
    rcl, msgl = py('t207_d5_agg.py', outl, 'pTT', '--ref', 'log', *(['main:%s=%s' % (t, Q(d)) for t, d in zip(LS[:3], L[:3])] + ['ctl:c1=' + Q(L[3])]))
    record('baseline_agg_gate_set_log', 'ok', rcl == 0, 'agg accepted log dev2-4 and the control run' if rcl == 0 else last(msgl))
    AL = json.load(open(outl))
    lspec = gspec(path='log', refs=QP.ref_names('log', 'nvh'), nul=list(NUL_LOG), prog_sha=prog_of('m_log_dev2'), binary_sha256=prov_of('m_log_dev2', 'binary_sha256'),
                  seeds={'main': [seeds_of(t) for t in LS[:3]], 'ctl_tauswap8': [seeds_of(LS[3])]},
                  env_roles={'main': gspec()['env_roles']['main'], 'ctl_tauswap8': {'extra_env': {'ALF_GPU_TAUSWAP': '8'}, 'omit_env': ['ALF_GPU_TAUDIAG_QVAR']}},
                  n_runs={'main': 3, 'ctl_tauswap8': 1}, members={'main': mlist(NUL_LOG), 'ctl_tauswap8': mlist(NUL_LOG, var=False)},
                  controls=[dict(name='var8', kind='variant', member='var8', trips='acc_pTT'), dict(name='tauswap8', kind='run', role='ctl_tauswap8', trips='acc_pTT')])
    lv = (AL['pooled']['dev']['S_ens']['mean'], AL['pooled']['var8']['S_ens']['mean'], AL['pooled']['ctl_dev']['S_ens']['mean'])
    expl = 4 if (lv[0] <= BOUND and lv[1] > BOUND and lv[2] > BOUND) else 1

    def same_as_agg_log(rep, o):
        c = {x['name']: x['checks'][0]['value'] for x in rep['controls']}
        dmax = max(abs(a - b) for a, b in zip((rep['gates'][0]['device_value'], c['var8'], c['tauswap8']), lv))
        au = rep.get('authenticity', [])
        au_ok = len(au) == 4 and all(a['dumps'] == 7 and a['equal'] == a['numbers'] > 0 for a in au)
        return dmax <= 1e-12 and au_ok, 'device %+.4f, var8 %+.4f, tauswap8 %+.4f; max |gate - agg| %.1e; re-processed %s' % (
            lv + (dmax, ' '.join('%d/%d' % (a['numbers'], a['equal']) for a in au)))
    gate('baseline_gate_test_only_kept_subset_log', L, lspec, expl, GW[expl], check=same_as_agg_log, cli=True)
    # ---- modes: cover, test_only, controls
    gate('gate_kept_dev5_subset_cover1_spec', B, smut(lambda s: s.pop('test_only'), cover=1), 2, 'processed with cover 40')
    gate('gate_cover40_without_test_only', B, smut(lambda s: s.pop('test_only')), 2, 'must cover every dump')
    gate('gate_test_only_empty_reason', B, gspec(test_only=' '), 2, 'test_only must be a non-empty reason')
    nocl = dict(controls=[], n_runs={'main': 3}, members={'main': mlist(NUL_DEF)}, env_roles={'main': gspec()['env_roles']['main']},
                seeds={'main': gspec()['seeds']['main']})
    gate('P6_no_controls_dev_only', M3, gspec(**nocl), 3, 'GATE DEV-ONLY: no control declared')
    gate('P6_no_controls_ctl_run_given', B, gspec(**nocl), 2, 'env allow-list')
    gate('R_C_control_declared_run_missing', M3, gspec(), 2, 'declared in the spec, but no run has its role')
    gate('P2_main_run_offered_as_control', M3, gspec(n_runs={'main': 2, 'ctl_tau_ns0': 1}), 2, 'declared in the spec, but no run has its role')
    two = lambda s: (s['gates'].append(dict(name='acc_emax', measure='emax', statistic='S_ens', pool='pooled', bound=BOUND)),     # noqa: E731
                     s['expect_points'].__setitem__('emax', 7))
    gate('gate_uncovered_gate_dev_only_emax_wired', B, smut(two), 3, 'gates without a control that trips them',
         check=lambda r, o: (len(r['gates']) == 2 and r['gates'][1]['passed'], 'emax device %+.4f' % r['gates'][1]['device_value']))
    for m in ('emax', 'eF'):
        gate('gate_%s_device_fail' % m, B, smut(lambda s: (s['gates'].append(dict(name='acc_' + m, measure=m, statistic='S_ens', pool='pooled', bound=-1.0)),
                                                         s['expect_points'].__setitem__(m, 7))), 1, 'gate acc_%s: device' % m)
    # ---- bounds and the decision
    gate('gate_bound_below_device', B, smut(lambda s: s['gates'][0].__setitem__('bound', -0.5)), 1, 'gate acc_pTT: device')
    gate('gate_bound_above_controls', B, smut(lambda s: s['gates'][0].__setitem__('bound', 0.5)), 1, 'NOT KILLED')
    mnt = lambda b: (lambda s: (s['gates'].append(dict(name='acc_pTT2', measure='pTT', statistic='S_ens', pool='pooled', bound=b)),   # noqa: E731
                                s['controls'][1].__setitem__('must_not_trip', ['acc_pTT2'])))
    gate('gate_control_trips_gate_it_must_not', B, smut(mnt(BOUND)), 1, 'control tau_ns0 trips gate acc_pTT2')
    gate('gate_control_must_not_trip_holds', B, smut(mnt(5.0)), 3, 'gates without a control that trips them: [\'acc_pTT2\']',
         check=lambda r, o: (r['controls'][1]['checks'][1]['ok'] is True and r['controls'][1]['checks'][1]['must_trip'] is False, 'must_not_trip check ok'))
    # ---- spec structure and values (each: INPUT_ERROR with its own message)
    rp =os.path.join(gdir(), 'rep_gate_missing_spec_file.json')
    rc, o = py('t207_d5_gate.py', '--spec', os.path.join(gdir(), 'does_not_exist.json'), '--spec-sha', '0' * 64, '--prog', qbin('def'), '--report', rp, *B)
    record('gate_missing_spec_file', 'rc 2: gate spec missing', rc == 2 and 'gate spec missing' in o and json.load(open(rp))['exit_code'] == 2, last(o))
    for name, raw, substr in (('gate_truncated_spec', '{"spec_version": 1, ', 'is not valid JSON'),
                              ('gate_spec_top_level_list', '[1, 2]', 'not a JSON object'),
                              ('gate_spec_duplicate_key', json.dumps(gspec(), indent=1).replace('"qs": 20,', '"qs": 20, "qs": 20,', 1), 'duplicate key'),
                              ('gate_bound_NaN_token', json.dumps(gspec(), indent=1).replace('"bound": 0.15', '"bound": NaN', 1), 'non-finite constant NaN'),
                              ('gate_bound_Infinity_token', json.dumps(gspec(), indent=1).replace('"bound": 0.15', '"bound": Infinity', 1), 'non-finite constant')):
        assert raw.count('NaN') + raw.count('Infinity') <= 1
        gate(name, B, None, 2, substr, raw=raw)
    agg_as_spec = json.dumps(A)
    gate('P8b_aggregate_passed_as_spec', B, None, 2, 'allow_v1 is refused', raw=agg_as_spec)       # the aggregate carries allow_v1: refused first
    gate('P8c_aggregate_without_allow_v1_as_spec', B, None, 2, 'unknown keys', raw=json.dumps({k: v for k, v in A.items() if k != 'allow_v1'}))
    spec_tests = [
        ('P8_allow_v1_true', lambda s: s.__setitem__('allow_v1', True), 'allow_v1 is refused'),
        ('P8_allow_v1_false', lambda s: s.__setitem__('allow_v1', False), 'allow_v1 is refused'),
        ('gate_spec_unknown_key', lambda s: s.__setitem__('allow_cover', 40), 'unknown keys'),
        ('gate_spec_version_2', lambda s: s.__setitem__('spec_version', 2), 'spec_version'),
        ('gate_spec_path_unknown', lambda s: s.__setitem__('path', 'foo'), 'path'),
        ('gate_no_expect_points', lambda s: s.pop('expect_points'), "missing keys ['expect_points']"),
        ('gate_no_qs', lambda s: s.pop('qs'), "missing keys ['qs']"),
        ('R_A2_prog_sha_absent', lambda s: s.pop('prog_sha'), "missing keys ['prog_sha']"),
        ('gate_no_controls_key', lambda s: s.pop('controls'), "missing keys ['controls']"),
        ('gate_no_n_runs', lambda s: s.pop('n_runs'), "missing keys ['n_runs']"),
        # gates
        ('gate_wrong_measure_pmax', lambda s: s['gates'][0].__setitem__('measure', 'pmax'), "unknown measure 'pmax'"),
        ('gate_unknown_measure', lambda s: s['gates'][0].__setitem__('measure', 'foo'), "unknown measure 'foo'"),
        ('P4_statistic_S_cpu', lambda s: s['gates'][0].__setitem__('statistic', 'S_cpu'), 'not a decision statistic'),
        ('P4c_statistic_missing', lambda s: s['gates'][0].pop('statistic'), 'exactly the keys'),
        ('gate_pool_unknown', lambda s: s['gates'][0].__setitem__('pool', 'median'), "pool 'median'"),
        ('gate_REAL_crash_bound_stats_None', lambda s: s['gates'][0].__setitem__('bound', None), 'bound None is not a finite number'),
        ('gate_bound_a_string', lambda s: s['gates'][0].__setitem__('bound', '0.15'), "bound '0.15' is not a finite number"),
        ('gate_bound_true', lambda s: s['gates'][0].__setitem__('bound', True), 'bound True is not a finite number'),
        ('gate_no_bound', lambda s: s['gates'][0].pop('bound'), 'exactly the keys'),
        ('gate_gates_empty', lambda s: s.__setitem__('gates', []), 'gates must be a non-empty list'),
        ('gate_gates_not_a_list', lambda s: s.__setitem__('gates', 'acc_pTT'), 'gates must be a non-empty list'),
        ('gate_gate_extra_key', lambda s: s['gates'][0].__setitem__('note', 'x'), 'exactly the keys'),
        ('gate_gate_names_duplicate', lambda s: s['gates'].append(dict(s['gates'][0])), 'is not a unique name'),
        # expect_points
        ('gate_n_used_vs_expect_points_8', lambda s: s['expect_points'].__setitem__('pTT', 8), 'points used for measure pTT'),
        ('gate_n_used_short_expect_6', lambda s: s['expect_points'].__setitem__('pTT', 6), 'points used for measure pTT'),
        ('P7_expect_points_float', lambda s: s['expect_points'].__setitem__('pTT', 7.0), 'pTT: 7.0 is not an integer >= 1'),
        ('gate_REAL_crash_runs_entry_n_a_list', lambda s: s['expect_points'].__setitem__('pTT', [7]), 'pTT: [7] is not an integer >= 1'),
        ('gate_expect_points_true', lambda s: s['expect_points'].__setitem__('pTT', True), 'pTT: True is not an integer >= 1'),
        ('gate_expect_points_zero', lambda s: s['expect_points'].__setitem__('pTT', 0), 'pTT: 0 is not an integer >= 1'),
        ('gate_expect_points_extra_measure', lambda s: s['expect_points'].__setitem__('emax', 7), 'exactly the measures of the gates'),
        # qs (one qs per gate run), nul, refs
        ('gate_qs_zero', lambda s: s.__setitem__('qs', 0), 'qs 0: one integer > 0'),
        ('gate_qs_a_list', lambda s: s.__setitem__('qs', [20]), 'qs [20]: one integer > 0'),
        ('gate_qs_two_values', lambda s: s.__setitem__('qs', [20, 40]), 'qs [20, 40]: one integer > 0'),
        ('gate_qs_not_an_integer', lambda s: s.__setitem__('qs', '20'), "qs '20': one integer > 0"),
        ('gate_qs_float', lambda s: s.__setitem__('qs', 20.0), 'qs 20.0: one integer > 0'),
        ('gate_qs_not_passed_40', lambda s: s.__setitem__('qs', 40), 'processed with qs 20, the spec says 40'),
        ('gate_nul_without_0', lambda s: s['nul'].remove(0), 'containing 0'),
        ('gate_nul_misses_ref_code', lambda s: s['nul'].remove(1608), 'processed with nul'),
        ('gate_nul_duplicate', lambda s: s['nul'].append(4), 'distinct integers'),
        ('gate_nul_extra_code_vs_runs', lambda s: s['nul'].append(7), 'processed with nul'),
        ('gate_refs_unknown_name', lambda s: s['refs'].append('foo'), 'every reference member'),
        ('gate_refs_dev_listed', lambda s: s['refs'].append('dev'), 'is not a reference member'),
        ('gate_refs_n0_listed', lambda s: s['refs'].append('n0'), 'is not a reference member'),
        ('gate_refs_qself_listed', lambda s: s['refs'].append('qself'), 'is not a reference member'),
        # fix round 4: var<m> are never references (refused before the control check, which named it before)
        ('gate_refs_variant_control_member', lambda s: (s['refs'].append('var8'), s['members']['ctl_tau_ns0'].insert(1, 'var8')),
         "refs: 'var8' is not a reference member"),
        ('gate_refs_one_name', lambda s: s.__setitem__('refs', ['cpu']), 'at least two names'),
        # member lists (names from the spec; strictly exactly those, in that order)
        ('gate_members_order_swapped', lambda s: s['members']['main'].__setitem__(slice(4, 6), s['members']['main'][5:3:-1]), 'differs from the member list of the spec'),
        ('gate_members_one_missing', lambda s: s['members']['main'].remove('u2'), 'differs from the member list of the spec'),
        ('gate_members_extra', lambda s: s['members']['main'].append('zz'), 'differs from the member list of the spec'),
        ('gate_members_ctl_list_with_var8', lambda s: s['members']['ctl_tau_ns0'].insert(1, 'var8'), 'differs from the member list of the spec'),
        ('gate_host_member_other_label', lambda s: s.update(host_member='host', members={r: ['host' if x == 'nvh' else x for x in l] for r, l in s['members'].items()},
                                                            refs=['host' if x == 'nvh' else x for x in s['refs']]), "member 'nvh' where 'host' was expected"),
        ('gate_host_member_not_in_list', lambda s: s.__setitem__('host_member', 'host'), 'does not contain host_member host'),
        ('gate_members_role_missing', lambda s: s['members'].pop('ctl_tau_ns0'), 'must have exactly the roles'),
        ('gate_members_no_dev', lambda s: s['members']['main'].remove('dev'), 'must contain dev'),
        ('gate_members_diag_listed', lambda s: s['members']['main'].append('qself'), 'no diagnostic member'),
        ('gate_members_not_a_dict', lambda s: s.__setitem__('members', ['dev']), 'members must be a dictionary'),
        ('gate_refs_duplicate', lambda s: s['refs'].append('n4'), 'duplicate name'),
        ('gate_refs_empty', lambda s: s.__setitem__('refs', []), 'at least two names'),
        ('gate_cover_zero', lambda s: s.__setitem__('cover', 0), 'cover 0: an integer >= 1'),
        # fix round 5 (M1): keep, hostin, min_reprocessed are required spec keys with exact integer values
        ('M1_spec_no_keep', lambda s: s.pop('keep'), "missing keys ['keep']"),
        ('M1_spec_no_hostin', lambda s: s.pop('hostin'), "missing keys ['hostin']"),
        ('M1_spec_no_min_reprocessed', lambda s: s.pop('min_reprocessed'), "missing keys ['min_reprocessed']"),
        ('M1_spec_keep_float', lambda s: s.__setitem__('keep', 0.0), 'keep 0.0: the processing KEEP'),
        ('M1_spec_keep_true', lambda s: s.__setitem__('keep', True), 'keep True: the processing KEEP'),
        ('M1_spec_keep_negative', lambda s: s.__setitem__('keep', -1), 'keep -1: the processing KEEP'),
        ('M1_spec_hostin_2', lambda s: s.__setitem__('hostin', 2), 'hostin 2: the processing HOSTIN'),
        ('M1_spec_hostin_string', lambda s: s.__setitem__('hostin', '0'), "hostin '0': the processing HOSTIN"),
        ('M1_spec_min_reprocessed_zero', lambda s: s.__setitem__('min_reprocessed', 0), 'min_reprocessed 0: an integer >= 1'),
        ('M1_spec_min_reprocessed_float', lambda s: s.__setitem__('min_reprocessed', 7.0), 'min_reprocessed 7.0: an integer >= 1'),
        # the processing arguments of every run must be the spec's: keep, hostin (consistent runs, other spec value)
        ('M1_spec_keep_differs_from_runs', lambda s: s.__setitem__('keep', 3), 'processed with keep 0, the spec says 3'),
        ('M1_spec_hostin_differs_from_runs', lambda s: s.__setitem__('hostin', 1), 'processed with hostin 0, the spec says 1'),
        ('M1_runs_fewer_dumps_than_min_reprocessed', lambda s: s.__setitem__('min_reprocessed', 8), 'the spec requires at least min_reprocessed 8'),
        ('gate_cover_float', lambda s: s.__setitem__('cover', 40.0), 'cover 40.0: an integer >= 1'),
        # program identity
        ('gate_prog_sha_mismatch', lambda s: s.__setitem__('prog_sha', '0' * 64), 'differs from the required'),
        ('gate_prog_sha_not_a_hash', lambda s: s.__setitem__('prog_sha', 'abc'), 'is not a sha256'),
        ('gate_prog_sha_upper_case', lambda s: s.__setitem__('prog_sha', s['prog_sha'].upper()), 'is not a sha256'),
        ('gate_src_sha_not_this_source', lambda s: s.__setitem__('src_sha', '0' * 64), 'of the t207_d5_qref.F90 next to the gate'),
        # env table
        ('gate_env_roles_no_main', lambda s: s['env_roles'].pop('main'), 'role main is required'),
        ('gate_env_roles_main_no_QREF', lambda s: s['env_roles']['main'].pop('ALF_GPU_TAUDIAG_QREF'), 'ALF_GPU_TAUDIAG_QREF missing'),
        ('gate_env_roles_override', lambda s: s['env_roles']['ctl_tau_ns0']['extra_env'].__setitem__('ALF_GPU_TAUDIAG_QHOST', '0'), 'no override'),
        ('gate_env_roles_omit_unknown', lambda s: s['env_roles']['ctl_tau_ns0']['omit_env'].append('ALF_GPU_FOO'), 'omit_env'),
        ('gate_env_roles_value_int', lambda s: s['env_roles']['main'].__setitem__('ALF_GPU_DEBUG', 1), 'non-empty string'),
        ('gate_env_roles_key_not_ALF_GPU', lambda s: s['env_roles']['main'].__setitem__('FOO', '1'), 'not an ALF_GPU_* variable'),
        ('gate_env_roles_two_roles_same_env', lambda s: (s['env_roles'].__setitem__('ctl2', json.loads(json.dumps(s['env_roles']['ctl_tau_ns0']))),
                                                         s['controls'].append(dict(name='c2', kind='run', role='ctl2', trips='acc_pTT')),
                                                         s['n_runs'].__setitem__('ctl2', 1)), 'have the same env'),
        ('gate_env_roles_unused_role', lambda s: s['env_roles'].__setitem__('ctl_swap', {'extra_env': {'ALF_GPU_TAUSWAP': '8'}, 'omit_env': []}), 'used by no control'),
        ('gate_env_roles_main_QVAR_not_a_list', lambda s: (s['env_roles']['main'].__setitem__('ALF_GPU_TAUDIAG_QVAR', 'x'),
                                                           s['env_roles']['ctl_tau_ns0'].__setitem__('omit_env', [])), 'not a comma list of integers'),
        ('gate_env_roles_entry_extra_key', lambda s: s['env_roles']['ctl_tau_ns0'].__setitem__('note', 'x'), 'exactly the keys extra_env and omit_env'),
        ('gate_env_roles_empty_control_role', lambda s: s['env_roles'].__setitem__('ctl_tau_ns0', {'extra_env': {}, 'omit_env': []}), 'both empty'),
        # run counts
        ('gate_n_runs_main_4', lambda s: s['n_runs'].__setitem__('main', 4), 'runs per role'),
        ('P7_n_runs_float', lambda s: s['n_runs'].__setitem__('main', 3.0), 'main: 3.0 is not an integer >= 1'),
        ('gate_n_runs_missing_role', lambda s: s['n_runs'].pop('ctl_tau_ns0'), 'must have exactly the roles'),
        ('gate_n_runs_extra_role', lambda s: s['n_runs'].__setitem__('ctl_x', 1), 'must have exactly the roles'),
        ('gate_n_runs_zero', lambda s: s['n_runs'].__setitem__('main', 0), 'main: 0 is not an integer >= 1'),
        # controls
        ('gate_controls_not_a_list', lambda s: s.__setitem__('controls', 'var8'), 'controls must be a list'),
        ('gate_control_entry_a_list', lambda s: s['controls'].__setitem__(0, ['var8']), 'is not a dictionary'),
        ('gate_control_kind_unknown', lambda s: s['controls'][0].__setitem__('kind', 'foo'), 'kind variant'),
        ('P4b_control_statistic_key', lambda s: s['controls'][1].__setitem__('statistic', 'S_cpu'), 'kind variant'),
        ('gate_control_member_var16', lambda s: s['controls'][0].__setitem__('member', 'var16'), 'must be a member of main'),
        ('gate_control_member_dev', lambda s: s['controls'][0].__setitem__('member', 'dev'), 'member dev is the device under test'),
        # fix round 5 (N3): the branch 'control member is a reference' (unreachable via var8 since var<m> are never references)
        ('N3_control_member_is_reference_n4', lambda s: s['controls'][0].__setitem__('member', 'n4'), "member 'n4' is a reference member"),
        ('gate_control_trips_unknown_gate', lambda s: s['controls'][0].__setitem__('trips', 'nope'), 'is not a gate of the spec'),
        ('gate_control_must_not_trip_its_own', lambda s: s['controls'][0].__setitem__('must_not_trip', ['acc_pTT']), 'list of distinct gates other than'),
        ('gate_control_run_role_main', lambda s: s['controls'][1].__setitem__('role', 'main'), 'is not a control role'),
        ('gate_control_run_role_unknown', lambda s: s['controls'][1].__setitem__('role', 'ghost'), 'is not a control role'),
        ('gate_control_names_duplicate', lambda s: s['controls'][1].__setitem__('name', 'var8'), 'is not a unique name'),
        # P1 / P1b: the default-path runs judged as LOG
        ('P1_def_runs_spec_says_log', lambda s: s.update(path='log', refs=QP.ref_names('log', 'nvh'), prog_sha=prog_of('m_log_dev2'),
                                                          env_roles=dict(s['env_roles'], ctl_tau_ns0={'extra_env': {'ALF_GPU_TAUSWAP': '8'}, 'omit_env': ['ALF_GPU_TAUDIAG_QVAR']})),
         'a def-path run (info)'),
        # fix round 4: the control role gets the LOG fault as in P1 (TAU_NS=0 is a default-path fault; on a LOG spec it is refused first)
        ('P1b_def_runs_spec_says_log_def_codes_removed', lambda s: s.update(path='log', refs=QP.ref_names('log', 'nvh'), nul=list(NUL_LOG), prog_sha=prog_of('m_log_dev2'),
                                                                            env_roles=dict(s['env_roles'], ctl_tau_ns0={'extra_env': {'ALF_GPU_TAUSWAP': '8'}, 'omit_env': ['ALF_GPU_TAUDIAG_QVAR']})),
         'a def-path run (info)'),
    ]
    for name, fn, substr in spec_tests:
        gate(name, B, smut(fn), 2, substr)
    # ---- the reference list is the spec's (manager note 2026-11-01: names come from the spec, no frozen list in the scripts): another list
    # is accepted and judged with it; the expected outcome and numbers come from agg's CLI with the same --ref list
    for name, refs in (('gate_refs_member_dropped', [r for r in QP.ref_names('def', 'nvh') if r != 'n1608']),
                       ('gate_refs_reordered', QP.ref_names('def', 'nvh')[::-1]),
                       ('gate_refs_other_path_list', QP.ref_names('log', 'nvh'))):
        oa = os.path.join(ST, 'agg_%s.json' % name)
        rcr, msgr = py('t207_d5_agg.py', oa, 'pTT', '--ref', ','.join(refs), *(['main:%s=%s' % (t, Q(d)) for t, d in zip(GSET[:3], M3)] + ['ctl:c1=' + Q(B[3])]))
        Ar = json.load(open(oa)) if rcr == 0 else None
        if Ar is None:
            record(name, 'agg ok', False, last(msgr))
            continue
        vr = (Ar['pooled']['dev']['S_ens']['mean'], Ar['pooled']['var8']['S_ens']['mean'], Ar['pooled']['ctl_dev']['S_ens']['mean'])
        er = 4 if (vr[0] <= BOUND and vr[1] > BOUND and vr[2] > BOUND) else 1

        def same(rep, o, vr=vr):
            c = {x['name']: x['checks'][0]['value'] for x in rep['controls']}
            dm = max(abs(a - b) for a, b in zip((rep['gates'][0]['device_value'], c['var8'], c['tau_ns0']), vr))
            return dm <= 1e-12, 'device %+.4f var8 %+.4f tau_ns0 %+.4f, max |gate - agg| %.1e' % (vr + (dm,))
        gate(name, B, gspec(refs=refs), er, GW[er], check=same)
    # a host member printed under another label (renamed in every qref.txt line, processing record kept consistent) is accepted when the
    # spec declares that label; the numbers are those of the baseline
    ren = e_qref(lambda t: t.replace(' | nvh ', ' | hostx '))
    hs = gspec(host_member='hostx', members={r: ['hostx' if x == 'nvh' else x for x in l] for r, l in gspec()['members'].items()},
               refs=['hostx' if x == 'nvh' else x for x in QP.ref_names('def', 'nvh')])
    hx = gcopy('hostx', edits={t: ren for t in GSET})
    gate('gate_host_member_relabelled_accepted', hx, hs, 4, 'GATE TEST-ONLY', noauth=True,
         check=lambda r, o: (abs(r['gates'][0]['device_value'] - rep0['gates'][0]['device_value']) <= 1e-12, 'device %+.4f as the baseline' % r['gates'][0]['device_value']))
    # fix round 4: the program prints 'nvh'; a relabelled qref.txt is not its output (the label is the spec's only for a program that prints it)
    gate('gate_host_member_relabelled_caught_by_reprocessing', hx, hs, 2, "member 'hostx ")
    gate('gate_host_member_relabelled_spec_says_nvh', gcopy('hostx2', edits={t: ren for t in GSET}), gspec(), 2, "member 'hostx' where 'nvh' was expected")
    # fix round 5, manager item X1 (Intel calibration): the host label comes from the spec's host_member only.  The NVIDIA records relabelled
    # consistently to the Intel label ifh (every qref.txt line, qref_sha256 following), spec host_member ifh with the frozen list
    # (refs = FROZEN_REFS with H = ifh, no test_only needed for the list): the same verdict and numbers as the baseline (reprocessing replaced:
    # the NVIDIA build of the program prints its own label, so the real step must refuse the relabelled lines, companion test); and agg's
    # CLI on the same files (--ref def, host label taken from the runs) gives the baseline numbers
    ren_i = e_qref(lambda t: t.replace(' | nvh ', ' | ifh '))
    hi_ = gspec(host_member='ifh', members={r: ['ifh' if x == 'nvh' else x for x in l] for r, l in gspec()['members'].items()},
                refs=QP.ref_names('def', 'ifh'))
    hxi = gcopy('host_ifh', edits={t: ren_i for t in GSET})
    gate('X1_host_member_ifh_relabelled_same_verdict', hxi, hi_, rep0['exit_code'], GW[rep0['exit_code']], noauth=True,
         check=lambda r, o: (r['verdict'] == rep0['verdict'] and max(abs(a['checks'][0]['value'] - b['checks'][0]['value']) for a, b in zip(r['controls'], rep0['controls'])) <= 1e-12
                             and abs(r['gates'][0]['device_value'] - rep0['gates'][0]['device_value']) <= 1e-12,
                             'verdict %s as the baseline, device %+.4f, controls %s' % (r['verdict'], r['gates'][0]['device_value'],
                                                                                         ' '.join('%+.4f' % c['checks'][0]['value'] for c in r['controls']))))
    gate('X1_host_member_ifh_relabelled_caught_by_reprocessing', hxi, hi_, 2, "member 'ifh ")
    oi = os.path.join(ST, 'agg_gate_ifh.json')
    rci, msgi = py('t207_d5_agg.py', oi, 'pTT', '--ref', 'def', *(['main:%s=%s' % (t, Q(d)) for t, d in zip(GSET[:3], hxi[:3])] + ['ctl:c1=' + Q(hxi[3])]))
    Ai = json.load(open(oi)) if rci == 0 else None
    dmi = None if Ai is None else max(abs(Ai['pooled'][m]['S_ens']['mean'] - A['pooled'][m if m != 'ifh' else 'nvh']['S_ens']['mean']) for m in Ai['pooled'])
    record('X1_agg_host_label_ifh_from_runs', 'agg ok, numbers = baseline', Ai is not None and Ai['host_member'] == 'ifh' and 'ifh' in Ai['refs'] and dmi <= 1e-12,
           last(msgi) if Ai is None else 'host %s, refs %s, max |agg ifh - agg nvh| %.1e' % (Ai['host_member'], ','.join(Ai['refs']), dmi))
    gate('X1_spec_host_member_reserved_name', B, gspec(host_member='cpu'), 2, "host_member")
    # 207-5 (round-5 review follow-ups): ctl_dev (agg's name of a control run's dev) is a reserved label; a kept dump that changes while the
    # program re-processes it is refused (re-hash at read time; the dump of the copy is first made a regular file with the same bytes, then
    # rewritten by a hook that runs right after the program, so nothing is written through a symbolic link into the dev5 data)
    gate('Y1_spec_host_member_ctl_dev_reserved', B, gspec(host_member='ctl_dev'), 2, "host_member")
    def to_regular(d, name):              # the symbolic link replaced by a regular file with the same bytes (manifest unchanged)
        data = dump_bytes(d, name)
        os.remove(os.path.join(d, 'qref', name))
        with open(os.path.join(d, 'qref', name), 'wb') as fh:
            fh.write(data)
    yd = gcopy('rehash')
    yname = sorted(os.listdir(os.path.join(yd[0], 'qref')))[0]
    to_regular(yd[0], yname)
    yfile = os.path.join(yd[0], 'qref', yname)
    ysaved = GATE.reprocess

    def yhook(x, s, prog):
        r = ysaved(x, s, prog)
        if x['dir'] == os.path.abspath(yd[0]):
            assert not os.path.islink(yfile)
            with open(yfile, 'r+b') as fh:
                fh.seek(os.path.getsize(yfile) - 8)
                fh.write(b'\x01' * 8)
        return r
    GATE.reprocess = yhook
    try:
        gate('Y2_dump_changed_during_reprocessing', yd, gspec(), 2, 'changed during the re-processing')
    finally:
        GATE.reprocess = ysaved
    gate('Y2b_same_copy_unchanged_baseline', gcopy('rehash_base', edits={GSET[0]: lambda d: to_regular(d, yname)}),
         gspec(), 4, 'GATE TEST-ONLY')
    # ---- command line (fix round 4: --spec-sha and --prog are required; every case below names them, so that only the defect named fails)
    sp0 = os.path.join(gdir(), 'spec_baseline_gate_test_only_kept_subset.json')
    hs0 = sha_file(sp0)
    pp = qbin('def')
    for name, args, substr in (('gate_cli_no_report', ['--spec', sp0, '--spec-sha', hs0, '--prog', pp] + B, 'usage'),
                               ('gate_cli_spec_twice', ['--spec', sp0, '--spec', sp0, '--spec-sha', hs0, '--prog', pp, '--report', os.path.join(gdir(), 'rep_cli2.json')] + B, 'given twice'),
                               ('gate_cli_unknown_option', ['--require-controls', '--spec', sp0, '--spec-sha', hs0, '--prog', pp, '--report', os.path.join(gdir(), 'rep_cli3.json')] + B, 'unknown option'),
                               ('gate_cli_no_run_dirs', ['--spec', sp0, '--spec-sha', hs0, '--prog', pp, '--report', os.path.join(gdir(), 'rep_cli4.json')], 'usage'),
                               ('gate_cli_no_spec_sha', ['--spec', sp0, '--prog', pp, '--report', os.path.join(gdir(), 'rep_cli5.json')] + B, 'usage'),
                               ('gate_cli_no_prog', ['--spec', sp0, '--spec-sha', hs0, '--report', os.path.join(gdir(), 'rep_cli6.json')] + B, 'usage'),
                               ('gate_cli_spec_sha_not_hex', ['--spec', sp0, '--spec-sha', 'abc', '--prog', pp, '--report', os.path.join(gdir(), 'rep_cli7.json')] + B, 'is not a sha256'),
                               ('gate_cli_spec_sha_upper_case', ['--spec', sp0, '--spec-sha', hs0.upper(), '--prog', pp, '--report', os.path.join(gdir(), 'rep_cli8.json')] + B, 'is not a sha256')):
        rc, o = py('t207_d5_gate.py', *args)
        record(name, 'rc 2: ' + substr, rc == 2 and substr in o and 'GATE FAIL' not in o and 'GATE PASS' not in o, last(o))
    rp = os.path.join(gdir(), 'rep_exists.json')
    open(rp, 'w').write('keep\n')
    rc, o = py('t207_d5_gate.py', '--spec', sp0, '--spec-sha', hs0, '--prog', pp, '--report', rp, *B)
    record('gate_report_exists_not_overwritten', 'rc 2, file unchanged', rc == 2 and 'exists' in o and open(rp).read() == 'keep\n', last(o))
    # ---- evidence: one defect per copy of the run set (INPUT_ERROR)
    D3 = 'm_def_dev3'
    ev = [
        ('P10_run_dir_nonexistent', None, B + [os.path.join(gdir(), 'nonexistent', 'hub12x12_cb_p')], 'not a directory'),
        ('gate_run_qref_txt_missing', {D3: e_file(lambda d: os.remove(os.path.join(d, 'qref.txt')))}, None, 'qref.txt missing'),
        ('gate_run_qref_edited_after_processing', {D3: e_qref_raw(lines_all(scale('u1', 0, 1.001), only=2))}, None, 'differs from qref_sha256'),
        ('gate_run_proc_json_missing', {D3: e_file(lambda d: os.remove(os.path.join(d, 'qref_proc.json')))}, None, 'qref_proc.json missing'),
        ('R_A_no_program_identity', {t: e_proc(lambda x: x.pop('prog')) for t in GSET}, None, 'no valid program sha256'),
        ('R_A3_prog_not_hex_in_run', {D3: e_proc(lambda x: x.__setitem__('prog', 'G' * 64))}, None, 'no valid program sha256'),
        ('gate_prog_sha_one_run_differs', {D3: e_proc(lambda x: x.__setitem__('prog', '1' * 64))}, None, 'program sha256 1111'),
        ('gate_run_prog_src_differs', {D3: e_proc(lambda x: x.__setitem__('prog_src', '2' * 64))}, None, 'program source sha256 2222'),
        ('gate_run_proc_rc_1', {D3: e_proc(lambda x: x.__setitem__('rc', 1))}, None, 'are not both 0'),
        ('gate_run_proc_check_rc_true', {D3: e_proc(lambda x: x.__setitem__('check_rc', True))}, None, 'are not both 0'),
        ('gate_run_proc_qref_sha256_missing', {D3: e_proc(lambda x: x.pop('qref_sha256'))}, None, 'no valid qref_sha256'),
        ('gate_run_proc_args_differ', {D3: e_proc(lambda x: x['args'].__setitem__('qs', 40))}, None, 'differ from the arguments in qref_spec.json'),
        ('P9_mixed_qs_one_run_40', {D3: e_args(qs=40)}, None, 'processed with qs 40, the spec says 20'),
        ('R_E_hostin_7', {D3: e_args(hostin=7)}, None, 'bad processing arguments'),
        ('R_E_hostin_true', {D3: e_args(hostin=True)}, None, 'must be integers'),
        ('R_E_keep_negative', {D3: e_args(keep=-5)}, None, 'bad processing arguments'),
        ('keep_rule_refused_by_gate', {D3: e_args(keep=1)}, None, 'keep rule'),
        # fix round 5: the processing arguments are checked before the manifest; a run with cover 2 among cover-40 runs is refused there
        # (gate_cover_2 itself, which tests the manifest check behind it, runs below with a cover-2 spec on that run alone)
        ('M1_one_run_cover_2_args_check', {D3: e_args(cover=2)}, None, 'processed with cover 2, the spec says 40'),
        ('gate_run_records_edited', {D3: e_file(lambda d: open(os.path.join(d, 'run.log'), 'a').write(' extra\n'))}, None, 'does not match the record run.log'),
        # fix round 4: manifest_sha256 of the processing record follows the edit (consistent), so that the dump hash check behind it is reached
        ('gate_run_manifest_hash_edited', {D3: chain(e_file(lambda d: rewrite(os.path.join(d, 'qref_manifest.txt'),
                                                                              lambda t: re.sub(r'^[0-9a-f]{64}', '0' * 64, t, count=1))), e_file(fix_msha))},
         None, 'sha256 differs from qref_manifest.txt'),
        ('gate_run_dump_removed', {D3: e_file(lambda d: os.remove(os.path.join(d, 'qref', sorted(os.listdir(os.path.join(d, 'qref')))[2])))}, None, 'qref/ holds 6 files'),
        ('gate_env_missing_stamp_line', {D3: e_file(lambda d: rewrite(os.path.join(d, 'stamp.txt'), lambda t: re.sub(r'^env .*\n', '', t, flags=re.M)))},
         None, 'stamp.txt: env line'),
        ('gate_run_point_dropped_ctl', {'c_def_dev1_q20': e_qref(lambda t: '\n'.join(t.split('\n')[:2] + t.split('\n')[3:]).replace(
            'QREFEND dumps 7 ok 7 err 0', 'QREFEND dumps 6 ok 6 err 0'))}, None, 'differs from dumps 7'),
        ('gate_run_dev_NaN', {D3: e_qref(lines_all(setnum('dev', 5, 'NaN'), only=1))}, None, 'not a number in the format ES12.4E3'),
        ('gate_run_dev_removed', {D3: e_qref(lines_all(drop_member('dev'), only=1))}, None, 'members, expected'),
        ('gate_run_var8_removed_all_points', {D3: e_qref(lines_all(drop_member('var8')))}, None, 'members, expected'),
        ('gate_variant_mixed_log_run', None, 'm_log_dev2', 'a log-path run (info)'),
        # env (consistent edits of info, stamp.txt, provenance.json and the stored spec)
        ('gate_env_TAUSWAP_main', {D3: e_env(lambda e: dict(e, ALF_GPU_TAUSWAP='8'))}, None, 'env allow-list'),
        ('gate_env_unknown_var_main', {D3: e_env(lambda e: dict(e, ALF_GPU_FOO='1'))}, None, 'env allow-list'),
        ('gate_env_DEBUG_removed_main', {D3: e_env(lambda e: {k: v for k, v in e.items() if k != 'ALF_GPU_DEBUG'})}, None, 'env allow-list'),
        ('gate_env_QVAR_removed_main', {D3: e_env(lambda e: {k: v for k, v in e.items() if k != 'ALF_GPU_TAUDIAG_QVAR'})}, None, 'differs from the env'),
        ('gate_env_QVAR_other_main', {D3: e_env(lambda e: dict(e, ALF_GPU_TAUDIAG_QVAR='8,16'))}, None, 'differs from the env'),
        ('gate_env_QHOST_0_main', {D3: e_env(lambda e: dict(e, ALF_GPU_TAUDIAG_QHOST='0'))}, None, 'differs from the env'),
        ('gate_env_QREF_not_a_number', {D3: e_env(lambda e: dict(e, ALF_GPU_TAUDIAG_QREF='x'))}, None, 'not a QREF run'),
        ('gate_env_FROM_other_main', {D3: e_env(lambda e: dict(e, ALF_GPU_TAUDIAG_FROM='330'))}, None, 'differs from the env'),
        ('gate_env_ctl_TAU_NS_1', {'c_def_dev1_q20': e_env(lambda e: dict(e, ALF_GPU_TAU_NS='1'))}, None, 'env allow-list'),
        ('gate_env_ctl_TAUSWAP_on_default_path', {'c_def_dev1_q20': e_env(lambda e: dict({k: v for k, v in e.items() if k != 'ALF_GPU_TAU_NS'}, ALF_GPU_TAUSWAP='8'))},
         None, 'env allow-list'),
        ('P2c_ctl_two_faults', {'c_def_dev1_q20': e_env(lambda e: dict(e, ALF_GPU_TAUSWAP='8'))}, None, 'env allow-list'),
        ('P2_ctl_run_with_main_env', {'c_def_dev1_q20': e_env(lambda e: dict({k: v for k, v in e.items() if k != 'ALF_GPU_TAU_NS'}, ALF_GPU_TAUDIAG_QVAR='8'))},
         None, 'differs from the env'),
        ('P11_ctl_same_seeds_as_main', {'c_def_dev1_q20': e_seeds(seeds_of('m_def_dev2'))}, None, 'same seeds'),
    ]
    for name, edits, extra, substr in ev:
        if isinstance(extra, list):
            dirs = extra
        elif isinstance(extra, str):
            dirs = gcopy(name, tags=GSET + [extra])
        else:
            dirs = gcopy(name, edits=edits)
        gate(name, dirs, gspec(), 2, substr)
    # the manifest check behind the arguments: run processed with cover 2 (args edited consistently, its manifest is the cover-40 one), spec cover 2
    gate('gate_cover_2', gcopy('gate_cover_2', tags=[D3], edits={D3: e_args(cover=2)}), gspec(cover=2), 2, 'covered dumps')
    gate('gate_runs_different_parameters', gcopy('gate_runs_different_parameters', edits={D3: chain(
        e_file(lambda d: open(os.path.join(d, 'parameters'), 'a').write('! a comment line: other parameters file\n')), e_file(fix_records))}),
        gspec(), 2, 'different parameters')
    gate('gate_no_main_run', B[3:], gspec(), 2, 'no main run')
    # a version 1 run (s_def: no B, Bi) through a spec that matches it in everything else (qs 1, cover 1, its env): refused by the parser
    v1env = dict(gspec()['env_roles']['main'], ALF_GPU_TAUDIAG_QREF='1', ALF_GPU_TAUDIAG_FROM='399')
    gate('gate_v1_run_in_set', gcopy('gate_v1_run_in_set', tags=['s_def']), gspec(qs=1, cover=1, env_roles={'main': v1env}, n_runs={'main': 1}, controls=[],
                                                                                  members={'main': mlist(NUL_DEF)}, seeds={'main': [seeds_of('s_def')]},
                                                                                  binary_sha256=prov_of('s_def', 'binary_sha256')), 2, 'only with allow_v1')
    # the catch-all: an unreadable record is an exception the gate does not test for; it must end as INPUT_ERROR naming the exception
    dirs = gcopy('gate_REAL_crash_unreadable_record')
    os.chmod(os.path.join(dirs[1], 'qref_spec.json'), 0)
    gate('gate_REAL_crash_unreadable_record', dirs, gspec(), 2, 'PermissionError')
    os.chmod(os.path.join(dirs[1], 'qref_spec.json'), 0o644)
    # ---- the same run twice (INPUT_ERROR)
    gate('gate_dup_same_dir_twice', B + [B[1]], gspec(n_runs={'main': 4, 'ctl_tau_ns0': 1}), 2, 'same directory')
    dd = gcopy('gate_dup_copied_dir', tags=GSET + ['dup3'], src={'dup3': D3})
    gate('gate_dup_copied_dir', dd, gspec(n_runs={'main': 4, 'ctl_tau_ns0': 1}), 2, 'same qref.txt sha256')
    dd = gcopy('gate_dup_edited_qref_same_records', tags=GSET + ['dup3'], src={'dup3': D3}, edits={'dup3': e_qref(lines_all(scale('u1', 0, 1.001), only=2))})
    gate('gate_dup_edited_qref_same_records', dd, gspec(n_runs={'main': 4, 'ctl_tau_ns0': 1}), 2, 'same run.log sha256')
    forge = chain(e_qref(lines_all(scale('u1', 0, 1.001), only=2)), e_file(lambda d: open(os.path.join(d, 'run.log'), 'a').write(' harmless extra line\n')),
                  e_file(fix_records), e_seeds('f' * 64))
    dd = gcopy('P5_dup_run_new_hashes', tags=GSET + ['dup3'], src={'dup3': D3}, edits={'dup3': forge})
    gate('P5_dup_run_new_hashes', dd, gspec(n_runs={'main': 4, 'ctl_tau_ns0': 1}), 2, 'same dump sha256')

    def own_dumps(d):                       # the dumps copied and changed by one byte each, manifest rewritten: only the QREF lines are shared
        # fix round 4: the changed byte is the lowest mantissa byte of the last double of the file (Bi), so that the header and size checks of
        # the dumps pass (appending a byte, as before, is now refused there) and the duplicate check behind them is reached; manifest_sha256 follows
        qd = os.path.join(d, 'qref')
        lines = []
        for f in sorted(os.listdir(qd)):
            p = os.path.join(qd, f)
            src = os.path.realpath(p)
            os.remove(p)
            shutil.copyfile(src, p)
            sz = os.path.getsize(p)
            with open(p, 'r+b') as g:
                g.seek(sz - 8)
                b0 = g.read(1)
                g.seek(sz - 8)
                g.write(bytes([b0[0] ^ 1]))
            lines.append('%s  qref/%s' % (sha_file(p), f))
        open(os.path.join(d, 'qref_manifest.txt'), 'w').write('\n'.join(lines) + '\n')
        fix_msha(d)
    dd = gcopy('P5b_dup_run_dumps_forged_too', tags=GSET + ['dup3'], src={'dup3': D3}, edits={'dup3': chain(forge, e_file(own_dumps))})
    gate('P5b_dup_run_dumps_forged_too', dd, gspec(n_runs={'main': 4, 'ctl_tau_ns0': 1}), 2, 'same QREF data line')
    shutil.rmtree(os.path.join(dd[-1], 'qref'))           # own temporary dump copies (40 MB); the record of the test stays
    # ---- planted defects in the evidence: the gate's criterion catches them (FAIL)
    # fix round 4: these planted defects are edits of qref.txt numbers.  The decision logic is tested with the authenticate step replaced
    # (noauth); every one has a companion '<name>_caught_by_reprocessing' with the real step, which must end in INPUT_ERROR (the numbers are
    # not what the program computes from the dumps).  Planted defects through the program itself: the dump tests (dev output overflow, NaN).
    RP = '(re-processing check)'
    up = {t: e_qref(lines_all(scale('dev', 5, 10.0))) for t in GSET[:3]}
    dd = gcopy('P3a', edits=up)
    gate('P3a_device_x10_all_runs', dd, gspec(), 1, 'gate acc_pTT: device', noauth=True)
    gate('P3a_device_x10_all_runs_caught_by_reprocessing', dd, gspec(), 2, RP)
    dd = gcopy('huge', edits={'m_def_dev2': e_qref(lines_all(setnum('dev', 5, '1.0000E+100'), only=1))})
    gate('huge_value_dev_pTT', dd, gspec(), 1, 'gate acc_pTT: device', noauth=True)
    gate('huge_value_dev_pTT_caught_by_reprocessing', dd, gspec(), 2, RP)
    dd = gcopy('swap', edits={t: e_qref(lines_all(swap_members('dev', 'var8'))) for t in GSET[:3]})
    gate('dev_var8_swapped_all_sets', dd, gspec(), 1,
         'control var8 not killed by gate acc_pTT', check=lambda r, o: ('gate acc_pTT: device' in o and 'FAIL' in o, 'device and var8 both caught'), noauth=True)
    gate('dev_var8_swapped_all_sets_caught_by_reprocessing', dd, gspec(), 2, RP)
    dd = gcopy('P2b', edits={'c_def_dev1_q20': e_qref(lines_all(copy_member('dev', 'nvh')))})
    gate('P2b_ctl_dev_is_a_correct_result', dd, gspec(), 1, 'control tau_ns0 not killed by gate acc_pTT', noauth=True)
    gate('P2b_ctl_dev_is_a_correct_result_caught_by_reprocessing', dd, gspec(), 2, RP)
    dd = gcopy('v8nk', edits={t: e_qref(lines_all(copy_member('var8', 'nvh'))) for t in GSET[:3]})
    gate('control_var8_not_killed', dd, gspec(), 1, 'control var8 not killed by gate acc_pTT', noauth=True)
    gate('control_var8_not_killed_caught_by_reprocessing', dd, gspec(), 2, RP)
    # P3b: failing evidence with a doctored aggregate (device mean -0.5) lying next to it: the gate never reads an aggregate
    dd = gcopy('P3b', edits=up)
    doc = json.loads(json.dumps(A))
    doc['pooled']['dev']['S_ens']['mean'] = -0.5
    for t in GSET[:3]:
        doc['runs'][t]['dev']['S_ens']['mean'] = -0.5
    for d in dd:
        json.dump(doc, open(os.path.join(d, 'agg.json'), 'w'))
    gate('P3b_doctored_aggregate_beside_failing_runs', dd, gspec(), 1, 'gate acc_pTT: device', noauth=True)
    gate('P3b_doctored_aggregate_beside_failing_runs_caught_by_reprocessing', dd, gspec(), 2, RP)
    # P3c: (a) the shuffled aggregate cannot be an input: beside the true runs the decision and the per-run means are the true ones
    dd = gcopy('P3c_a')
    doc = json.loads(json.dumps(A))
    doc['runs']['m_def_dev2']['dev']['S_ens']['mean'] = 5.0
    for d in dd:
        json.dump(doc, open(os.path.join(d, 'agg.json'), 'w'))
    gate('P3c_shuffled_aggregate_beside_true_runs', dd, gspec(), 4, 'GATE TEST-ONLY',
         check=lambda r, o: (max(abs(a - b) for a, b in zip(r['gates'][0]['device_run_means'], [A['runs'][t]['dev']['S_ens']['mean'] for t in GSET[:3]])) <= 1e-12,
                             'per-run device means are the true ones'))
    # (b) the evidence analogue: one run's device errors x1e5, another's x1e-5, pooled mean unchanged.  pool 'pooled' cannot see it (a known
    # limit of a pooled criterion: rc 4, recorded as such, not a detection); pool 'each_run' catches it (FAIL)
    ud = {'m_def_dev2': e_qref(lines_all(scale('dev', 5, 1e5))), 'm_def_dev3': e_qref(lines_all(scale('dev', 5, 1e-5)))}
    dd = gcopy('P3c_b', edits=ud)
    gate('P3c_evidence_up_down_pooled_LIMIT_not_detected', dd, gspec(), 4, 'GATE TEST-ONLY', noauth=True,
         check=lambda r, o: (abs(r['gates'][0]['device_value'] - rep0['gates'][0]['device_value']) < 1e-9 and max(r['gates'][0]['device_run_means']) > 4,
                             'pooled unchanged %+.4f, run means %s' % (r['gates'][0]['device_value'], ' '.join('%+.2f' % v for v in r['gates'][0]['device_run_means']))))
    gate('P3c_evidence_up_down_each_run', dd, smut(lambda s: s['gates'][0].__setitem__('pool', 'each_run')), 1, 'gate acc_pTT: device', noauth=True)
    gate('P3c_evidence_up_down_caught_by_reprocessing', dd, gspec(), 2, RP)
    # fix round 5 (M1): hostin is pinned by the spec.  One run processed with hostin 1 among hostin-0 runs: refused (was R_E_hostin_1_ok,
    # rc 4, before hostin was a spec key).  hostin 1 stays a valid argument (no din members without ALF_GPU_TAUDIAG): every run with
    # hostin 1 and a spec with hostin 1 is accepted, the re-processing runs with hostin 1
    gate('M1_one_run_hostin_1_spec_hostin_0', gcopy('hostin1', edits={D3: e_args(hostin=1)}), gspec(), 2, 'processed with hostin 1, the spec says 0')
    gate('R_E_hostin_1_all_runs_spec_hostin_1_ok', gcopy('hostin1all', edits={t: e_args(hostin=1) for t in GSET}), gspec(hostin=1), 4, 'GATE TEST-ONLY',
         check=lambda r, o: (all(a['source'] == 'program' and a['equal'] == a['numbers'] > 0 for a in r['authenticity']), 're-processed with hostin 1, all equal'))
    # fix round 5 (N4): the program runs in a minimal env; LD_PRELOAD in the gate's env is refused (CLI, the variable set for that call only)
    sp0 = os.path.join(gdir(), 'spec_baseline_gate_test_only_kept_subset.json')
    rp = os.path.join(gdir(), 'rep_N4_ld_preload_refused.json')
    rc, o = sh([PY, '-I', os.path.join(HERE, 't207_d5_gate.py'), '--spec', sp0, '--spec-sha', sha_file(sp0), '--prog', qbin('def'), '--report', rp] + B,
               env=dict(os.environ, LD_PRELOAD='/nonexistent/libx.so'))
    rj = json.load(open(rp)) if os.path.exists(rp) else {}
    record('N4_ld_preload_refused', 'rc 2: LD_PRELOAD is set', rc == 2 and 'LD_PRELOAD is set' in o and rj.get('exit_code') == 2 and 'GATE FAIL' not in o, ' '.join([l for l in o.split('\n') if l.startswith('GATE ')][-1:]))
    pe = GATE.program_env()
    allowed = {'PATH', 'HOME', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS'}
    LAST['rc'] = None
    record('N4_program_env_minimal', 'only PATH HOME OMP OPENBLAS', set(pe) <= allowed and 'LD_LIBRARY_PATH' not in pe and pe['OPENBLAS_NUM_THREADS'] == '1',
           'program env keys %s' % sorted(pe))
    fix4_tests(B, A, rep0)


# ---------------------------------------------------------------- fix round 4 (review 207-2d)
NUMT = re.compile(r'^(-?)([0-9])\.([0-9]{4})E([+-][0-9]{3})$')


def bump(tok, k):
    """an ES12.4E3 token changed by k units in its last printed digit (no carry across the exponent: asserted)"""
    m = NUMT.match(tok)
    v = int(m.group(2) + m.group(3)) + k
    assert 10000 <= v <= 99999, tok
    s = '%05d' % v
    return '%s%s.%sE%s' % (m.group(1), s[0], s[1:], m.group(4))


def bump_num(member, idx, k):
    def g(line):
        segs = line.split(' | ')
        for j, s in enumerate(segs):
            t = s.split()
            if j and t[0] == member:
                t[1 + idx] = bump(t[1 + idx], k)
                segs[j] = ' '.join(t)
        return ' | '.join(segs)
    return g


def scale_all(member, f):
    """every non-negative number of member times f (the -1 placeholders kept), as the reviewer's probes do"""
    def g(line):
        segs = line.split(' | ')
        for j, s in enumerate(segs):
            t = s.split()
            if j and t[0] == member:
                segs[j] = ' '.join([t[0]] + [fmt(float(x) * f) if float(x) >= 0 else x for x in t[1:]])
        return ' | '.join(segs)
    return g


def hdr_of(line):
    h = line.split(' | ')[0].split()
    return {h[j]: h[j + 1] for j in range(1, len(h) - 1, 2)}


def lines_where(fn, pred):
    """fn applied to the QREF data lines whose header satisfies pred(header dict)"""
    def g(t):
        return '\n'.join(fn(l) if (l.startswith('QREF file ') and pred(hdr_of(l))) else l for l in t.split('\n'))
    return g


def indep_points(dirs, member='dev', refs=None):
    """INDEPENDENT of the gate and of t207_d5_agg.py (read from the qref.txt text here): per point with B (run index, (nst, nfe, br), S_ens of
    member for pTT = log10 pTT(member) - mean over the references other than member of log10 pTT(ref))"""
    refs = refs or QP.ref_names('def', 'nvh')
    out = []
    for i, d in enumerate(dirs):
        for ln in open(os.path.join(d, 'qref.txt')).read().split('\n'):
            if not ln.startswith('QREF file '):
                continue
            hd = hdr_of(ln)
            if hd['hasb'] != '1':
                continue
            v = {s.split()[0]: float(s.split()[6]) for s in ln.split(' | ')[1:] if len(s.split()) == 9}
            r = [k for k in refs if k != member]
            out.append((i, (int(hd['nst']), int(hd['nfe']), int(hd['br'])), math.log10(v[member]) - sum(math.log10(v[k]) for k in r) / len(r)))
    return out


def indep_q(v, p):
    s = sorted(v)
    k = (len(s) - 1) * p
    i = int(math.floor(k))
    j = min(i + 1, len(s) - 1)
    return s[i] + (s[j] - s[i]) * (k - i)


def indep_strata(pts):
    st = {}
    for i, p, v in pts:
        st.setdefault(p, []).append(v)
    return {p: (len(v), sum(v) / len(v)) for p, v in st.items()}


def e_binary(full, stamp=True):
    """the run's ALF binary sha256 set to FULL in provenance.json (and in the binary line of stamp.txt unless stamp is False), records consistent"""
    def e(d):
        pj = os.path.join(d, 'provenance.json')
        pv = json.load(open(pj))
        pv['binary_sha256'] = full
        json.dump(pv, open(pj, 'w'), indent=1)
        if stamp:
            rewrite(os.path.join(d, 'stamp.txt'), lambda t: re.sub(r'^binary [0-9a-f]{64}', 'binary ' + full, t, flags=re.M))
        fix_records(d)
    return e


def fabricate_cover1(d, ri):
    """the reviewer's F1 (rev207-2d p3.py) on a copied run directory: qref.txt invented for every dump of the run (cover 1), the dumps
    replaced by short text files, manifest, processing arguments (cover 1, keep 0), qref_sha256 and manifest_sha256 all made consistent"""
    tmpl = [ln for ln in open(os.path.join(d, 'qref.txt')).read().split('\n') if ln.startswith('QREF file ')]
    tseg = [{sg.split()[0]: sg for sg in ln.split(' | ')[1:]} for ln in tmpl]
    st0 = json.load(open(os.path.join(d, 'qref_spec.json')))
    a = st0['args']
    sp = QP.make_spec(d, a['path'], a['qs'], a['nul'], a['hostin'], 1, 0)
    json.dump(sp, open(os.path.join(d, 'qref_spec.json'), 'w'), indent=1)
    qd = os.path.join(d, 'qref')
    for f in os.listdir(qd):
        os.remove(os.path.join(qd, f))
    out, man = [], []
    for (k, p, nst, nfe, br) in QP.cover_dumps(sp):
        hasb = 1 if nst < sp['run']['nstm'] else 0
        src = tseg[k % len(tseg)]
        f = 1 + 1e-3 * (k + 1000 * ri)
        sg = []
        for nm in QP.expected_names(sp, k, hasb, 'nvh'):
            if nm in ('qself', 'qbr'):
                sg.append(src[nm])
                continue
            tok = [fmt(float(x) * f) for x in src[nm].split()[1:]]
            if not hasb:
                tok[5:8] = ['-1.0000E+000'] * 3
            sg.append(nm + ' ' + ' '.join(tok))
        out.append('QREF file %d pass %d nst %d nfe %d br %d hbr %d n %d ver 2 hasb %d dg 0 tq 1.80 | ' % (k, p, nst, nfe, br, br, sp['run']['n'], hasb) + ' | '.join(sg))
        fn = 'qref/q%06d.bin' % k
        open(os.path.join(d, fn), 'wb').write(('fake %d\n' % k).encode())
        man.append('%s  %s' % (sha_file(os.path.join(d, fn)), fn))
    out.append('QREFEND dumps %d ok %d err 0' % (len(out), len(out)))
    open(os.path.join(d, 'qref.txt'), 'w').write('\n'.join(out) + '\n')
    open(os.path.join(d, 'qref_manifest.txt'), 'w').write('\n'.join(sorted(man)) + '\n')
    pj = os.path.join(d, 'qref_proc.json')
    x = json.load(open(pj))
    x['args']['cover'] = 1
    x['args']['keep'] = 0
    json.dump(x, open(pj, 'w'), indent=1)
    fix_qsha(d)
    fix_msha(d)


def fix4_tests(B, A, rep0):
    """the reviewer's probes D5, E1, S1/S2, F1, L1 and the new checks of fix round 4: each planted defect must give INPUT_ERROR or FAIL as
    named; every new check is also shown silent on correct input (the baselines, or the cases marked 'silent')"""
    D3 = 'm_def_dev3'
    M = GSET[:3]
    SET4 = ['m_def_dev2', 'm_def_dev3', 'm_def_dev4', 'c_def_dev1_q20']
    # ---- 6. spec integrity: only the pre-registered spec (its sha256) is evaluated
    s030 = smut(lambda s: s['gates'][0].__setitem__('bound', 0.30))
    gate('F6_spec_edited_after_preregistration', B, s030, 2, 'differs from --spec-sha', cli=True,
         spec_sha=hashlib.sha256(json.dumps(gspec(), indent=1).encode()).hexdigest())
    gate('F6_spec_sha_of_another_file', B, gspec(), 2, 'differs from --spec-sha', spec_sha='0' * 64)
    # ---- 2. pinning: seeds (S1/S2 cherry-picking), ALF binary (old binary), parameters (easier model)
    S = gcopy('S', edits={t: e_qref(lines_all(scale_all('dev', 10 ** 0.14))) for t in M})
    gate('S1_x1.38_cherry_pick_best_run', [S[0], S[3]], gspec(n_runs={'main': 1, 'ctl_tau_ns0': 1}), 2, 'the spec requires exactly', cli=True)
    gate('S2_x1.38_cherry_pick_two_runs', [S[0], S[2], S[3]], gspec(n_runs={'main': 2, 'ctl_tau_ns0': 1}), 2, 'the spec requires exactly', cli=True)
    sd = gspec()['seeds']
    gate('S3_spec_seed_set_replaced', B, gspec(seeds={'main': sd['main'][:2] + ['a' * 64], 'ctl_tau_ns0': sd['ctl_tau_ns0']}), 2, 'the spec requires exactly')
    gate('S4_spec_main_seed_set_shared_with_dev1_ctl_run', B, gspec(seeds={'main': sd['main'][:2] + [seeds_of('m_def_dev1')], 'ctl_tau_ns0': sd['ctl_tau_ns0']}), 2,
         'is listed for role')
    gate('S5_ctl_seed_equals_main_seed_in_spec', B, gspec(seeds={'main': sd['main'], 'ctl_tau_ns0': [sd['main'][0]]}), 2, 'is listed for role')
    gate('S6_seeds_role_missing', B, gspec(seeds={'main': sd['main']}), 2, 'must be a dictionary with exactly the roles')
    gate('S7_seeds_not_hex', B, gspec(seeds={'main': sd['main'][:2] + ['xyz'], 'ctl_tau_ns0': sd['ctl_tau_ns0']}), 2, 'distinct seeds_sha256')
    gate('S8_seeds_duplicate', B, gspec(seeds={'main': sd['main'][:2] + sd['main'][:1], 'ctl_tau_ns0': sd['ctl_tau_ns0']}), 2, 'distinct seeds_sha256')
    gate('S9_no_seeds_key', B, smut(lambda s: s.pop('seeds')), 2, "missing keys ['seeds']")
    gate('B1_spec_binary_other', B, gspec(binary_sha256='0' * 64), 2, 'ALF binary sha256')
    gate('B2_spec_binary_not_hex', B, gspec(binary_sha256='abc'), 2, 'binary_sha256 \'abc\' is not a sha256')
    gate('B3_one_run_old_binary', gcopy('B3', edits={D3: e_binary(prov_of('m_log_dev2', 'binary_sha256'))}), gspec(), 2, 'ALF binary sha256')
    gate('B4_provenance_binary_not_stamp', gcopy('B4', edits={D3: e_binary('b' * 64, stamp=False)}), gspec(), 2, 'binary line of stamp.txt')
    gate('B5_no_binary_key', B, smut(lambda s: s.pop('binary_sha256')), 2, "missing keys ['binary_sha256']")
    gate('P1_spec_parameters_other', B, gspec(parameters_sha256='0' * 64), 2, 'parameters_sha256')
    easier = chain(e_file(lambda d: open(os.path.join(d, 'parameters'), 'a').write('! an easier model\n')), e_file(fix_records))
    gate('P2_easier_model_all_runs_parameters_changed', gcopy('P2easy', edits={t: easier for t in SET4}), gspec(), 2, 'the spec requires parameters_sha256')
    # ---- 3. main env restricted, control faults from the fixed table
    E = gcopy('E1', edits={t: e_env(lambda e: dict(e, ALF_GPU_TAUSWAP='255')) for t in SET4})
    gate('E1_main_env_TAUSWAP_255', E, smut(lambda s: s['env_roles']['main'].__setitem__('ALF_GPU_TAUSWAP', '255')), 2, 'env_roles.main', cli=True)
    em = lambda k, v: smut(lambda s: s['env_roles']['main'].__setitem__(k, v))        # noqa: E731
    gate('E2_main_QVAR_16', B, em('ALF_GPU_TAUDIAG_QVAR', '16'), 2, 'env_roles.main')
    gate('E3_main_QHOST_0', B, em('ALF_GPU_TAUDIAG_QHOST', '0'), 2, 'env_roles.main')
    gate('E4_main_without_DEBUG', B, smut(lambda s: s['env_roles']['main'].pop('ALF_GPU_DEBUG')), 2, 'env_roles.main')
    gate('E5_main_FROM_not_int', B, em('ALF_GPU_TAUDIAG_FROM', '33x'), 2, 'ALF_GPU_TAUDIAG_FROM')
    ce = lambda extra, omit=('ALF_GPU_TAUDIAG_QVAR',), **kw: gspec(env_roles={'main': gspec()['env_roles']['main'], 'ctl_tau_ns0': {'extra_env': extra, 'omit_env': list(omit)}}, **kw)  # noqa: E731
    gate('E6_ctl_TAUSWAP_255', B, ce({'ALF_GPU_TAUSWAP': '255'}), 2, 'is not a fault of the fixed list')
    gate('E7_ctl_TAU_NS_1', B, ce({'ALF_GPU_TAU_NS': '1'}), 2, 'is not a fault of the fixed list')
    gate('E8_ctl_TAUSWAP_8_on_def_path', B, ce({'ALF_GPU_TAUSWAP': '8'}), 2, 'is not a fault of the fixed list for the def path')
    gate('E9_ctl_TAU_NS_0_on_log_path', B, ce({'ALF_GPU_TAU_NS': '0'}, path='log', refs=QP.ref_names('log', 'nvh'), nul=list(NUL_LOG), prog_sha=prog_of('m_log_dev2'),
                                               members={'main': mlist(NUL_LOG), 'ctl_tau_ns0': mlist(NUL_LOG, var=False)}), 2, 'is not a fault of the fixed list for the log path')
    gate('E10_ctl_omits_QHOST', B, ce({'ALF_GPU_TAU_NS': '0'}, omit=('ALF_GPU_TAUDIAG_QVAR', 'ALF_GPU_TAUDIAG_QHOST')), 2, 'a control may omit only')
    gate('E11_ctl_adds_no_fault', B, ce({}), 2, 'adds no fault')
    for i, v in enumerate(('perm', 'eps:abc', 'eps:-1', 'eps:0', 'eps:1e300', 'eps:1e-9@0', 'eps:1e-9@1:0', 'eps:1e-9@a', 'eps:1e-9@1:2:3',
                           'eps:2.5e-10@2', 'eps:1e-9@3:2')):     # 207-5: the last two (the dropped '@' form) moved here from E13
        gate('E12_ctl_FAULT_TAU_invalid_%d' % i, B, ce({'ALF_GPU_TAU_NS': '0', 'ALF_GPU_FAULT_TAU': v}), 2, 'is not a fault of the fixed list')
    # silent: valid fault values pass the spec check; the dev5 control run has no such env, so its role is not found (env allow-list)
    for i, v in enumerate(('eps:1e-9', 'eps:.5d-9')):      # 207-5: the '@' forms are no longer valid (see E12_9, E12_10)
        gate('E13_ctl_FAULT_TAU_valid_silent_%d' % i, B, ce({'ALF_GPU_TAU_NS': '0', 'ALF_GPU_FAULT_TAU': v}), 2, 'env allow-list',
             check=lambda r, o: ('fixed list' not in o, 'spec accepted, refused at the run'))
    # 207-5 (manager, 2026-11-01): 207-3's localized fault is the separate variable ALF_GPU_FAULT_TAU_AT=<nst>:<nfe>, only with ALF_GPU_FAULT_TAU
    for i, (ex, om) in enumerate((({'ALF_GPU_FAULT_TAU': 'eps:1e-9'}, ()), ({'ALF_GPU_FAULT_TAU': 'eps:1e-9', 'ALF_GPU_FAULT_TAU_AT': '2:1'}, ()),
                                  ({'ALF_GPU_FAULT_TAU': 'eps:3e-10', 'ALF_GPU_FAULT_TAU_AT': '12:2'}, ('ALF_GPU_TAUDIAG_QVAR',)))):
        gate('E14_ctl_FAULT_TAU_with_or_without_AT_valid_silent_%d' % i, B, ce(ex, omit=om), 2, 'env allow-list',
             check=lambda r, o: ('fixed list' not in o and 'without ALF_GPU_FAULT_TAU' not in o, 'spec accepted, refused at the run'))
    gate('E15_ctl_FAULT_TAU_AT_without_FAULT_TAU', B, ce({'ALF_GPU_TAU_NS': '0', 'ALF_GPU_FAULT_TAU_AT': '2:1'}), 2, 'ALF_GPU_FAULT_TAU_AT without ALF_GPU_FAULT_TAU')
    gate('E15b_ctl_only_FAULT_TAU_AT', B, ce({'ALF_GPU_FAULT_TAU_AT': '2:1'}), 2, 'ALF_GPU_FAULT_TAU_AT without ALF_GPU_FAULT_TAU')
    for i, v in enumerate(('0:1', '2:0', '2', 'a:1', '1:2:3', '2:1 ', '-1:1')):
        gate('E16_ctl_FAULT_TAU_AT_invalid_%d' % i, B, ce({'ALF_GPU_FAULT_TAU': 'eps:1e-9', 'ALF_GPU_FAULT_TAU_AT': v}), 2,
             'is not a fault of the fixed list' if v.strip() == v else 'value of ALF_GPU_FAULT_TAU_AT')
    gate('E17_main_with_FAULT_TAU', B, em('ALF_GPU_FAULT_TAU', 'eps:1e-9'), 2, 'env_roles.main')
    # ---- 4. refs frozen (D5: refs = three noisiest variants let the swapped device pass)
    Dd = gcopy('D', edits={t: e_qref(lines_all(swap_members('dev', 'var8'))) for t in M})
    nov8 = lambda s: s.__setitem__('controls', [c for c in s['controls'] if c['name'] != 'var8'])       # noqa: E731
    d5 = smut(lambda s: (nov8(s), s.__setitem__('refs', ['n5', 'n100004', 'n100005']), s.pop('test_only')))
    gate('D5_refs_noisy_nulls_not_frozen', Dd, d5, 2, 'differ from the frozen reference list', cli=True)
    d5t = smut(lambda s: (nov8(s), s.__setitem__('refs', ['n5', 'n100004', 'n100005'])))
    gate('D5_refs_noisy_nulls_test_only_LIMIT', Dd, d5t, 4, 'GATE TEST-ONLY', noauth=True,
         check=lambda r, o: (r['verdict'] == 'TEST-ONLY', 'another refs list is allowed only under test_only, which never gives PASS'))
    gate('D5_refs_noisy_nulls_test_only_caught_by_reprocessing', Dd, d5t, 2, '(re-processing check)')
    gate('D4_refs_u1_u2_only', B, smut(lambda s: (nov8(s), s.__setitem__('refs', ['u1', 'u2']))), 2, "refs: 'u1' is not a reference member")
    gate('D2_refs_plus_u1_u2', B, smut(lambda s: s['refs'].extend(['u1', 'u2'])), 2, "refs: 'u1' is not a reference member")
    gate('H2_refs_plus_var8', B, smut(lambda s: (nov8(s), s['refs'].append('var8'))), 2, "refs: 'var8' is not a reference member")
    gate('D6_refs_reordered_without_test_only', B, smut(lambda s: (s.pop('test_only'), s.__setitem__('refs', s['refs'][::-1]))), 2, 'differ from the frozen reference list')
    gate('D7_refs_log_list_on_def_without_test_only', B, smut(lambda s: (s.pop('test_only'), s.__setitem__('refs', QP.ref_names('log', 'nvh')))), 2,
         'differ from the frozen reference list')
    gate('D8_frozen_refs_silent_without_test_only', B, smut(lambda s: s.pop('test_only')), 2, 'must cover every dump',
         check=lambda r, o: ('frozen' not in o, 'frozen list accepted; refused only for cover 40'))
    # ---- 1. authenticity
    F = gcopy('F1', edits={t: (lambda i: e_file(lambda d: fabricate_cover1(d, i)))(i) for i, t in enumerate(SET4)})
    f1 = smut(lambda s: (s.__setitem__('cover', 1), s.pop('test_only'), s['expect_points'].__setitem__('pTT', 210)))
    gate('F1_fabricated_cover1_text_dumps', F, f1, 2, 'not a dump', cli=True)
    for d in F:                                        # 4 x 280 small fake files: the record of the test stays, the files go
        shutil.rmtree(os.path.join(d, 'qref'))
    kept = sorted(os.listdir(os.path.join(base(D3), 'qref')))
    gate('F1c_text_file_posing_as_kept_dump', gcopy('F1c', edits={D3: e_file(lambda d: replace_dump(d, kept[2], b'fake dump\n'))}), gspec(), 2, 'not a dump')
    gate('F1c2_zero_bytes_of_dump_size', gcopy('F1c2', edits={D3: e_file(lambda d: replace_dump(d, kept[2], b'\0' * len(dump_bytes(d, kept[2]))))}), gspec(), 2,
         'bad dump header')

    def patch_hdr(i, val):
        def f(d):
            b = bytearray(dump_bytes(d, kept[2]))
            b[4 * i:4 * i + 4] = struct.pack('<i', val)
            replace_dump(d, kept[2], bytes(b))
        return e_file(f)
    h3 = struct.unpack('<13i', dump_bytes(base(D3), kept[2])[:52])
    gate('F1d_dump_header_pass_patched', gcopy('F1d', edits={D3: patch_hdr(3, h3[3] + 1)}), gspec(), 2, 'header pass')
    gate('F1d2_dump_header_N_patched', gcopy('F1d2', edits={D3: patch_hdr(7, h3[7] + 1)}), gspec(), 2, 'header N')
    gate('F1d3_dump_padded', gcopy('F1d3', edits={D3: e_file(lambda d: replace_dump(d, kept[2], dump_bytes(d, kept[2]) + b'\0' * 8))}), gspec(), 2,
         'the header implies')

    def other_position(d):                             # a genuine dump of the run at another place, header patched to this place
        b = bytearray(dump_bytes(d, kept[3]))
        b[8:28] = dump_bytes(d, kept[2])[8:28]         # file, pass, nst, nfe, br of dump kept[2]
        replace_dump(d, kept[2], bytes(b))
    gate('F1b_real_dump_from_another_position_header_patched', gcopy('F1b', edits={D3: e_file(other_position)}), gspec(), 2, '(re-processing check)')

    def dev_changed(d):                                # the device output in a kept dump changed (one entry x (1 + 1e-6)); qref.txt kept
        b = bytearray(dump_bytes(d, kept[2]))
        h, n, off = dump_layout(os.path.realpath(os.path.join(base(D3), 'qref', kept[2])), False)
        x = struct.unpack('<d', bytes(b[off['dev']:off['dev'] + 8]))[0]
        b[off['dev']:off['dev'] + 8] = struct.pack('<d', x * (1 + 1e-6))
        replace_dump(d, kept[2], bytes(b))
    gate('F1e_dump_dev_output_changed_qref_kept', gcopy('F1e', edits={D3: e_file(dev_changed)}), gspec(), 2, '(re-processing check)')
    gate('F1f_qref_number_2_digits_off', gcopy('F1f', edits={D3: e_qref(lines_all(bump_num('u1', 0, 2), only=3))}), gspec(), 2,
         'exact equality required', cli=True)
    gate('F1f2_qref_number_2_digits_down', gcopy('F1f2', edits={D3: e_qref(lines_all(bump_num('dev', 5, -2), only=5))}), gspec(), 2,
         'exact equality required')
    # fix round 5 (N1): exact equality.  CHANGED ASSERTION: was F1g_qref_number_1_digit_off_silent (rc 4, a one-digit difference allowed);
    # the pinned program is bit-identical (reviewer probe P1), so a one-digit difference is now INPUT_ERROR
    gate('F1g_qref_number_1_digit_off', gcopy('F1g', edits={D3: e_qref(lines_all(bump_num('u1', 0, 1), only=3))}), gspec(), 2, 'exact equality required')
    gate('F1g2_qref_number_1_digit_down_dev', gcopy('F1g2', edits={D3: e_qref(lines_all(bump_num('dev', 5, -1), only=5))}), gspec(), 2, 'exact equality required')
    gate('F1h_qref_tq_changed_silent', gcopy('F1h', edits={D3: e_qref(lines_all(lambda l: re.sub(r' tq +[0-9.]+ \|', ' tq 99.99 |', l), only=2))}), gspec(), 4,
         'GATE TEST-ONLY')
    gate('F1i_prog_other_binary', B, gspec(), 2, 'differs from the required prog_sha', prog=qbin('log'))
    gate('F1j_prog_missing', B, gspec(), 2, 'not an executable file', prog=os.path.join(gdir(), 'no_such_program'))
    pc = os.path.join(gdir(), 'progcopy')
    os.makedirs(pc)
    shutil.copy2(qbin('def'), os.path.join(pc, 'qref_pp.out'))
    gate('F1k_prog_copy_same_hash_silent', B, gspec(), 4, 'GATE TEST-ONLY', prog=os.path.join(pc, 'qref_pp.out'))
    gate('F1l_manifest_changed_proc_sha_stale', gcopy('F1l', edits={D3: e_file(lambda d: rewrite(os.path.join(d, 'qref_manifest.txt'),
                                                                                                     lambda t: '\n'.join(t.split('\n')[:-1][::-1]) + '\n'))}),
         gspec(), 2, 'differs from manifest_sha256')
    gate('F1m_proc_without_manifest_sha256', gcopy('F1m', edits={D3: e_proc(lambda x: x.pop('manifest_sha256'))}), gspec(), 2, 'no valid manifest_sha256')
    # ---- 5. localized defects: L1 (x5 at the 4th kept dump of every main run) and the new statistics
    L1 = gcopy('L1', edits={t: e_qref(lines_all(scale_all('dev', 5.0), only=4)) for t in M})
    gate('L1_x5_one_point_in_7_pooled_LIMIT', L1, gspec(), 4, 'GATE TEST-ONLY', noauth=True,
         check=lambda r, o: (r['gates'][0]['device_value'] <= BOUND, 'pooled mean %+.4f <= %.2f: a pooled mean misses it' % (r['gates'][0]['device_value'], BOUND)))
    gate('L1_x5_one_point_in_7_caught_by_reprocessing', L1, gspec(), 2, '(re-processing check)', cli=True)
    # the new statistics on the baseline: every number equals the independent computation from the qref.txt text (indep_points)
    P = indep_points(B[:3])
    pv = [v for i, p, v in P]
    Pst = indep_strata(P)
    gs = [dict(name='acc_pTT', measure='pTT', statistic='S_ens', pool='pooled', bound=BOUND),
          dict(name='strat', measure='pTT', statistic='S_ens', pool='by_stratum', bound=5.0),
          dict(name='mx', measure='pTT', statistic='S_ens_max', pool='pooled', bound=5.0),
          dict(name='q95', measure='pTT', statistic='S_ens_q95', pool='pooled', bound=5.0),
          dict(name='mxrun', measure='pTT', statistic='S_ens_max', pool='each_run', bound=5.0)]

    def stats_match(r, o):
        g = {x['name']: x for x in r['gates']}
        exp = dict(acc_pTT=sum(pv) / len(pv), strat=max(t[1] for t in Pst.values()), mx=max(pv), q95=indep_q(pv, 0.95),
                   mxrun=max(max(v for i, p, v in P if i == j) for j in range(3)))
        dm = max(abs(g[k]['device_value'] - v) for k, v in exp.items())
        cnt = sorted((s['nst'], s['nfe'], s['br'], s['n']) for s in g['strat']['device_strata']) == sorted(p + (t[0],) for p, t in Pst.items())
        loc = r['localized']['pTT'][0]
        lok = loc['member'] == 'dev' and abs(loc['pooled']['max'] - max(pv)) < 1e-12 and len(loc['strata']) == len(Pst)
        return dm <= 1e-12 and cnt and lok, 'max |gate - independent| %.1e; strata %s; localized section %s' % (
            dm, ' '.join('(%d,%d,%d):%d:%+.3f' % (p + t) for p, t in sorted(Pst.items())), 'ok' if lok else 'WRONG')
    gate('L0_new_statistics_match_independent_baseline', B, gspec(gates=gs), 3, 'gates without a control that trips them', check=stats_match)
    # by_stratum at the self-test bound 0.15: the expected outcome on the baseline from the independent numbers (dev strata, controls' strata)
    sp_bs = smut(lambda s: s['gates'][0].__setitem__('pool', 'by_stratum'))
    v8 = indep_strata(indep_points(B[:3], member='var8'))
    cdv = indep_strata(indep_points([B[3]]))
    exp_bs = 4 if (max(t[1] for t in Pst.values()) <= BOUND and max(t[1] for t in v8.values()) > BOUND and max(t[1] for t in cdv.values()) > BOUND) else 1
    gate('L2_by_stratum_baseline', B, sp_bs, exp_bs, GW[exp_bs], check=lambda r, o: (True, 'expected from the independent strata: rc %d' % exp_bs))
    # The kept dev5 subset has ONE stratum: every kept dump is at nst 1, nfe 1, br 1 (stride 2 x cover 40 = 80 = 0 mod NSTM*N_FL = 8: the
    # bias the keep rule refuses for thinning), so by_stratum equals pooled on it (L2 above, L8 below) and cannot be shown end to end here.
    # Unit-level test of the stratum code instead (GATE.values / GATE.decide, the functions the gate runs): the real per-point values of the
    # baseline runs (GATE.load_run on B), with positions assigned as a cover-1 run would visit them (nst 1..3 x nfe 1..2 x br 1..2 cycling over
    # the points of each run).  Bound of this test only: the baseline's largest stratum mean + 0.35 (not a calibration).
    record_unit_by_stratum(B)
    # tail statistics: bounds of this test only = baseline device value + 0.35 (half of log10 5; NOT a calibration, 207-3 sets bounds)
    vmax, vq = max(pv), indep_q(pv, 0.95)
    tg = [dict(name='acc_pTT', measure='pTT', statistic='S_ens', pool='pooled', bound=BOUND),
          dict(name='mx', measure='pTT', statistic='S_ens_max', pool='pooled', bound=round(vmax + 0.35, 6)),
          dict(name='q95', measure='pTT', statistic='S_ens_q95', pool='pooled', bound=round(vq + 0.35, 6))]
    gate('L4_tail_bounds_baseline_silent', B, gspec(gates=tg), 3, 'gates without a control that trips them',
         check=lambda r, o: (all(g['passed'] for g in r['gates']), 'mx %+.4f <= %+.4f, q95 %+.4f <= %+.4f' % (vmax, vmax + 0.35, vq, vq + 0.35)))
    order = sorted(range(len(P)), key=lambda j: -P[j][2])
    top = [(P[j][0], P[j][1]) for j in order[:2]]

    def plant_top(nmax):
        """x5 (pTT of dev) at the nmax points with the largest baseline values (located by run and dump number)"""
        eds = {}
        for j in order[:nmax]:
            ri = P[j][0]
            # the j-th point of run ri in the order of its qref.txt: the dump number of that line
            k = [hdr_of(l)['file'] for l in open(os.path.join(B[ri], 'qref.txt')).read().split('\n') if l.startswith('QREF file ') and hdr_of(l)['hasb'] == '1'][
                [jj for jj in range(len(P)) if P[jj][0] == ri].index(j)]
            eds.setdefault(M[ri], []).append(k)
        return {t: e_qref(lines_where(scale('dev', 5, 5.0), lambda h, ks=tuple(ks): h['file'] in ks)) for t, ks in eds.items()}
    Lm = gcopy('Lmax', edits=plant_top(1))
    gate('L5_x5_at_the_worst_point_max_FAIL', Lm, gspec(gates=tg), 1, 'gate mx: device', noauth=True,
         check=lambda r, o: ({g['name']: g['passed'] for g in r['gates']}['acc_pTT'], 'pooled mean still passes'))
    Lq = gcopy('Lq95', edits=plant_top(2))
    gate('L6_x5_at_the_two_worst_points_q95_FAIL', Lq, gspec(gates=tg), 1, 'gate q95: device', noauth=True,
         check=lambda r, o: ({g['name']: g['passed'] for g in r['gates']}['acc_pTT'], 'pooled mean still passes'))
    gate('L5_caught_by_reprocessing', Lm, gspec(gates=tg), 2, '(re-processing check)')
    # the reviewer's L1 under the new statistics: the expected outcome from the independent numbers of the L1 data
    PL = indep_points(L1[:3])
    plv = [v for i, p, v in PL]
    e_mx = 1 if max(plv) > vmax + 0.35 else 3
    e_q = 1 if indep_q(plv, 0.95) > vq + 0.35 else 3
    e_l1 = 1 if 1 in (e_mx, e_q) else 3
    gate('L7_reviewer_L1_under_tail_statistics', L1, gspec(gates=tg), e_l1, GW[e_l1], noauth=True,
         check=lambda r, o: (True, 'independent: max %+.4f (bound %+.4f), q95 %+.4f (bound %+.4f)' % (max(plv), vmax + 0.35, indep_q(plv, 0.95), vq + 0.35)))
    plst = indep_strata(PL)
    e_bs = 1 if max(t[1] for t in plst.values()) > BOUND else 4
    gate('L8_reviewer_L1_under_by_stratum_0.15', L1, sp_bs, e_bs if exp_bs == 4 else 1, GW[e_bs if exp_bs == 4 else 1], noauth=True,
         check=lambda r, o: (True, 'independent strata %s' % ' '.join('(%d,%d,%d):%d:%+.3f' % (p + t) for p, t in sorted(plst.items()))))


def record_unit_by_stratum(B):
    s, table = GATE.check_spec(gspec())
    POS = [(1 + j % 3, 1 + (j // 3) % 2, 1 + (j // 6) % 2) for j in range(7)]      # 7 points per run -> 7 positions of the 12

    def runs_with_pos(plant=None, f=0.0):
        runs = [GATE.load_run(d, s, table) for d in B]
        for x in runs:
            for m, d in x['res']['pTT'].items():
                d['pos'] = [POS[j % len(POS)] for j in range(d['n'])]
                if plant is not None and m == 'dev' and x['role'] == 'main':
                    d['re'] = [v + (f if p == plant else 0.0) for p, v in zip(d['pos'], d['re'])]
        return runs
    base_runs = runs_with_pos()
    main = [x for x in base_runs if x['role'] == 'main']
    # independent: per-stratum means of the dev values with the same positions
    st = {}
    for x in main:
        for p, v in zip(x['res']['pTT']['dev']['pos'], x['res']['pTT']['dev']['re']):
            st.setdefault(p, []).append(v)
    ind = {p: (len(v), sum(v) / len(v)) for p, v in st.items()}
    smax = max(ind, key=lambda p: ind[p][1])
    b = round(ind[smax][1] + 0.35, 6)
    g = dict(name='strat', measure='pTT', statistic='S_ens', pool='by_stratum', bound=b)
    pooled, per, strata, allv = GATE.values(main, 'dev', g)
    dm = max(abs(strata[p][1] - ind[p][1]) for p in ind)
    ok = dm <= 1e-12 and set(strata) == set(ind) and all(strata[p][0] == ind[p][0] for p in ind) and GATE.judge(pooled, per, strata, g, False) <= b
    LAST['rc'] = None
    record('L3_unit_by_stratum_baseline_silent', 'strata = independent, <= bound', ok, '%d strata (points %s), largest %+.4f at %s, bound %+.4f, max |gate - independent| %.1e' % (
        len(ind), '/'.join(str(ind[p][0]) for p in sorted(ind)), ind[smax][1], smax, b, dm))
    # planted: x5 (+log10 5) on the dev values of the points of the stratum with the largest mean, in every main run
    f = math.log10(5.0)
    pr = runs_with_pos(plant=smax, f=f)
    sp = dict(s, gates=[g], controls=[dict(name='var8', kind='variant', member='var8', trips='strat')], n_runs={'main': 3, 'ctl_tau_ns0': 1})
    verdict, reasons, lines, rg, rc_ = GATE.decide(pr, sp)
    pooled_new = statistics_mean([v for x in pr if x['role'] == 'main' for v in x['res']['pTT']['dev']['re']])
    record('L3_unit_by_stratum_planted_FAIL', 'FAIL: gate strat device', verdict == 'FAIL' and any('gate strat: device' in r for r in reasons),
           '%s: %s; pooled mean after planting %+.4f (bound 0.15 for a pooled gate: %s)' % (verdict, '; '.join(reasons)[:120], pooled_new,
                                                                                           'misses it' if pooled_new <= BOUND else 'also catches it'))
    # the same planted runs under pool 'pooled' with the same bound: not detected (the limit by_stratum removes)
    gp = dict(g, pool='pooled')
    v2, r2, l2, rg2, rc2 = GATE.decide(pr, dict(sp, gates=[gp], controls=[dict(name='var8', kind='variant', member='var8', trips='strat')]))
    record('L3_unit_pooled_same_bound_LIMIT', 'device passes pooled', rg2[0]['passed'], 'pooled device %+.4f <= bound %+.4f' % (rg2[0]['device_value'], b))


def statistics_mean(v):
    return sum(v) / len(v)


def reproc_tests():
    """the gate's re-processing check (GATE.authenticate, the function the gate runs) on every kept dev5 dump of every repro run directory:
    every number of the re-computed lines exactly equal to qref.txt (fix round 5)"""
    for tag in sorted(os.listdir(os.path.join(W, 'runs'))):
        d = base(tag)
        if not os.path.isfile(os.path.join(d, 'qref.txt')):
            continue
        v1 = tag.startswith('s_')
        try:
            rspec = QP.load_spec(d)
            prog = qbin(rspec['args']['path'])
            run = QP.load(os.path.join(d, 'qref.txt'), allow_v1=v1)
            man, msha = QP.verify_dumps(d, rspec)
            QP.proc_record(d, rspec, prog=sha_file(prog), qref_sha256=run.qref_sha256, manifest_sha256=msha)
            nh = QP.check_dump_headers(d, rspec, allow_v1=v1)
            x = dict(dir=d, run=run, lines=run.lines, kept=QP.kept_dumps(rspec), manifest=man)
            au = GATE.authenticate([x], {'prog_sha': sha_file(prog)}, prog)[0]
            ok = au['dumps'] == len(x['kept']) == nh and au['numbers'] > 0 and au['equal'] == au['numbers']
            det = '%d dumps (headers checked %d), %d numbers, %d equal' % (au['dumps'], nh, au['numbers'], au['equal'])
        except Exception as e:
            ok, det = False, '%s: %s' % (type(e).__name__, str(e)[:200])
        LAST['rc'] = 0 if ok else 2
        record('reproc_' + tag, 'all numbers equal', ok, det)


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
    # fix round 5 (M1), the reviewer's probe in small: a thinned run (s_def, cover 1, KEEP 3: 6 of 16 dumps kept for re-processing) against a
    # spec that pins keep and min_reprocessed.  Refused for keep (spec 0) and for too few kept dumps (spec keep 3, min_reprocessed 7); both
    # checks come before the version 1 refusal, which ends the silent case (spec keep 3, min_reprocessed 6)
    os.makedirs(gdir(), exist_ok=True)
    v1env = dict(gspec()['env_roles']['main'], ALF_GPU_TAUDIAG_QREF='1', ALF_GPU_TAUDIAG_FROM='399')
    sv1 = lambda **kw: gspec(qs=1, cover=1, env_roles={'main': v1env}, n_runs={'main': 1}, controls=[], members={'main': mlist(NUL_DEF)},     # noqa: E731
                             seeds={'main': [seeds_of('s_def')]}, binary_sha256=prov_of('s_def', 'binary_sha256'), **kw)
    gate('M1_thinned_run_keep_3_spec_keep_0', [dd], sv1(keep=0, min_reprocessed=6), 2, 'processed with keep 3, the spec says 0')
    gate('M1_thinned_run_6_kept_min_reprocessed_7', [dd], sv1(keep=3, min_reprocessed=7), 2, '6 kept dumps to re-process, the spec requires at least min_reprocessed 7')
    gate('M1_thinned_run_6_kept_min_reprocessed_6_silent', [dd], sv1(keep=3, min_reprocessed=6), 2, 'only with allow_v1',
         check=lambda r, o: ('min_reprocessed' not in o and 'processed with' not in o, 'keep and min_reprocessed accepted; refused later as version 1'))
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
    # fix round 2: the binary's .src.sha256 must name the source next to the script (prog_src is a verified claim, not a copied one)
    w2 = os.path.join(ST, 'w2'); os.makedirs(os.path.join(w2, 'bin')); os.makedirs(os.path.join(w2, 'tmp'))
    for f in ('qref_pp.out', 'qref_pp.out.sha256'):
        os.symlink(os.path.join(W, 'bin', f), os.path.join(w2, 'bin', f)) if f.endswith('.out') else shutil.copy2(os.path.join(W, 'bin', f), os.path.join(w2, 'bin', f))
    with open(os.path.join(w2, 'bin', 'qref_pp.src.sha256'), 'w') as f:
        f.write('0' * 64 + '  t207_d5_qref.F90\n')
    d, files = mkdumpdir('wrong_src_sha', 'def')
    env = dict(os.environ, D5W=w2, D5RUNS=os.path.join(ST, 'dumpruns'), DUMP_COVER='40', TMPDIR=os.path.join(w2, 'tmp'))
    rc, out = sh(['bash', os.path.join(HERE, 't207_d5_proc.sh'), d, 'def', '0', NU, '0', '0'], env=env)
    record('dump_src_sha256_not_this_source', 'rc 1: refusing, nothing written', rc == 1 and 'built from this source' in out and not os.path.exists(os.path.join(d, 'qref.txt'))
           and not os.path.exists(os.path.join(d, 'qref_spec.json')), out.strip().split('\n')[-1][:100])
    # fix round 5 (N2): end-to-end device FAIL through the real processing.  The device output (every entry of the dev block) of every kept
    # dump of the three main runs scaled by (1 + 1e-6), t207_d5_proc.sh run on the copies with the arguments of the repro runs (qs 20, the
    # default null codes, hostin 0, cover 40, keep 0), the gate CLI on them plus the unchanged control run: FAIL (device above the bound)
    # with the real re-processing step (which agrees with the new qref.txt)
    n2dirs = []
    n2ok = True
    for tag in GSET[:3]:
        d, files = mkdumpdir('N2_' + tag, 'def', src_tag=tag)
        for f in files:
            fp = os.path.join(d, 'qref', f)
            h, n, off = dump_layout(fp, False)
            with open(fp, 'r+b') as g:
                g.seek(off['dev'])
                blk = g.read(64 * n * n)
                v = struct.unpack('<%dd' % (8 * n * n), blk)
                g.seek(off['dev'])
                g.write(struct.pack('<%dd' % (8 * n * n), *[x * (1 + 1e-6) for x in v]))
        rc, out = proc(d, 'def', 0, qs='20', nul=','.join(str(c) for c in NUL_DEF))
        n2ok = n2ok and rc == 0
        n2dirs.append(d)
    record('N2_proc_on_scaled_device_dumps', 'rc 0 for the 3 runs', n2ok, 'processing of the scaled dumps: %s' % ('ok' if n2ok else last(out)))
    os.makedirs(gdir(), exist_ok=True)
    gate('N2_device_scaled_1e-6_every_kept_dump_end_to_end_FAIL', n2dirs + [base('c_def_dev1_q20')], gspec(), 1, 'gate acc_pTT: device', cli=True,
         check=lambda r, o: (all(a['source'] == 'program' and a['equal'] == a['numbers'] > 0 for a in r['authenticity']) and r['gates'][0]['device_value'] > BOUND,
                             'device %+.4f > %.2f, re-processing all equal' % (r['gates'][0]['device_value'], BOUND)))
    for d in n2dirs:
        shutil.rmtree(os.path.join(d, 'qref'))          # own dump copies (21 x 6 MB); the records of the test stay
    # the LOG path baseline (one run through the program with the LOG build); its proc record carries the program hash
    d, files = mkdumpdir('log_baseline', 'log')
    rc, out = proc(d, 'log', 0)
    pr = json.load(open(os.path.join(d, 'qref_proc.json')))
    want = open(os.path.join(W, 'bin', 'qref_log.out.sha256')).read().split()[0]
    record('dump_log_baseline', 'rc 0, prog hash recorded', rc == 0 and pr['prog'] == want, 'rc %d, prog %s' % (rc, pr['prog'][:16]))
    # 207-5: the interpreter of the parser steps is D5_PYTHON (default ~/alf-gpu/venv/bin/python3), recorded in qref_proc.json with the machine
    import socket
    mc = pr.get('machine') or {}
    record('Z1_proc_records_python_and_machine', 'default python, this host', pr.get('python') == os.path.expanduser('~/alf-gpu/venv/bin/python3')
           and mc.get('host') == socket.gethostname() and re.match(r'^[0-9a-f]{64}$', mc.get('openblas_sha256') or '') is not None,
           'python %s, host %s, openblas %s' % (pr.get('python'), mc.get('host'), (mc.get('openblas') or '')[-40:]))
    d, files = mkdumpdir('python_missing', 'def')
    env = dict(os.environ, D5W=W, D5RUNS=os.path.join(ST, 'dumpruns'), DUMP_COVER='40', TMPDIR=os.path.join(W, 'tmp'), D5_PYTHON='/nonexistent/python3')
    rc, out = sh(['bash', os.path.join(HERE, 't207_d5_proc.sh'), d, 'def', '0', NU, '0', '0'], env=env)
    record('Z2_proc_python_not_executable_refused', 'rc 1, nothing written', rc == 1 and 'D5_PYTHON' in out and not os.path.exists(os.path.join(d, 'qref_spec.json'))
           and files_left(d) == files, out.strip().split('\n')[-1][:100])
    wrap = os.path.join(ST, 'pywrap.sh')
    with open(wrap, 'w') as f:
        f.write('#!/bin/bash\necho "$@" >> %s.calls\nexec %s "$@"\n' % (wrap, os.path.expanduser('~/alf-gpu/venv/bin/python3')))
    os.chmod(wrap, 0o755)
    d, files = mkdumpdir('python_wrapper', 'def')
    env = dict(env, D5_PYTHON=wrap)
    rc, out = sh(['bash', os.path.join(HERE, 't207_d5_proc.sh'), d, 'def', '0', NU, '0', '0'], env=env)
    calls = open(wrap + '.calls').read().count('t207_d5_qparse.py') if os.path.exists(wrap + '.calls') else 0
    pw = json.load(open(os.path.join(d, 'qref_proc.json'))) if os.path.exists(os.path.join(d, 'qref_proc.json')) else {}
    record('Z3_proc_python_from_D5_PYTHON', 'rc 0, wrapper used and recorded', rc == 0 and calls >= 2 and pw.get('python') == wrap,
           'rc %d, %d parser calls through the wrapper, recorded %s' % (rc, calls, pw.get('python')))


def main():
    if os.path.exists(ST):
        sys.exit('refusing: %s exists' % ST)
    os.makedirs(ST)
    if WHAT not in ('text', 'gate', 'dump', 'reproc', 'all'):
        sys.exit('unknown mode %s' % WHAT)
    if WHAT in ('text', 'gate', 'all'):
        gate_tests()
    if WHAT in ('text', 'all'):
        text_tests()
        agg_tests()
        pairs_tests()
    if WHAT in ('reproc', 'all'):
        reproc_tests()
    if WHAT in ('dump', 'all'):
        dump_tests()
    bad = [r for r in RESULTS if not r[2]]
    print('\nself-test: %d tests, %d with the expected outcome, %d WRONG' % (len(RESULTS), len(RESULTS) - len(bad), len(bad)))
    with open(os.path.join(ST, 'results.txt'), 'w') as f:
        f.write('name\texpected\toutcome\tlast_rc\tdetail\n')
        for name, expect, ok, detail, rc in RESULTS:
            f.write('%s\t%s\t%s\t%s\t%s\n' % (name, expect, 'ok' if ok else 'WRONG', rc, detail.replace('\n', ' ')))
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
