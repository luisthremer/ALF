#!/usr/bin/env python3
"""Ticket 302c-2 round 2 (development only): summary of the T0 runs t302c-2b/runs/t0<RUN>_* generated from their run.log / info /
stamp.txt (no hand-made table).  Per run: rc, self-test verdict line (PASSED / FAILED), the largest self-test value with its tolerance
line ('right after CGR2_2', 'CGR2_2 alone ...'), refusal / error lines, Precision tau (mean, max) from info.
usage: t302c2b_t0eval.py RUN [WORKDIR]   (round 3: WORKDIR ~/alf-gpu/work/t302c-2c; default t302c-2b)"""
import glob, os, re, sys

W = os.path.expanduser(sys.argv[2] if len(sys.argv) > 2 else '~/alf-gpu/work/t302c-2b')
R = sys.argv[1]
num = re.compile(r'[0-9]\.[0-9]{3}E[+-][0-9]{2}')
for d in sorted(glob.glob(os.path.join(W, 'runs', 't0%s_*' % R, '*'))):
    rid = os.path.basename(os.path.dirname(d))[len('t0%s_' % R):]
    st = open(os.path.join(d, 'stamp.txt')).read() if os.path.exists(os.path.join(d, 'stamp.txt')) else ''
    m = re.search(r' rc=([0-9]+) wall=([0-9.]+)s lockwait=([0-9.]+)s .* env \[(.*)\]', st)
    rc, wall, lw, env = (m.group(1), m.group(2), m.group(3), m.group(4)) if m else ('?', '?', '?', '?')
    log = open(os.path.join(d, 'run.log'), errors='replace').read() if os.path.exists(os.path.join(d, 'run.log')) else ''
    verdict = sorted(set(re.findall(r'GPU tau self-test (PASSED|FAILED)', log)))
    stvals = []
    for line in log.splitlines():
        if 'GPU tau self-test right after CGR2_2' in line or 'GPU tau self-test CGR2_2 alone' in line:
            stvals += [float(x) for x in num.findall(line.split('(tol')[0].split(':')[-1])]
    errs = [l.strip()[:150] for l in log.splitlines() if re.search(r'refus|invalid|ERROR|error stop|not allowed|FAILED', l)][:2]
    prec = ''
    ip = os.path.join(d, 'info')
    if os.path.exists(ip):
        t = open(ip).read()
        p = re.search(r'Precision tau[^\n]*\n[^\n]*', t)
        prec = ' '.join(p.group(0).split()) if p else ''
    print('%-40s rc=%-3s wall=%-6s lockwait=%-5s %s%s%s%s  env[%s]' % (
        rid, rc, wall, lw, ('selftest ' + '/'.join(verdict) + ' ') if verdict else '', ('max %.3e ' % max(stvals)) if stvals else '',
        ('| ' + ' || '.join(errs) + ' ') if errs else '', ('| ' + prec) if prec else '', env))
