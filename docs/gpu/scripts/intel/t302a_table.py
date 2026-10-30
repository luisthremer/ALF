#!/usr/bin/env python3
"""Generate the per-routine GPU/host/hybrid table (acceptance 3) from the raw outputs of t302a_table.sh.
usage: t302a_table.py DIR   (DIR has table_timing.txt, table_timing_rep2.txt, table_timing_pinned1core.txt, table_trace.txt)
Method: (1) MKL_VERBOSE trace of one call per routine: oneapi::mkl::...GPU0 lines = call made on the device,
plain MKL_VERBOSE Z.../D... lines = host MKL call; (2) timing of the alf_la call on the device, default run (mean of two)
vs the same run pinned to ONE host core (taskset, 1 thread): a call whose time grows >= 1.5x at N=1024 or N=400 does
host work (hybrid); (3) the same operation as a direct host MKL call (host-pinned reference column)."""
import re, sys, collections
D = sys.argv[1] + '/'
def rd(f):
    d = {}
    for l in open(D + f, errors='replace'):
        if l.startswith('TABLE'):
            m = re.match(r'TABLE\s+(.+?)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)', l)
            d[(m[1], int(m[2]))] = (float(m[3]), float(m[5]))
    return d
a, b, p = rd('table_timing.txt'), rd('table_timing_rep2.txt'), rd('table_timing_pinned1core.txt')
tr = collections.defaultdict(lambda: [0, 0, []]); cur = None
for l in open(D + 'table_trace.txt', errors='replace'):
    if l.startswith('SECTION_END'): cur = None; continue
    if l.startswith('SECTION'): cur = re.match(r'SECTION N=(\d+) (.*)', l.strip()).groups(); continue
    if cur and l.startswith('MKL_VERBOSE') and 'oneMKL 2026' not in l and 'Detected GPU' not in l:
        if 'GPU0' in l and 'oneapi::mkl' in l: tr[cur][0] += 1; tr[cur][2].append(re.match(r'MKL_VERBOSE (\S+?)[\[\(]', l)[1])
        else: tr[cur][1] += 1; tr[cur][2].append('host:' + re.match(r'MKL_VERBOSE (\S+?)\(', l)[1])
print('%-24s %5s %9s %9s %9s %8s  %-22s %s' % ('routine', 'N', 'gpu_ms', 'host_ms', 'pin1c_ms', 'pin/gpu', 'trace (dev,host calls)', 'class'))
byr = collections.defaultdict(dict)
for (r, n) in a:
    g = (a[(r, n)][0] + b[(r, n)][0]) / 2
    byr[r][n] = (g, a[(r, n)][1], p[(r, n)][0], p[(r, n)][0] / g)
for r in byr:
    t = tr.get(('400', r)) or [0, 0, []]
    ratios = [byr[r][n][3] for n in (400, 1024)]
    if r == 'la_geqp3_host': cls = 'host (ZGEQP3 + transfers)'
    elif t[1] > 0: cls = 'host'
    elif t[0] == 0: cls = '?'
    elif max(ratios) >= 1.5: cls = 'hybrid (device call, host work: 1-core run slower)'
    else: cls = 'GPU'
    for n in (144, 400, 1024):
        g, h, pn, ra = byr[r][n]
        print('%-24s %5d %9.3f %9.3f %9.3f %8.2f  %-22s %s' % (r, n, g, h, pn, ra, '%d dev %s, %d host' % (t[0], ','.join(sorted(set(t[2]))[:2]), t[1]) if n == 400 else '', cls if n == 400 else ''))
