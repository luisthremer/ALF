#!/usr/bin/env python3
"""P7: build the per-kernel table, per-phase busy/idle split and round-trip counts from data/ (see NOTES.md).  usage: tables.py > data/tables.md"""
import sys, re, bisect, collections
sys.path.insert(0, '.')
from ana import *
D = 'data/'
def shim(path):
    o = {}
    for l in open(path, errors='replace'):
        m = re.match(r"ALFSH_TIMING (\S+)\s+calls\s+(\d+)\s+total\s+([\d.]+) s", l)
        if m: o[m.group(1)] = (int(m.group(2)), float(m.group(3)))
    return o
def off(path): return parse_tag(path, r"ALFP7_OFF (.+?)\s+calls\s+(\d+)\s+dev_s\s+([\d.]+)\s+host_s\s+([\d.]+)")
# --- symbol table for the call sites
syms = []
for l in open(D + 'ALF.out.nm.txt'):
    p = l.split()
    if len(p) >= 3: syms.append((int(p[0], 16), p[2]))
syms.sort(); addrs = [a for a, _ in syms]
def sym(a):
    i = bisect.bisect_right(addrs, a) - 1
    return syms[i][1] if i >= 0 else '?'
def site_phase(sname):
    s = sname.lower()
    if 'dev_decompose' in s or 'dev_testscale' in s or 'wrapul' in s or 'wrapur' in s: return 'S'
    if 'cgr' in s or 'dev_cgr' in s or 'cgr_solve' in s or 'dev_prec' in s or 'dev_orth' in s: return 'C'
    if 'gpu_obs' in s or 'obs' in s: return 'O'
    if 'gpu_wrap' in s or 'gpu_lowrank' in s or 'gpu_update' in s or 'wrapgr' in s or 'gpu_lr' in s: return 'U'
    return '?'
def kphase(k):
    if k.startswith('gpu_wrap'): return 'U'
    if k.startswith('gpu_obs'): return 'O'
    cg = ('k_cgr_phase', 'k_prec', 'k_det_phase', 'k_tpup', 'k_ct', 'gpu_stab_cgr', 'om_getrf', 'k_rowdiv', 'k_colinv', 'k_iperm')
    st = ('k_phase', 'k_scale_cols', 'k_scale_row1', 'k_scale_col1', 'k_triu', 'k_ns', 'k_testscale', 'k_perm_rows', 'k_perm_cols')
    if any(x in k for x in cg): return 'C'
    if any(x in k for x in st): return 'S'
    return 'SC'
kinfo = {}
for l in open(D + 'kinfo_V7.txt'):
    m = re.match(r"KINFO\s+(\d+) (\S+) \| (.*)", l)
    if m: kinfo[short(m.group(2))] = m.group(3).split('  ')[0]
out = []
P = lambda *a: out.append(' '.join(str(x) for x in a))
for lay, d, ranks, v in (('1 rank / tile 0', D + 'j1/', [0], 'V7'),):
    pp = diff(parse_pp(d + v + '_pp3/rank0.log'), parse_pp(d + v + '_pp1/rank0.log'))
    blk3 = parse_blk(D + 'j4/V7_of3/rank0.log'); blk1 = parse_blk(D + 'j4/V7_of1/rank0.log')
    of3 = off(D + 'j4/V7_of3/rank0.log'); of1 = off(D + 'j4/V7_of1/rank0.log')
    cs3 = off(D + 'j6/V7_cs3/rank0.log'); cs1 = off(D + 'j6/V7_cs1/rank0.log')
    sh3 = shim(D + 'j1/V7_pf3/rank0.log'); sh1 = shim(D + 'j1/V7_pf1/rank0.log')
    rows = []   # (phase, name, calls/sweep, dev s/sweep, launch)
    for k, (c, dev, host) in pp.items():
        if k.startswith('Data'): continue
        rows.append((kphase(k), 'omp ' + k, c, dev / 1e3, kinfo.get(k, '')))
    for k in blk3:
        n3, d3, h3 = blk3[k]; n1, d1, h1 = blk1.get(k, (0, 0, 0))
        rows.append(('U', 'sycl ' + k, (n3 - n1) / 2, (d3 - d1) / 2, ''))
    # dispatched BLAS by call site
    ph_blas = collections.defaultdict(lambda: [0.0, 0.0])
    for k in cs3:
        n3, d3, h3 = cs3[k]; n1, d1, h1 = cs1.get(k, (0, 0, 0))
        shape, site = k.rsplit(' @', 1)
        sn = '?'
        for a_ in site.split(','):
            nm_ = sym(int(a_, 16))
            if 'alf_la' not in nm_: sn = nm_; break
        ph = site_phase(sn)
        rows.append((ph, 'mkl %s [%s]' % (shape, sn[:34]), (n3 - n1) / 2, (d3 - d1) / 2, ''))
    for r in rows: ph_blas[r[0]][0] += r[3]
    for kk, nm in (('geqrf_dev', 'devqr geqrf (sync call)'), ('unmqr_dev', 'devqr unmqr (sync call)'), ('getrf', 'getrf oneMKL (hybrid, sync call)')):
        c = (sh3[kk][0] - sh1[kk][0]) / 2; t = (sh3[kk][1] - sh1[kk][1]) / 2
        rows.append(('SC' if kk != 'getrf' else 'C', nm, c, t, ''))
    rows.sort(key=lambda r: -r[3])
    P('## Per-kernel table, V7, 1 rank on tile 0, per sweep (C: jobs 560831 / 560845 / 560846)')
    P('| # | phase | kernel / op | calls/sweep | device s/sweep | mean us | launch (teams x threads) |'); P('|--|--|--|--|--|--|--|')
    for i, r in enumerate(rows[:30]):
        P('| %d | %s | %s | %.0f | %.4f | %.1f | %s |' % (i + 1, r[0], r[1], r[2], r[3], 1e6 * r[3] / r[2] if r[2] else 0, r[4]))
    # phase totals
    tot = collections.defaultdict(float)
    for r in rows:
        if r[1].startswith('devqr') or r[1].startswith('getrf'): continue
        tot[r[0]] += r[3]
    P(); P('device busy per phase from event timestamps (omp kernels + sycl block ops + dispatched BLAS), excluding QR/getrf sync calls:', {k: round(v, 3) for k, v in tot.items()})
print('\n'.join(out))
