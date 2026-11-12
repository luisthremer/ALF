#!/usr/bin/env python3
"""P7 report tables (stdout, markdown).  Inputs: data/ (see NOTES.md)."""
import sys, re, bisect, collections, glob
sys.path.insert(0, '.')
exec(open('tables.py').read().split("out = []")[0])    # helpers: shim, off, sym, site_phase, kphase, kinfo, parse_*
def agg(rows):
    a = collections.OrderedDict()
    for ph, nm, c, d, l in rows:
        k = (ph, nm); x = a.setdefault(k, [0, 0, l]); x[0] += c; x[1] += d
        if l and not x[2]: x[2] = l
    return [(k[0], k[1], v[0], v[1], v[2]) for k, v in a.items()]
d1 = D + 'j1/'
pp = diff(parse_pp(d1 + 'V7_pp3/rank0.log'), parse_pp(d1 + 'V7_pp1/rank0.log'))
blk3 = parse_blk(D + 'j4/V7_of3/rank0.log'); blk1 = parse_blk(D + 'j4/V7_of1/rank0.log')
cs3 = off(D + 'j6/V7_cs3/rank0.log'); cs1 = off(D + 'j6/V7_cs1/rank0.log')
sh3 = shim(d1 + 'V7_pf3/rank0.log'); sh1 = shim(d1 + 'V7_pf1/rank0.log')
rows = []
for k, (c, dev, host) in pp.items():
    if k.startswith('Data'): continue
    rows.append((kphase(k), 'omp ' + k.replace('__', ' ').replace('.', ':'), c, dev / 1e3, kinfo.get(k, '')))
for k in blk3:
    n3, d3, h3 = blk3[k]; n1, d1_, h1 = blk1.get(k, (0, 0, 0))
    ls = {'gather': 'global 2x32x32 = 2048 WI (~64 sub-groups)', 'scale': 'global 2xbx576, mean 36,864 WI', 'fgather': 'global kx576, mean 6,045 WI', 'copy': 'DMA copy (mean 13.7 kB)', 'gemm m576': 'oneMKL 576x576xk, k=2..23 (mostly 10-14)'}.get(k, '')
    rows.append(('U', 'sycl ' + k, (n3 - n1) / 2, (d3 - d1_) / 2, ls))
for k in cs3:
    n3, d3, h3 = cs3[k]; n1, d1_, h1 = cs1.get(k, (0, 0, 0))
    shape, site = k.rsplit(' @', 1); sn = '?'
    for a_ in site.split(','):
        nm_ = sym(int(a_, 16))
        if 'alf_la' not in nm_: sn = nm_; break
    ph = site_phase(sn)
    rows.append((ph, 'mkl ' + shape, (n3 - n1) / 2, (d3 - d1_) / 2, 'oneMKL kernel(s), full tile'))
qr = {}
for kk in ('geqrf_dev', 'unmqr_dev', 'getrf'):
    qr[kk] = ((sh3[kk][0] - sh1[kk][0]) / 2, (sh3[kk][1] - sh1[kk][1]) / 2)
rows = agg(rows)
rows_q = [('S+C', 'devqr geqrf call (sync; panel 1 work-group 60%, apply 40%)', qr['geqrf_dev'][0], qr['geqrf_dev'][1], 'panel: 1 WG x 512 WI; apply: ~70 WG x 256 WI'),
          ('S+C', 'devqr unmqr call (sync)', qr['unmqr_dev'][0], qr['unmqr_dev'][1], '~72 WG x 256 WI'),
          ('C', 'getrf oneMKL call (hybrid host/device, sync)', qr['getrf'][0], qr['getrf'][1], 'oneMKL, partly host')]
allr = sorted(rows + rows_q, key=lambda r: -r[3])
print('### Top 20 kernels / ops, V7, one rank on tile 0, per sweep (C; jobs 560831 pp, 560845 off, 560847 callsites)\n')
print('| # | phase | kernel / op | calls/sweep | device s/sweep | mean us | launch size |'); print('|--|--|--|--|--|--|--|')
for i, r in enumerate(allr[:20]):
    print('| %d | %s | %s | %.0f | %.3f | %.1f | %s |' % (i + 1, r[0], r[1], r[2], r[3], 1e6 * r[3] / r[2] if r[2] else 0, r[4]))
print('\n### Per phase, top 6 each\n')
for ph in ('U', 'S', 'C', 'O'):
    print('**phase %s**\n' % ph); print('| kernel / op | calls/sweep | dev s/sweep | mean us | launch |'); print('|--|--|--|--|--|')
    sel = [r for r in allr if r[0] == ph or (ph in r[0] and '+' in r[0] and ph in ('S', 'C'))]
    for r in sel[:6]: print('| %s | %.0f | %.3f | %.1f | %s |' % (r[1], r[2], r[3], 1e6 * r[3] / r[2] if r[2] else 0, r[4]))
    print()
tot = collections.defaultdict(float)
for r in rows: tot[r[0]] += r[3]
print('device time from event timestamps per phase [s/sweep] (omp kernels + sycl block ops + oneMKL BLAS, QR/getrf calls excluded):', {k: round(v, 3) for k, v in tot.items()})
print('SC (shared stab/cgr kernels) split ~50/50 by call count.')
# node layout: kernel time inflation
def node(v, ranks=range(32)):
    acc = collections.defaultdict(lambda: [0, 0])
    for r in ranks:
        x = diff(parse_pp(D + 'j2/%s_pp3/rank%d.log' % (v, r)), parse_pp(D + 'j2/%s_pp1/rank%d.log' % (v, r)))
        for k, (c, dv, h) in x.items(): acc[k][0] += c / 32; acc[k][1] += dv / 32
    return acc
nd = node('V7')
print('\n### Node layout (32 ranks, 4/tile, CCS 4), V7: omp kernels, mean over the 32 ranks vs one rank alone (C: job 560841 vs 560831)\n')
print('| kernel | calls/sweep | alone ms/sweep | 4/tile ms/sweep | slow-down | launch |'); print('|--|--|--|--|--|--|')
lst = sorted(((k, v) for k, v in nd.items() if not k.startswith('Data')), key=lambda kv: -kv[1][1])[:12]
for k, (c, dv) in lst:
    a = pp.get(k, (0, 0, 0))
    print('| %s | %.0f | %.1f | %.1f | %.2fx | %s |' % (k, c, a[1], dv, dv / a[1] if a[1] else 0, kinfo.get(k, '')))
tn = collections.defaultdict(float)
for k, (c, dv) in nd.items():
    if not k.startswith('Data'): tn[kphase(k)] += dv / 1e3
print('\nomp-kernel device time per phase and sweep, per rank at 4/tile:', {k: round(v, 3) for k, v in tn.items()})
