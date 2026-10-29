#!/usr/bin/env python3
# Ticket 207 (LOG precision tail): paired statistics of the ALF_GPU_TAUDIAG=2 lines of a run.log.
#   usage: t207_diagstat.py run.log [FROM_PASS]
# Per point the GTT precision number (max |old - new|, the pair that carries the tail) of: the device pass (dev), the full host
# chain on the same trajectory (host = host TAU_M), and the host-propagated chains of diag_post (P0: device output; P1, P2: host
# CGR2_2 on the host udvr with U perturbed by ~1 ulp; P3, P4: the same on the device udvr).  Prints exceedance counts and the
# mean / sd / max of log10 ratios: if dev/host has the same spread as P2/P1 (two 1-ulp perturbations of the host chain), the
# device-host difference per event is rounding-level noise amplified by the conditioning, not a bias.
import sys, math, statistics as st
frm = int(sys.argv[2]) if len(sys.argv) > 2 else 0
rows = []
for l in open(sys.argv[1]):
    if 'TAUDIAG' not in l:
        continue
    f = l.split()
    ps = int(f[f.index('pass') + 1])
    i = f.index('dev'); j = f.index('hostfull')
    P = [float(x) for x in f[f.index('P') + 1:f.index('P') + 6]] if 'P' in f else [0.0] * 5
    rows.append(dict(ps=ps, nst=int(f[f.index('nst') + 1]), dev=max(map(float, f[i + 1:i + 5])), host=max(map(float, f[j + 1:j + 5])),
                     dgtt=float(f[i + 2]), hgtt=float(f[j + 2]), P=P))
vals = [v for r in rows for v in [r['dev'], r['host'], r['dgtt'], r['hgtt']] + r['P']]
if not rows or not all(math.isfinite(v) for v in vals):
    print('FAIL: no TAUDIAG lines or non-finite values (%d rows)' % len(rows))
    sys.exit(2)
print('all points: %d   max (all pairs): dev %.3e  host %.3e' % (len(rows), max(r['dev'] for r in rows), max(r['host'] for r in rows)))
for t in (1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 2e-6, 1e-5, 3e-5):
    print('  > %.0e: dev %5d  host %5d' % (t, sum(r['dev'] > t for r in rows), sum(r['host'] > t for r in rows)))
def lr(sel, a, b, name):
    x = [math.log10(a(r) / b(r)) for r in rows if sel(r) and a(r) > 0 and b(r) > 0]
    if len(x) > 1:
        print('  %-34s n %5d  mean %+.3f  sd %.3f  max|.| %.2f' % (name, len(x), st.mean(x), st.pstdev(x), max(abs(v) for v in x)))
print('log10 ratios of the per-point maxima (all pairs), all passes:')
lr(lambda r: True, lambda r: r['dev'], lambda r: r['host'], 'dev / host')
print('GTT pair, passes >= %d (perturbation chains active):' % frm)
sel = lambda r: r['ps'] >= frm and r['P'][1] > 0
lr(sel, lambda r: r['dgtt'], lambda r: r['hgtt'], 'dev / host')
lr(sel, lambda r: r['P'][1], lambda r: r['hgtt'], 'host+1ulp / host')
lr(sel, lambda r: r['P'][2], lambda r: r['P'][1], 'host+1ulp(b) / host+1ulp(a)')
lr(sel, lambda r: r['P'][3], lambda r: r['dgtt'], 'dev+1ulp / dev')
lr(sel, lambda r: r['P'][4], lambda r: r['P'][3], 'dev+1ulp(b) / dev+1ulp(a)')
lr(sel, lambda r: r['P'][0], lambda r: r['dgtt'], 'dev output host-propagated / dev')
