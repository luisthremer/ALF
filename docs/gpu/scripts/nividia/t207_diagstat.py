#!/usr/bin/env python3
# Ticket 207 (LOG precision tail): paired statistics of the ALF_GPU_TAUDIAG=2 lines of a run.log.
#   usage: t207_diagstat.py run.log FROM_PASS [--expect-points N] [--expect-used N] [--check MEAN SDRATIO K]
# Per point the GTT precision number (max |old - new|, the pair that carries the tail) of: the device pass (dev), the full host
# chain on the same trajectory (host = host TAU_M), and the host-propagated chains of diag_post (P0: device output; P1, P2: host
# CGR2_2 on the host udvr with U perturbed by ~1 ulp; P3, P4: the same on the device udvr).  Prints exceedance counts and the
# mean / sd / max of log10 ratios: if dev/host has the same spread as P2/P1 (two 1-ulp perturbations of the host chain), the
# device-host difference per event is rounding-level noise amplified by the conditioning, not a bias.
# Fail-closed (2026-10-29): exit 2 on no lines, a malformed line, a missing P block, a non-finite value, a non-positive value in
# a used row, an empty selection or a count different from --expect-*; with --check, exit 1 if a criterion fails (|mean| > MEAN,
# sd(dev/host)/sd(host+1ulp(b)/host+1ulp(a)) > SDRATIO, counts above 1e-7 / 3e-7 differ by more than K*sqrt(n_dev + n_host)).
import sys, math, statistics as st
args = sys.argv[1:]
if len(args) < 2:
    sys.exit(__doc__ if __doc__ else "usage: t207_diagstat.py run.log FROM_PASS [...]")
fn, frm = args[0], int(args[1])
opt = {a: args[i + 1:i + 4] for i, a in enumerate(args) if a.startswith('--')}
def die(msg):
    print('FAIL:', msg); sys.exit(2)
rows = []
for l in open(fn):
    if 'TAUDIAG' not in l:
        continue
    f = l.split()
    try:
        ps = int(f[f.index('pass') + 1])
        i = f.index('dev'); j = f.index('hostfull'); k = f.index('P')
        dev = [float(x) for x in f[i + 1:i + 5]]; host = [float(x) for x in f[j + 1:j + 5]]; P = [float(x) for x in f[k + 1:k + 6]]
    except (ValueError, IndexError):
        die('malformed TAUDIAG line: ' + l[:80])
    if len(dev) != 4 or len(host) != 4 or len(P) != 5 or not all(math.isfinite(v) for v in dev + host + P):
        die('non-finite or short TAUDIAG line: ' + l[:80])
    rows.append(dict(ps=ps, nst=int(f[f.index('nst') + 1]), dev=max(dev), host=max(host), dgtt=dev[1], hgtt=host[1], P=P))
if not rows:
    die('no TAUDIAG lines')
if '--expect-points' in opt and len(rows) != int(opt['--expect-points'][0]):
    die('%d TAUDIAG lines, expected %s (incomplete run?)' % (len(rows), opt['--expect-points'][0]))
print('all points: %d   max (all pairs): dev %.3e  host %.3e' % (len(rows), max(r['dev'] for r in rows), max(r['host'] for r in rows)))
for t in (1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 2e-6, 1e-5, 3e-5):
    print('  > %.0e: dev %5d  host %5d' % (t, sum(r['dev'] > t for r in rows), sum(r['host'] > t for r in rows)))
def lr(sel, a, b, name):
    used = [r for r in rows if sel(r)]
    if len(used) < 2:
        die('row "%s": only %d rows selected' % (name, len(used)))
    bad = [r for r in used if not (a(r) > 0 and b(r) > 0)]
    if bad:
        die('row "%s": %d selected rows with non-positive values (pass %d)' % (name, len(bad), bad[0]['ps']))
    x = [math.log10(a(r) / b(r)) for r in used]
    print('  %-34s n %5d  mean %+.4f  sd %.4f  max|.| %.2f' % (name, len(x), st.mean(x), st.pstdev(x), max(abs(v) for v in x)))
    return len(x), st.mean(x), st.pstdev(x)
print('log10 ratios of the per-point maxima (all pairs), all passes:')
lr(lambda r: True, lambda r: r['dev'], lambda r: r['host'], 'dev / host')
print('GTT pair, passes >= %d (perturbation chains active):' % frm)
sel = lambda r: r['ps'] >= frm
n, m, sd = lr(sel, lambda r: r['dgtt'], lambda r: r['hgtt'], 'dev / host')
lr(sel, lambda r: r['P'][1], lambda r: r['hgtt'], 'host+1ulp / host')
_, _, sref = lr(sel, lambda r: r['P'][2], lambda r: r['P'][1], 'host+1ulp(b) / host+1ulp(a)')
lr(sel, lambda r: r['P'][3], lambda r: r['dgtt'], 'dev+1ulp / dev')
lr(sel, lambda r: r['P'][4], lambda r: r['P'][3], 'dev+1ulp(b) / dev+1ulp(a)')
lr(sel, lambda r: r['P'][0], lambda r: r['dgtt'], 'dev output host-propagated / dev')
if '--expect-used' in opt and n != int(opt['--expect-used'][0]):
    die('%d rows used, expected %s' % (n, opt['--expect-used'][0]))
if '--check' in opt:
    mx, rx, kk = (float(v) for v in opt['--check'])
    fails = []
    if abs(m) > mx: fails.append('|mean| %.4f > %g' % (abs(m), mx))
    if sd / sref > rx: fails.append('sd ratio %.3f > %g' % (sd / sref, rx))
    for t in (1e-7, 3e-7):
        nd, nh = sum(r['dgtt'] > t for r in rows), sum(r['hgtt'] > t for r in rows)
        if abs(nd - nh) > kk * math.sqrt(nd + nh): fails.append('counts > %g: %d vs %d' % (t, nd, nh))
    print('CHECK:', 'FAIL ' + '; '.join(fails) if fails else 'PASS', '(sd ratio %.3f, mean %+.4f)' % (sd / sref, m))
    sys.exit(1 if fails else 0)
