#!/usr/bin/env python3
# Ticket 207 dev3 (diagnosis of the default-path paired spread, development only): calibration of the paired statistic against the
# size of the perturbation.  Reads the ALF_GPU_TAUDIAG_AMP chains (" | A" GTT precision numbers, " | Ad" max rel. differences of
# the new GTT from the host chain's new GTT; last Ad value: the device's) of the "TAUDIAG pt " lines.
#   usage: t207_d3_ampstat.py run.log FROM_PASS AMP1,AMP2,...   (amplitudes in the order of ALF_GPU_TAUDIAG_AMP, in pairs)
# Prints per chain: median / 90 % quantile of Ad, sd of log10(A_k/host) and its ratio to the reference sd(host+1ulp(b)/host+1ulp(a));
# per amplitude pair: sd of log10(A_b/A_a) and the ratio; and the device row (dev/host, Ad of the device).  Fail-closed (exit 2).
import sys, math, statistics as st
fn, frm, amps = sys.argv[1], int(sys.argv[2]), [float(a) for a in sys.argv[3].split(',')]
na = len(amps)
def die(msg):
    print('FAIL:', msg); sys.exit(2)
rows = []
for l in open(fn):
    if not l.lstrip().startswith('TAUDIAG pt '):
        continue
    f = l.split()
    try:
        ps = int(f[f.index('pass') + 1])
        i = f.index('dev'); j = f.index('hostfull'); k = f.index('P'); a = f.index('A'); d = f.index('Ad')
        dev = [float(x) for x in f[i + 1:i + 5]]; host = [float(x) for x in f[j + 1:j + 5]]; P = [float(x) for x in f[k + 1:k + 6]]
        A = [float(x) for x in f[a + 1:a + 1 + na]]; Ad = [float(x) for x in f[d + 1:d + 2 + na]]
    except (ValueError, IndexError):
        die('malformed TAUDIAG line: ' + l[:80])
    if len(A) != na or len(Ad) != na + 1:
        die('A / Ad block length')
    rows.append(dict(ps=ps, dgtt=dev[1], hgtt=host[1], P=P, A=A, Ad=Ad))
used = [r for r in rows if r['ps'] >= frm]
if len(used) < 2:
    die('only %d rows from pass %d' % (len(used), frm))
def sd(a, b):
    x = []
    for r in used:
        u, v = a(r), b(r)
        if not (u > 0 and v > 0):
            die('non-positive value in pass %d' % r['ps'])
        x.append(math.log10(u / v))
    return st.mean(x), st.pstdev(x)
def q(v, p):
    v = sorted(v); return v[min(len(v) - 1, int(p * len(v)))]
_, sref = sd(lambda r: r['P'][2], lambda r: r['P'][1])
print('file %s  used (pass >= %d) %d  reference sd(host+1ulp(b)/host+1ulp(a)) %.4f' % (fn, frm, len(used), sref))
print('  %-30s %10s %10s %8s %7s %7s' % ('chain', 'Ad median', 'Ad q90', 'mean', 'sd', 'ratio'))
m, s = sd(lambda r: r['P'][1], lambda r: r['hgtt'])
print('  %-30s %10s %10s %+8.4f %7.4f %7.3f' % ('host+1ulp(a) / host (P1)', '-', '-', m, s, s / sref))
for k in range(na):
    m, s = sd(lambda r, k=k: r['A'][k], lambda r: r['hgtt'])
    ad = [r['Ad'][k] for r in used]
    print('  %-30s %10.2e %10.2e %+8.4f %7.4f %7.3f' % ('amp %g [%d] / host' % (amps[k], k), st.median(ad), q(ad, 0.9), m, s, s / sref))
m, s = sd(lambda r: r['dgtt'], lambda r: r['hgtt'])
ad = [r['Ad'][na] for r in used]
print('  %-30s %10.2e %10.2e %+8.4f %7.4f %7.3f' % ('dev / host', st.median(ad), q(ad, 0.9), m, s, s / sref))
print('  pairs (same amplitude, independent perturbations):')
for k in range(0, na - 1, 2):
    if amps[k] != amps[k + 1]:
        continue
    m, s = sd(lambda r, k=k: r['A'][k + 1], lambda r, k=k: r['A'][k])
    print('  %-30s %10s %10s %+8.4f %7.4f %7.3f' % ('amp %g (b)/(a)' % amps[k], '', '', m, s, s / sref))
