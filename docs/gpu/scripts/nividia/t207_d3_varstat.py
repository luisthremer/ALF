#!/usr/bin/env python3
# Ticket 207 dev3 (diagnosis of the default-path paired spread, development only): paired statistics of the variant chains of
# ALF_GPU_TAUDIAG_VAR (" | V" fields of the "TAUDIAG pt " lines) next to the production device chain.
#   usage: t207_d3_varstat.py run.log FROM_PASS MASK1,MASK2,...   (the masks in the order of ALF_GPU_TAUDIAG_VAR)
# Per row, over the points of passes >= FROM_PASS, GTT pair: n, mean, sd of log10(a/b) and the sd ratio to the reference
# sd(host+1ulp(b)/host+1ulp(a)) (= t207_diagstat.py's criterion b), plus the counts above 1e-7 and 3e-7 (all passes for dev/host,
# chain passes for the variants).  Rows: dev/host (production, device-propagated), P0/host (device output host-propagated),
# V_k/host (variant k, host-propagated), V_k/P0 (variant vs device output, both host-propagated).  Fail-closed: exit 2 on malformed
# input, a missing V block or a non-positive value.
import sys, math, statistics as st
fn, frm, masks = sys.argv[1], int(sys.argv[2]), [int(m) for m in sys.argv[3].split(',')]
nv = len(masks)
def die(msg):
    print('FAIL:', msg); sys.exit(2)
rows = []
for l in open(fn):
    if not l.lstrip().startswith('TAUDIAG pt '):
        continue
    f = l.split()
    try:
        ps = int(f[f.index('pass') + 1])
        i = f.index('dev'); j = f.index('hostfull'); k = f.index('P')
        dev = [float(x) for x in f[i + 1:i + 5]]; host = [float(x) for x in f[j + 1:j + 5]]; P = [float(x) for x in f[k + 1:k + 6]]
        V = [float(x) for x in f[f.index('V') + 1:f.index('V') + 1 + nv]] if 'V' in f else []
    except (ValueError, IndexError):
        die('malformed TAUDIAG line: ' + l[:80])
    if len(V) != nv:
        die('V block with %d values, expected %d' % (len(V), nv))
    rows.append(dict(ps=ps, dgtt=dev[1], hgtt=host[1], P=P, V=V))
if not rows:
    die('no TAUDIAG lines')
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
_, sref = sd(lambda r: r['P'][2], lambda r: r['P'][1])
print('file %s  points %d  used (pass >= %d) %d  reference sd(host+1ulp(b)/host+1ulp(a)) %.4f' % (fn, len(rows), frm, len(used), sref))
print('  %-28s %6s %8s %7s %7s %8s %8s' % ('row', 'n', 'mean', 'sd', 'ratio', '>1e-7', '>3e-7'))
def row(name, a, b, cnt):
    m, s = sd(a, b)
    c1 = sum(cnt(r) > 1e-7 for r in used); c3 = sum(cnt(r) > 3e-7 for r in used)
    print('  %-28s %6d %+8.4f %7.4f %7.3f %8d %8d' % (name, len(used), m, s, s / sref, c1, c3))
row('host (count only)', lambda r: r['hgtt'], lambda r: r['hgtt'], lambda r: r['hgtt'])
row('dev / host', lambda r: r['dgtt'], lambda r: r['hgtt'], lambda r: r['dgtt'])
row('P0 (dev, host-prop) / host', lambda r: r['P'][0], lambda r: r['hgtt'], lambda r: r['P'][0])
row('host+1ulp(a) / host', lambda r: r['P'][1], lambda r: r['hgtt'], lambda r: r['P'][1])
row('dev+1ulp(a) / host', lambda r: r['P'][3], lambda r: r['hgtt'], lambda r: r['P'][3])
for k, m in enumerate(masks):
    row('V[%d] / host' % m, lambda r, k=k: r['V'][k], lambda r: r['hgtt'], lambda r, k=k: r['V'][k])
for k, m in enumerate(masks):
    row('V[%d] / P0' % m, lambda r, k=k: r['V'][k], lambda r: r['P'][0], lambda r, k=k: r['V'][k])
same = [sum(r['V'][k] == r['P'][0] for r in used) for k in range(nv)]
print('  points with V_k == P0 exactly:', ', '.join('%d: %d' % (m, c) for m, c in zip(masks, same)))
