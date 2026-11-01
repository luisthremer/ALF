#!/usr/bin/env python3
# Ticket 207 dev4 (development only): secondary tables from the stat.json files of t207_d4_stat.py (generated, no hand-made numbers).
#   usage: t207_d4_share.py OUT.json LABEL=stat.json [...]      (one path per call)
# (1) Reference noise: the device sd ratio against each in-run 1-ulp reference (P2/P1 = the gate's, A2/A1, A4/A3).
# (2) Share of the 206 input: log(dev/host) = log(dev/P0) + log(P0/N16) + log(N16/host) per point (propagation, device CGR2_2 on the
#     device udvr, host CGR2_2 device-vs-host udvr); per set the variance of each term and of the total, and the covariance part.
# (3) Statistics side by side for the device, the in-run control V[8] (stand-in for TAU_NS=0 / TAUSWAP=8) and the nulls:
#     sd ratio to P2/P1 (gate), |mean| (gate), and the sd ratio to a fully decorrelated symmetric host null (row-permutation pair
#     N[200004] / N[100004]).
import json, math, statistics as st, sys
out, args = sys.argv[1], sys.argv[2:]
sets = [(a.split('=', 1)[0], json.load(open(a.split('=', 1)[1]))) for a in args]
def key(s, prefix, code):
    for k in s['series'][0]:
        if k.startswith(prefix) and k.endswith('[%s]' % code):
            return k
    return None
def lr(s, a, b):
    return [math.log10(p[a] / p[b]) for p in s['series']]
res = dict(refnoise={}, share={}, stats={})
print('(1) device sd ratio against each in-run 1-ulp reference')
for l, s in sets:
    sd = st.pstdev(lr(s, 'dev', 'host'))
    res['refnoise'][l] = {k: sd / v for k, v in s['refs'].items()}
    print('  %s  dev/host sd %.4f  ' % (l, sd) + '  '.join('%s %.4f -> %.3f' % (k, v, sd / v) for k, v in s['refs'].items()))
print('(2) variance decomposition of log10(dev/host): var(total) = var(prop) + var(cgr) + var(input) + 2 cov')
for l, s in sets:
    n16 = key(s, 'N', 16)
    if n16 is None:
        continue
    a, b, c = lr(s, 'dev', 'P0'), lr(s, 'P0', n16), lr(s, n16, 'host')
    t = [x + y + z for x, y, z in zip(a, b, c)]
    va, vb, vc, vt = st.pvariance(a), st.pvariance(b), st.pvariance(c), st.pvariance(t)
    sref = s['refs']['P2/P1']
    res['share'][l] = dict(var_total=vt, var_prop=va, var_cgr=vb, var_input=vc, cov_part=vt - va - vb - vc,
                           ratio_total=math.sqrt(vt) / sref, ratio_prop=math.sqrt(va) / sref, ratio_cgr=math.sqrt(vb) / sref,
                           ratio_input=math.sqrt(vc) / sref)
    print('  %s  sd ratios: total %.3f = prop %.3f (+) device CGR2_2 %.3f (+) 206 input %.3f;  variance fractions: prop %.2f  cgr %.2f  '
          'input %.2f  covariances %+.2f' % (l, math.sqrt(vt) / sref, math.sqrt(va) / sref, math.sqrt(vb) / sref, math.sqrt(vc) / sref,
                                              va / vt, vb / vt, vc / vt, (vt - va - vb - vc) / vt))
print('(3) statistics side by side (per set: sd ratio to P2/P1 | mean | sd ratio to the row-permutation null pair)')
rows = [('device dev/host', 'dev', 'host')]
for l, s in sets[:1]:
    v8 = key(s, 'V', 8)
    if v8:
        rows.append(('control V[8]/host', v8, 'host'))
    for c in s['codes']:
        rows.append(('null %d/host' % c, key(s, 'N', c), 'host'))
    if key(s, 'N', 100004) and key(s, 'N', 200004):
        rows.append(('null pair 200004/100004', key(s, 'N', 200004), key(s, 'N', 100004)))
for name, a, b in rows:
    cells, allx, allr, allp = [], [], [], []
    for l, s in sets:
        x = lr(s, a, b); sd = st.pstdev(x); m = st.mean(x)
        rp = st.pstdev(lr(s, key(s, 'N', 200004), key(s, 'N', 100004))) if key(s, 'N', 100004) else float('nan')
        cells.append((sd / s['refs']['P2/P1'], m, sd / rp))
        allx += x; allr += lr(s, 'P2', 'P1')
        if key(s, 'N', 100004):
            allp += lr(s, key(s, 'N', 200004), key(s, 'N', 100004))
    pooled = (st.pstdev(allx) / st.pstdev(allr), st.mean(allx), st.pstdev(allx) / st.pstdev(allp) if allp else float('nan'))
    res['stats'][name] = dict(per_set={l: c for (l, _), c in zip(sets, cells)}, pooled=pooled)
    print('  %-26s ' % name + '  '.join('%5.3f %+.3f %5.3f' % c for c in cells) + '   pooled %5.3f %+.3f %5.3f' % pooled)
json.dump(res, open(out, 'w'), indent=1)
