#!/usr/bin/env python3
# Ticket 207 dev4 (rounding-order nulls of the paired test, development only): paired statistics of one ALF_GPU_TAUDIAG=2 run.log
# with variant (V), amplitude (A, amplitude 1 = extra 1-ulp chains) and null (N) chains.
#   usage: t207_d4_stat.py run.log FROM_PASS OUT.json
# Every row is "a / b" over the GTT precision numbers of the points of passes >= FROM_PASS: n, mean and sd (population) of
# log10(a/b), the sd ratio to the gate's reference sd(P2/P1) = sd(host+1ulp(b) / host+1ulp(a)) and to the other 1-ulp references,
# and the three conditions of acceptance 3 applied descriptively (|mean| <= 0.05, sd ratio <= 1.25, counts above 1e-7 and 3e-7
# within 2 sqrt(n_a + n_b) over the same points).  The gate row "dev / host" also gets the counts over all passes (as
# t207_diagstat.py).  Correctness of the nulls: per chain the median / max / count-zero of Nd (max rel. difference of the chain's
# new GTT from the host chain's new GTT, same point).  Fail-closed: exit 2 on a malformed line, a missing block or a non-positive value.
import json, math, re, statistics as st, sys
fn, frm, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
def die(msg):
    print('FAIL:', msg); sys.exit(2)
masks, amps, codes = [], [], []
seed = None
rows = []
for l in open(fn):
    if 'GPU tau pass: DEBUG ALF_GPU_TAUDIAG=' in l:
        m = re.search(r'ALF_GPU_TAUDIAG_VAR masks:([^)]*)\)', l); masks = [int(x) for x in m.group(1).split()] if m else []
        m = re.search(r'ALF_GPU_TAUDIAG_AMP:([^)]*)\)', l); amps = [float(x) for x in m.group(1).split()] if m else []
        m = re.search(r'ALF_GPU_TAUDIAG_NUL:([^)]*)\)', l); codes = [int(x) for x in m.group(1).split()] if m else []
        m = re.search(r'perturbation seed (\d+)', l); seed = int(m.group(1)) if m else None
        continue
    if not l.lstrip().startswith('TAUDIAG pt '):
        continue
    f = l.split()
    try:
        ps = int(f[f.index('pass') + 1])
        i = f.index('dev'); j = f.index('hostfull'); k = f.index('P')
        r = dict(ps=ps, dev=[float(x) for x in f[i + 1:i + 5]], host=[float(x) for x in f[j + 1:j + 5]],
                 P=[float(x) for x in f[k + 1:k + 6]])
        def blk(name, n):
            if name not in f:
                return None
            a = f.index(name)
            v = [float(x) for x in f[a + 1:a + 1 + n]]
            if len(v) != n:
                die('block %s with %d values, expected %d' % (name, len(v), n))
            return v
        r['V'] = blk('V', len(masks)); r['A'] = blk('A', len(amps)); r['Ad'] = blk('Ad', len(amps) + 1)
        r['N'] = blk('N', len(codes)); r['Nd'] = blk('Nd', len(codes))
    except (ValueError, IndexError):
        die('malformed TAUDIAG line: ' + l[:100])
    if not all(math.isfinite(x) for x in r['dev'] + r['host'] + r['P']):
        die('non-finite value in pass %d' % ps)
    rows.append(r)
if not rows:
    die('no TAUDIAG lines')
used = [r for r in rows if r['ps'] >= frm]
if len(used) < 2:
    die('only %d rows from pass %d' % (len(used), frm))
for r in used:
    for b, n in (('V', len(masks)), ('A', len(amps)), ('N', len(codes))):
        if n and r[b] is None:
            die('block %s missing in pass %d' % (b, r['ps']))
# GTT precision numbers of every chain, by name
def series(r):
    s = {'dev': r['dev'][1], 'host': r['host'][1]}
    for q in range(5):
        s['P%d' % q] = r['P'][q]
    for q, m in enumerate(masks):
        s['V%d[%d]' % (q + 1, m)] = r['V'][q]
    for q, a in enumerate(amps):
        s['A%d[%g]' % (q + 1, a)] = r['A'][q]
    for q, c in enumerate(codes):
        s['N%d[%d]' % (q + 1, c)] = r['N'][q]
    return s
S = [series(r) for r in used]
def lr(a, b):
    x = []
    for s in S:
        if not (s[a] > 0 and s[b] > 0):
            die('non-positive value in %s / %s' % (a, b))
        x.append(math.log10(s[a] / s[b]))
    return x
def sdof(a, b):
    return st.pstdev(lr(a, b))
refs = {'P2/P1': sdof('P2', 'P1')}
ones = [k for k in S[0] if k.startswith('A') and k.endswith('[1]')]
for q in range(0, len(ones) - 1, 2):
    refs['%s/%s' % (ones[q + 1], ones[q])] = sdof(ones[q + 1], ones[q])
sref = refs['P2/P1']
def row(a, b, label):
    x = lr(a, b)
    m, sd = st.mean(x), st.pstdev(x)
    cnt = {}
    for t in (1e-7, 3e-7):
        na, nb = sum(s[a] > t for s in S), sum(s[b] > t for s in S)
        cnt['%g' % t] = [na, nb, abs(na - nb) <= 2 * math.sqrt(na + nb)]
    return dict(label=label, a=a, b=b, n=len(x), mean=m, sd=sd, ratio=sd / sref,
                ratio_alt={k: sd / v for k, v in refs.items()}, max_abs=max(abs(v) for v in x),
                cond_mean=abs(m) <= 0.05, cond_sd=sd / sref <= 1.25, counts=cnt,
                cond_counts=all(c[2] for c in cnt.values()))
R = []
R.append(row('dev', 'host', 'device (gate row)'))
R.append(row('P1', 'host', '1-ulp(a) / host'))
R.append(row('P2', 'P1', '1-ulp(b) / 1-ulp(a) (gate reference)'))
for k in refs:
    if k != 'P2/P1':
        a, b = k.split('/'); R.append(row(a, b, 'extra 1-ulp pair'))
R.append(row('P0', 'host', 'device output host-propagated / host'))
R.append(row('dev', 'P0', 'device propagation share (dev / P0)'))
R.append(row('P4', 'P3', '1-ulp pair on the device udvr'))
nk = {c: 'N%d[%d]' % (q + 1, c) for q, c in enumerate(codes)}
vk = {m: 'V%d[%d]' % (q + 1, m) for q, m in enumerate(masks)}
for c, k in nk.items():
    R.append(row(k, 'host', 'null %d / host' % c))
# symmetric pairs: same flags and block size, different stream; and every pair of compact-WY block sizes
byfl = {}
for c in codes:
    byfl.setdefault(c % 100000, []).append(c)
for base, cs in byfl.items():
    for q in range(len(cs) - 1):
        R.append(row(nk[cs[q + 1]], nk[cs[q]], 'null pair %d vs %d' % (cs[q + 1], cs[q])))
wy = [c for c in codes if (c % 100) & 8 and c // 100000 == 0]
for q in range(len(wy) - 1):
    R.append(row(nk[wy[q + 1]], nk[wy[q]], 'null pair %d vs %d' % (wy[q + 1], wy[q])))
if 16 in nk:
    R.append(row(nk[16], 'host', '206 input share (host CGR2_2 on the device udvr / host)'))
    R.append(row('P0', nk[16], 'device CGR2_2 share (P0 / host CGR2_2 on the device udvr)'))
for m, k in vk.items():
    R.append(row(k, 'host', 'variant %d / host' % m))
    R.append(row(k, 'P0', 'variant %d / P0' % m))
# gate counts over all passes (t207_diagstat.py --check)
gc = {}
for t in (1e-7, 3e-7):
    nd, nh = sum(r['dev'][1] > t for r in rows), sum(r['host'][1] > t for r in rows)
    gc['%g' % t] = [nd, nh, abs(nd - nh) <= 2 * math.sqrt(nd + nh)]
corr = {}
for q, c in enumerate(codes):
    d = [r['Nd'][q] for r in used]
    corr[str(c)] = dict(median=st.median(d), max=max(d), zero=sum(v == 0 for v in d), n=len(d),
                        identical_precision_to_host=sum(r['N'][q] == r['host'][1] for r in used))
if amps:
    d = [r['Ad'][len(amps)] for r in used]
    corr['device'] = dict(median=st.median(d), max=max(d), zero=sum(v == 0 for v in d), n=len(d))
    for q, a in enumerate(amps):
        d = [r['Ad'][q] for r in used]
        corr['amp%d[%g]' % (q + 1, a)] = dict(median=st.median(d), max=max(d), zero=sum(v == 0 for v in d), n=len(d))
eq = {k: sum(s[k] == s['P0'] for s in S) for k in vk.values()}
res = dict(file=fn, from_pass=frm, points=len(rows), used=len(used), seed=seed, masks=masks, amps=amps, codes=codes, refs=refs,
           rows=R, gate_counts_all_passes=gc, correctness=corr, V_equal_P0=eq, series=S)
json.dump(res, open(out, 'w'), indent=1)
print('file %s  points %d  used (pass >= %d) %d  perturbation seed %s' % (fn, len(rows), frm, len(used), seed))
print('references: ' + ', '.join('%s %.4f' % kv for kv in refs.items()))
print('  %-58s %5s %8s %7s %7s %s' % ('row', 'n', 'mean', 'sd', 'ratio', 'cnt>1e-7  cnt>3e-7  [mean sd cnt]'))
for x in R:
    c = x['counts']
    print('  %-58s %5d %+8.4f %7.4f %7.3f %4d/%-4d %4d/%-4d [%s %s %s]' % (
        (x['label'] + ': ' + x['a'] + '/' + x['b'])[:58], x['n'], x['mean'], x['sd'], x['ratio'], c['1e-07'][0], c['1e-07'][1],
        c['3e-07'][0], c['3e-07'][1], 'ok' if x['cond_mean'] else 'FAIL', 'ok' if x['cond_sd'] else 'FAIL',
        'ok' if x['cond_counts'] else 'FAIL'))
print('gate counts (all passes, dev vs host): ' + ', '.join('>%s %d/%d %s' % (t, v[0], v[1], 'ok' if v[2] else 'FAIL') for t, v in gc.items()))
print('correctness (rel. difference of the new GTT from the host CGR2_2 GTT, same point):')
for k, v in corr.items():
    print('  %-14s median %.2e  max %.2e  bitwise equal %d/%d%s' % (k, v['median'], v['max'], v['zero'], v['n'],
          ('  precision number identical to host %d' % v['identical_precision_to_host']) if 'identical_precision_to_host' in v else ''))
if eq:
    print('V chains with precision number == P0: ' + ', '.join('%s %d' % kv for kv in eq.items()))
