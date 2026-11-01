#!/usr/bin/env python3
# Ticket 207 dev4 (development only): is the in-run control chain V[8] (dev_cgr22 without the Newton-Schulz step on the production
# inputs, host-propagated) a faithful stand-in for a production control run (ALF_GPU_TAUSWAP=8 resp. ALF_GPU_TAU_NS=0: the same
# computation, device-propagated)?  Compares, on the same seeds, a dev4 run with an older paired run and an older control run.
#   usage: t207_d4_ctlcmp.py NEW_run.log OLD_pd_run.log OLD_ctl_run.log FROM_PASS
# Prints: points whose host chain (hostfull, 4 numbers) and device numbers agree exactly between the runs (same trajectory and same
# production device path), and n / mean / sd of log10 ratios over the GTT pair, passes >= FROM_PASS: ctl dev / host (old control run),
# V[8] / host and V[8] / P0 (new run), with the per-point correlation of ctl-dev/host and V[8]/host.  Ratios to the new run's P2/P1.
import math, re, statistics as st, sys
new, oldpd, oldctl, frm = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
def load(fn):
    rows, masks = {}, []
    for l in open(fn):
        if 'GPU tau pass: DEBUG ALF_GPU_TAUDIAG=' in l:
            m = re.search(r'ALF_GPU_TAUDIAG_VAR masks:([^)]*)\)', l); masks = [int(x) for x in m.group(1).split()] if m else []
        if not l.lstrip().startswith('TAUDIAG pt '):
            continue
        f = l.split()
        pt = int(f[2]); ps = int(f[f.index('pass') + 1])
        i = f.index('dev'); j = f.index('hostfull'); k = f.index('P')
        r = dict(ps=ps, dev=f[i + 1:i + 5], host=f[j + 1:j + 5], P=[float(x) for x in f[k + 1:k + 6]])
        if 'V' in f:
            a = f.index('V'); r['V'] = [float(x) for x in f[a + 1:a + 1 + len(masks)]]
        rows[pt] = r
    return rows, masks
N, mN = load(new); O, _ = load(oldpd); C, _ = load(oldctl)
if 8 not in mN:
    sys.exit('FAIL: the new run has no V[8] chain')
q8 = mN.index(8)
common = sorted(set(N) & set(O) & set(C))
print('points new %d old pd %d old ctl %d common %d' % (len(N), len(O), len(C), len(common)))
print('hostfull identical new = old pd: %d, new = old ctl: %d; device numbers identical new = old pd: %d' % (
    sum(N[p]['host'] == O[p]['host'] for p in common), sum(N[p]['host'] == C[p]['host'] for p in common),
    sum(N[p]['dev'] == O[p]['dev'] for p in common)))
used = [p for p in common if N[p]['ps'] >= frm]
def lg(a, b):
    return math.log10(a / b)
x_ctl = [lg(float(C[p]['dev'][1]), float(C[p]['host'][1])) for p in used]
x_v8 = [lg(N[p]['V'][q8], float(N[p]['host'][1])) for p in used]
x_v8p0 = [lg(N[p]['V'][q8], N[p]['P'][0]) for p in used]
x_dev = [lg(float(N[p]['dev'][1]), float(N[p]['host'][1])) for p in used]
sref = st.pstdev([lg(N[p]['P'][2], N[p]['P'][1]) for p in used])
for name, x in (('old ctl: dev / host', x_ctl), ('new: V[8] / host', x_v8), ('new: V[8] / P0', x_v8p0), ('new: dev / host', x_dev)):
    print('  %-22s n %4d  mean %+.4f  sd %.4f  ratio to new P2/P1 (%.4f) %.3f' % (name, len(x), st.mean(x), st.pstdev(x), sref, st.pstdev(x) / sref))
mc, mv = st.mean(x_ctl), st.mean(x_v8)
cov = sum((a - mc) * (b - mv) for a, b in zip(x_ctl, x_v8)) / len(x_ctl)
print('  per-point correlation of old ctl dev/host and new V[8]/host: %.3f' % (cov / (st.pstdev(x_ctl) * st.pstdev(x_v8))))
