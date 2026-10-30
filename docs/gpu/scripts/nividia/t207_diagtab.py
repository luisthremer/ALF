#!/usr/bin/env python3
# Ticket 207 (LOG precision tail): summary of the ALF_GPU_TAUDIAG=1 lines of a run.log.
#   usage: t207_diagtab.py run.log [NTOP]
# Prints the overall maxima per column (device, full host chain, host CGR2_2 on device inputs (hC), host CGR2_2 on the host
# one-block WRAPUR (hD), dev_cgr22 with host twins E1 V1INV, E2 ZTRSM, E4 host ZUNMQR, E8 device ZUNMQR w/o NS, E7 = 1+2+4)
# (P: GTT chains propagated on the host: device output; host CGR2_2 on the host / device udvr with U perturbed by 1 ulp)
# and the NTOP points with the largest device precision number.
import re, sys
rows = []
for l in open(sys.argv[1]):
    if 'TAUDIAG' not in l:
        continue
    f = l.split()
    g = lambda k: f.index(k)
    dev = [float(x) for x in f[g('dev') + 1:g('dev') + 5]]
    hf = [float(x) for x in f[g('hostfull') + 1:g('hostfull') + 5]]
    hc = float(f[g('hC') + 1]); hd = float(f[g('hD') + 1])
    e = [float(x) for x in f[g('E1,2,4,8,7') + 1:g('E1,2,4,8,7') + 6]]
    P = [float(x) for x in f[g('P') + 1:g('P') + 6]] if 'P' in f else [0.0] * 5
    rest = l.split('| dCd')[1].split()
    rows.append(dict(pt=int(f[g('pt') + 1]), ps=int(f[g('pass') + 1]), nst=int(f[g('nst') + 1]), nfe=int(f[g('nfe') + 1]),
                     br=int(f[g('br') + 1]), P=P, dev=dev, hf=hf, hc=hc, hd=hd, e=e, tail=' '.join(rest)))
if not rows:
    sys.exit('no TAUDIAG lines')
cols = [('dev', lambda r: max(r['dev'])), ('hostfull', lambda r: max(r['hf'])), ('hC', lambda r: r['hc']), ('hD', lambda r: r['hd'])] + \
       [('E%s' % m, (lambda i: lambda r: r['e'][i])(i)) for i, m in enumerate(['1', '2', '4', '8', '7'])]
print('points %d;  overall max:' % len(rows), '  '.join('%s %.3e' % (n, max(fn(r) for r in rows)) for n, fn in cols))
print('mean of per-point max:  ', '  '.join('%s %.3e' % (n, sum(fn(r) for r in rows) / len(rows)) for n, fn in cols))
for k, nm in enumerate(['G00', 'GTT', 'GT0', 'G0T']):
    print('pair %s: dev max %.3e  hostfull max %.3e' % (nm, max(r['dev'][k] for r in rows), max(r['hf'][k] for r in rows)))
ntop = int(sys.argv[2]) if len(sys.argv) > 2 else 15
print('top %d points by device max:' % ntop)
for r in sorted(rows, key=lambda r: -max(r['dev']))[:ntop]:
    print('pt %5d pass %4d nst %d nfe %d br %d dev %s | hf %.2e hC %.2e hD %.2e | E %s | dCd %s' % (
        r['pt'], r['ps'], r['nst'], r['nfe'], r['br'], ' '.join('%.2e' % x for x in r['dev']), max(r['hf']), r['hc'], r['hd'],
        ' '.join('%.2e' % x for x in r['e']), r['tail']))
    print('      GTT: dev %.2e hostfull %.2e | chains: device output host-propagated %.2e; host udvr + 1 ulp %s; device udvr + 1 ulp %s' % (
        r['dev'][1], r['hf'][1], r['P'][0], ' '.join('%.2e' % x for x in r['P'][1:3]), ' '.join('%.2e' % x for x in r['P'][3:5])))
