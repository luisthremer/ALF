#!/usr/bin/env python3
"""Ticket 207-1 (development only): compare two qref.txt files (QREF lines of the offline quad-reference program) by point.
usage: t207_1_qcmp.py A B [--files-of-A-only-in-B]   A is the reference (dev5), B the new file; points are matched by 'file N pass P nst S nfe F' (the
header up to 'br'), the timing field 'tq' is ignored.  Prints: lines, how many are textually identical after removing tq, and for every
method segment the max relative difference of the numbers (0 for identical lines)."""
import re, sys
NUM = re.compile(r'[-+]?\d\.\d+E[-+]\d+')
def parse(fn):
    out = {}
    for line in open(fn):
        if not line.startswith('QREF'): continue
        segs = line.rstrip('\n').split(' | ')
        h = segs[0].split(); key = tuple(h[1:9])
        hd = ' '.join(h[1:]); hd = re.sub(r' tq +[0-9.]+', '', hd)
        out[key] = (hd, {s.split()[0]: [float(x) for x in NUM.findall(s[len(s.split()[0]):])] for s in segs[1:]}, [' '.join(s.split()) for s in segs[1:]])
    return out
a, b = parse(sys.argv[1]), parse(sys.argv[2])
common = [k for k in b if k in a]
print(f'A {sys.argv[1]}: {len(a)} points; B {sys.argv[2]}: {len(b)} points; matched {len(common)}; B not in A {len([k for k in b if k not in a])}')
ident = 0; worst = {}; hdiff = 0
for k in common:
    ha, ma, sa = a[k]; hb, mb, sb = b[k]
    if ha != hb: hdiff += 1
    if sa == sb and ha == hb: ident += 1
    for m in set(ma) | set(mb):
        va, vb = ma.get(m), mb.get(m)
        if va is None or vb is None or len(va) != len(vb): worst[m] = float('inf'); continue
        for x, y in zip(va, vb):
            d = abs(x - y) / max(abs(x), abs(y), 1e-300) if (x != y) else 0.0
            worst[m] = max(worst.get(m, 0.0), d)
print(f'  identical lines (tq and spacing ignored): {ident} of {len(common)}; headers differing: {hdiff}')
print('  max relative difference per segment: ' + ', '.join(f'{m} {w:.2e}' for m, w in sorted(worst.items())))
