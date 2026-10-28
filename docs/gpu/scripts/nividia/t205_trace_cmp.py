#!/usr/bin/env python3
"""Ticket 205: compare two ALF_GPU_TRACE files (one line per WRAPGRUP/WRAPGRDO call: direction, slice, flipped-site
string, Phase real, Phase imag).  usage: t205_trace_cmp.py TRACE_A TRACE_B [LTROT=40]
Reports the first call whose flip pattern differs (sweep, direction, slice, first differing site) or 'none', and the
maximum Phase difference before that point."""
import sys

a = [l.split() for l in open(sys.argv[1])]
b = [l.split() for l in open(sys.argv[2])]
ltrot = int(sys.argv[3]) if len(sys.argv) > 3 else 40
per_sweep = 2 * ltrot
print(f"calls: A {len(a)}, B {len(b)} (sweep = {per_sweep} calls: {ltrot} WRAPGRUP + {ltrot} WRAPGRDO)")
first = None
dphase = 0.0
nflip = 0
for i, (x, y) in enumerate(zip(a, b)):
    if x[0] != y[0] or x[1] != y[1]:
        first = (i, "call sequence differs"); break
    if x[2] != y[2]:
        site = next(k for k, (p, q) in enumerate(zip(x[2], y[2])) if p != q) + 1
        first = (i, f"flip pattern differs, first site {site}"); break
    nflip += x[2].count("1")
    d = abs(complex(float(x[3]), float(x[4])) - complex(float(y[3]), float(y[4])))
    dphase = max(dphase, d)
if first is None and len(a) != len(b):
    first = (min(len(a), len(b)), "different number of calls")
if first is None:
    print(f"first differing decision: none; flips compared {nflip}; max |Phase diff| {dphase:.3e}")
else:
    i, why = first
    x = a[i]
    print(f"first differing decision: call {i+1}: sweep {i//per_sweep+1}, direction {x[0]}, slice {x[1]}: {why}; "
          f"flips agreed before: {nflip}; max |Phase diff| before: {dphase:.3e}")
