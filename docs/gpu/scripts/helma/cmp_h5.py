#!/usr/bin/env python3
"""Compare all datasets of Temp_0/data.h5 (and confout) of two PP run dirs.
usage: cmp_h5.py DIR_A DIR_B   (run with /hnvme/workspace/b133ae29-CLAUDE/venv_h5/bin/python: needs h5py, numpy)
deviation = max|a-b| / max|a| per dataset (file-scale norm as in compare.py)."""
import sys, os, glob, re, h5py, numpy as np
a, b = sys.argv[1], sys.argv[2]
def items(f):
    out = {}
    def v(n, o):
        if isinstance(o, h5py.Dataset) and o.dtype.kind in "fc": out[n] = o[()]
    f.visititems(v); return out
worst = 0.0
for pat in ("Temp_0/data.h5", "Temp_0/confout_0.h5"):
    fa, fb = os.path.join(a, pat), os.path.join(b, pat)
    if not (os.path.exists(fa) and os.path.exists(fb)): print("missing", pat); worst = 1; continue
    with h5py.File(fa) as A, h5py.File(fb) as B:
        ia, ib = items(A), items(B)
        for n in sorted(ia):
            if n not in ib: print("only in A:", n); worst = 1; continue
            x, y = np.asarray(ia[n]), np.asarray(ib[n])
            if x.shape != y.shape: print("shape", n); worst = 1; continue
            sc = np.max(np.abs(x)) or 1.0
            d = np.max(np.abs(x - y)) / sc if x.size else 0.0
            worst = max(worst, d)
            if d > 1e-12: print("%-60s %.3e" % (n, d))
print("worst rel dev %.3e over data.h5 + confout" % worst)
def grab(d):
    t = open(os.path.join(d, "Temp_0/info")).read()
    acc = re.search(r"Effective Acceptance\s*:\s*(\S+)", t); pg = re.search(r"Precision Green.*\n\s*(.*)\n", t)
    return acc.group(1) if acc else None, pg.group(1).strip() if pg else None
print("acceptance / Precision Green line:", grab(a), grab(b))
