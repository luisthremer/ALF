#!/usr/bin/env python3
"""Ticket 207 (development only): bit identity of two runs.
usage: t207_d2_ident.py RUN_A RUN_B [--eq]
All datasets of data.h5 (with --eq: only equal-time/scalar ones, i.e. names without '_tau/') compared byte for byte (union of
both files: a dataset missing on one side counts as differing); Precision lines of info compared as text (with --eq: only
'Precision Green' and 'Precision Phase'). Exit 0 identical, 1 differing, 2 input error."""
import sys, h5py, numpy as np
a, b = sys.argv[1], sys.argv[2]
eq = "--eq" in sys.argv
def ds(p):
    o = {}
    try:
        with h5py.File(p + "/data.h5", "r") as f:
            f.visititems(lambda n, x: o.__setitem__(n, np.asarray(x[()])) if isinstance(x, h5py.Dataset) else None)
    except OSError as e:
        print("INPUT ERROR", e); sys.exit(2)
    return o
da, db = ds(a), ds(b)
keys = sorted(set(da) | set(db))
if eq:
    keys = [k for k in keys if "_tau/" not in k]
if not keys:
    print("INPUT ERROR: no datasets"); sys.exit(2)
diff = [k for k in keys if k not in da or k not in db or da[k].tobytes() != db[k].tobytes()]
sel = (lambda l: "Precision" in l and ("Green" in l or "Phase" in l)) if eq else (lambda l: "Precision" in l)
def prec(p):  # a Precision line plus its continuation line (gfortran/nvfortran wrap the Max value onto the next line)
    L = open(p + "/info").read().splitlines(); out = []
    for i, l in enumerate(L):
        if sel(l):
            out.append(l.strip() + (" " + L[i + 1].strip() if i + 1 < len(L) and ":" not in L[i + 1] else ""))
    return out
pa, pb = prec(a), prec(b)
print(f"datasets compared {len(keys)} differing {len(diff)} {diff[:4]}; Precision lines {len(pa)} {'identical' if pa == pb and pa else 'DIFFER'}")
for l in pa: print("   ", l.strip()[:120])
sys.exit(0 if (not diff and pa == pb and pa) else 1)
