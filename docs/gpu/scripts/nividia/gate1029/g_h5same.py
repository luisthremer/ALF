#!/usr/bin/env python3
"""Gate 2026-10-29: compare two ALF runs bit for bit. Exit 0 only if every dataset of data.h5 (union of both files)
and every full-precision 'Precision' line of info are identical; 1 if they differ; 2 on missing input or NaN.
usage: g_h5same.py RUN_A RUN_B"""
import sys, os, re, hashlib, math
import h5py, numpy as np
a, b = sys.argv[1], sys.argv[2]
for d in (a, b):
    for f in ("data.h5", "info"):
        if not os.path.isfile(os.path.join(d, f)):
            print("FAIL: missing", os.path.join(d, f)); sys.exit(2)
def dsets(p):
    out = {}
    with h5py.File(p, "r") as h:
        h.visititems(lambda n, o: out.__setitem__(n, o[()]) if isinstance(o, h5py.Dataset) else None)
    return out
A, B = dsets(os.path.join(a, "data.h5")), dsets(os.path.join(b, "data.h5"))
if not A or not B:
    print("FAIL: empty data.h5"); sys.exit(2)
names = sorted(set(A) | set(B))
diff, first = [], None
for n in names:
    if n not in A or n not in B:
        diff.append(n); first = first or (n, "missing in one file"); continue
    x, y = np.asarray(A[n]), np.asarray(B[n])
    if np.issubdtype(x.dtype, np.number) and not (np.isfinite(x).all() and np.isfinite(y).all()):
        print("FAIL: non-finite values in", n); sys.exit(2)
    if x.shape != y.shape or x.tobytes() != y.tobytes():
        diff.append(n)
        if first is None and x.shape == y.shape and np.issubdtype(x.dtype, np.number):
            idx = np.argwhere(x != y)[0]
            first = (n, "first differing element %s: %r vs %r" % (tuple(idx), x[tuple(idx)], y[tuple(idx)]))
        elif first is None:
            first = (n, "shape/type differs")
def prec(d):
    return [l.strip() for l in open(os.path.join(d, "info"), errors="replace") if "Precision" in l]
pa, pb = prec(a), prec(b)
if not pa:
    print("FAIL: no Precision lines in info"); sys.exit(2)
print("datasets: %d (A %d, B %d), differing: %d" % (len(names), len(A), len(B), len(diff)))
print("info Precision lines identical:", pa == pb)
for x, y in zip(pa, pb):
    if x != y:
        print("  A:", x); print("  B:", y)
if first:
    print("first differing dataset:", first[0], "-", first[1])
sys.exit(0 if not diff and pa == pb else 1)
