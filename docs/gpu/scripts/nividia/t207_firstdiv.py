#!/usr/bin/env python3
"""Ticket 207: first divergent bin of two runs with identical seeds (NSweep = 1 per bin, so bin = sweep).
Ticket 207 addition: the tau-resolved datasets (*_tau/obser) are also reported separately, so that a regression of the device tau pass is visible.

usage: t207_firstdiv.py RUN_A/data.h5 RUN_B/data.h5
For every dataset named '*/obser' (first axis = bins) the per-bin deviation
   dev[b] = max |a-b| / max(max|a over all bins|, tiny)       (relative to the largest value of that observable)
is computed; the maximum over all datasets is printed per bin (every 10th bin and the first ones), together with
  first bin with dev > 1e-9 ("rounding level exceeded") and first bin with dev > 1e-4 ("trajectories diverged").
Rounding-level noise stays below ~1e-10; after the first accepted/rejected decision that differs, deviations are O(1e-2).
"""
import sys

import h5py
import numpy as np

fa, fb = h5py.File(sys.argv[1], "r"), h5py.File(sys.argv[2], "r")
names = []
fa.visititems(lambda n, o: names.append(n) if isinstance(o, h5py.Dataset) and n.endswith("obser") else None)
nb = None
dev = None
dev_tau = None
dev_eq = None
ntau = 0
for n in sorted(names):
    a, b = fa[n][()], fb[n][()]
    if a.shape != b.shape:
        print("shape mismatch", n)
        sys.exit(2)
    nbin = a.shape[0]
    d = np.abs(a - b).reshape(nbin, -1).max(axis=1) / max(np.abs(a).max(), 1e-300)
    dev = d if dev is None else np.maximum(dev, d)
    if "_tau/" in n:
        dev_tau = d if dev_tau is None else np.maximum(dev_tau, d)
        ntau += 1
    else:
        dev_eq = d if dev_eq is None else np.maximum(dev_eq, d)
    nb = nbin
print("bins:", nb, "datasets:", len(names))
# a NaN/Inf in either run would make every "dev > tol" test false and look like a pass: fail explicitly instead
bad = [n for n in sorted(names) if not (np.isfinite(fa[n][()]).all() and np.isfinite(fb[n][()]).all())]
if bad or len(names) == 0 or not np.isfinite(dev).all():
    print("FAIL: non-finite values or no datasets:", bad[:5])
    sys.exit(2)
for k in list(range(min(5, nb))) + list(range(9, nb, max(1, nb // 12))):
    print(f"  bin {k+1:5d}  max dev {dev[k]:.3e}")
f9 = int(np.argmax(dev > 1e-9)) + 1 if (dev > 1e-9).any() else None
f4 = int(np.argmax(dev > 1e-4)) + 1 if (dev > 1e-4).any() else None
print("max dev before the first 1e-4 bin:", f"{dev[:(f4 - 1) if f4 else nb].max():.3e}")
print("first bin with dev > 1e-9:", f9, " first bin with dev > 1e-4:", f4, " overall max dev:", f"{dev.max():.3e}")

if dev_tau is not None:
    print(f"tau-resolved datasets ({ntau}): overall max dev {dev_tau.max():.3e}, first bin > 1e-9: "
          f"{int(np.argmax(dev_tau > 1e-9)) + 1 if (dev_tau > 1e-9).any() else None}, first bin > 1e-4: "
          f"{int(np.argmax(dev_tau > 1e-4)) + 1 if (dev_tau > 1e-4).any() else None}")
    print(f"equal-time/scalar datasets ({len(names) - ntau}): overall max dev {dev_eq.max():.3e}")
