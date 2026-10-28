#!/usr/bin/env python3
"""Ticket 206: first divergent bin of two runs with identical seeds (NSweep = 1 per bin, so bin = sweep).

usage: t206_firstdiv.py RUN_A/data.h5 RUN_B/data.h5
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
for n in sorted(names):
    a, b = fa[n][()], fb[n][()]
    if a.shape != b.shape:
        print("shape mismatch", n)
        sys.exit(2)
    nbin = a.shape[0]
    d = np.abs(a - b).reshape(nbin, -1).max(axis=1) / max(np.abs(a).max(), 1e-300)
    dev = d if dev is None else np.maximum(dev, d)
    nb = nbin
print("bins:", nb, "datasets:", len(names))
for k in list(range(min(5, nb))) + list(range(9, nb, max(1, nb // 12))):
    print(f"  bin {k+1:5d}  max dev {dev[k]:.3e}")
f9 = int(np.argmax(dev > 1e-9)) + 1 if (dev > 1e-9).any() else None
f4 = int(np.argmax(dev > 1e-4)) + 1 if (dev > 1e-4).any() else None
print("max dev before the first 1e-4 bin:", f"{dev[:(f4 - 1) if f4 else nb].max():.3e}")
print("first bin with dev > 1e-9:", f9, " first bin with dev > 1e-4:", f4, " overall max dev:", f"{dev.max():.3e}")
