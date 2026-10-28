#!/usr/bin/env python3
"""Ticket 205: per-bin Ener_scal of two runs (data.h5, 1 sweep per bin) next to the first-divergence point.
usage: t205_ener.py H5_A H5_B [TRACE_A TRACE_B]   (traces only to translate the divergence into a bin number)
Prints, for every bin, the relative difference of Ener_scal/obser (all columns) and sign; NaN-safe (NaN counts as a
failure and is printed as nan)."""
import sys

import h5py
import numpy as np

a = h5py.File(sys.argv[1], "r")
b = h5py.File(sys.argv[2], "r")
ea, eb = a["Ener_scal/obser"][()], b["Ener_scal/obser"][()]
sa, sb = a["Ener_scal/sign"][()], b["Ener_scal/sign"][()]
rel = np.abs(ea - eb) / np.maximum(np.abs(ea), np.finfo(float).tiny)
rel = rel.reshape(rel.shape[0], -1).max(axis=1)
bad = ~np.isfinite(rel)
print("bin  max rel diff Ener_scal (E, Kin, Pot, Part columns present)   sign equal")
for i, r in enumerate(rel):
    print(f"{i+1:3d}  {r:.3e}   {sa[i] == sb[i]}")
print(f"max over bins {np.nanmax(rel):.3e}; non-finite bins: {int(bad.sum())}; identical bins: {int((rel == 0).sum())}/{len(rel)}")
if len(sys.argv) > 4:
    ta = [l.split() for l in open(sys.argv[3])]
    tb = [l.split() for l in open(sys.argv[4])]
    for i, (x, y) in enumerate(zip(ta, tb)):
        if x[2] != y[2] or x[0] != y[0] or x[1] != y[1]:
            print(f"first differing decision in bin (= sweep) {i // 80 + 1}; bins before it: max rel diff "
                  f"{np.nanmax(rel[:i // 80]) if i // 80 > 0 else 0.0:.3e}")
            break
    else:
        print("no differing decision in the traces")
