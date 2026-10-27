#!/usr/bin/env python3
"""Statistical comparison of two reference.json files (e.g. GPU run vs CPU reference).

usage: compare.py A.json B.json [--zmax 3] [--obs PREFIX ...]
For every observable X with an error column X_err, all entries are compared with
z = (a - b) / sqrt(da^2 + db^2). Reports per observable: n, chi2/n, max|z|, #|z|>zmax.
Entries with zero combined error are compared exactly (|a-b| <= 1e-12 * max(|a|,|b|,1)).
Exit 1 if any observable has max|z| > zmax beyond what multiple testing explains
(expected count of |z|>zmax is n*p, p=0.0027 for 3 sigma; fails if count > n*p + 3*sqrt(n*p) + 1)
or chi2/n > 1 + 4*sqrt(2/n).

SUPERSEDED by compare2.py (ticket 106: curated observables, autocorrelation-aware errors, Student-t,
family-wise alpha). Keep this script only for the old full-array chi2/n numbers.

PROVISIONAL (Phase 0, see docs/gpu/BASELINE.md): the thresholds assume independent entries with
correct error bars. Correlation functions violate both (all R/k/tau entries share the slow AF mode,
and short runs underestimate errors), so a FAIL on lattice observables is not yet evidence of a bug.
"""
import argparse
import json
import math
import sys

import numpy as np


def flat(x):
    return np.asarray(x, dtype=float).ravel()


def compare(a, b, zmax, prefixes):
    rows, ok = [], True
    p_tail = math.erfc(zmax / math.sqrt(2))
    for key in sorted(a):
        if key.endswith("_err") or key + "_err" not in a or key not in b or key + "_err" not in b:
            continue
        if prefixes and not any(key.startswith(p) for p in prefixes):
            continue
        try:
            va, ea, vb, eb = (flat(v) for v in (a[key], a[key + "_err"], b[key], b[key + "_err"]))
        except (TypeError, ValueError):
            continue
        if va.shape != vb.shape or va.size == 0:
            rows.append((key, "shape mismatch", va.shape, vb.shape))
            ok = False
            continue
        sig = np.sqrt(ea**2 + eb**2)
        exact = sig == 0
        if np.any(exact):
            scale = np.maximum(np.maximum(abs(va[exact]), abs(vb[exact])), 1)
            if np.any(abs(va[exact] - vb[exact]) > 1e-12 * scale):
                ok = False
                rows.append((key, "zero-error entries differ"))
        z = (va[~exact] - vb[~exact]) / sig[~exact]
        n = z.size
        if n == 0:
            continue
        chi2n = float(np.mean(z**2))
        nbad = int(np.sum(abs(z) > zmax))
        allowed = n * p_tail + 3 * math.sqrt(n * p_tail) + 1
        bad = nbad > allowed or chi2n > 1 + 4 * math.sqrt(2 / n)
        ok &= not bad
        rows.append((key, n, round(chi2n, 3), round(float(abs(z).max()), 2), nbad, "FAIL" if bad else "ok"))
    return rows, ok


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--zmax", type=float, default=3.0)
    ap.add_argument("--obs", nargs="*", default=[])
    args = ap.parse_args()
    with open(args.a, encoding="utf-8") as f:
        A = json.load(f)
    with open(args.b, encoding="utf-8") as f:
        B = json.load(f)
    rows, ok = compare(A, B, args.zmax, args.obs)
    print(f"{'observable':28s} {'n':>6s} {'chi2/n':>8s} {'max|z|':>7s} {'#>z':>4s}")
    for r in rows:
        if len(r) == 6:
            print(f"{r[0]:28s} {r[1]:6d} {r[2]:8.3f} {r[3]:7.2f} {r[4]:4d}  {r[5]}")
        else:
            print(*r)
    for k in ("precision_green_mean_max", "precision_phase_max", "precision_tau_mean_max"):
        print(k, A.get("info", {}).get(k), B.get("info", {}).get(k))
    print("RESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
