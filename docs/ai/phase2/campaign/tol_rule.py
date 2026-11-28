#!/usr/bin/env python3
"""Evaluates candidate tolerance rules on tol/tol_*.tsv.  A cell = (machine, model, scheme, L, class).  Prints, per rule, the number of
cells above the tolerance and the largest diff / tol.  Rules: A = current e2e rule max(1e-11, 100 x PG mean of the CPU run);
B = A or 10 x null (CPU Nwrap 10 vs 5); C = size scaled c * 1e-12 * Ndim^p with Ndim = sites*orbitals*N_FL-independent count.  usage: tol_rule.py"""
import csv, glob, math, re
NDIM = {"hub-mz": lambda L: L * L, "hub-su2": lambda L: L * L, "kondo": lambda L: 2 * L * L, "tV": lambda L: L * L, "sp": lambda L: L * L}
cells = []
for f in sorted(glob.glob("tol/tol_*.tsv")):
    if "_files" in f or "refused" in f: continue
    for r in csv.DictReader(open(f), delimiter="\t"):
        try: g = float(r["gpu_cpu_diff"]); n = float(r["null_nwrap"])
        except ValueError: continue
        m = re.search(r"PG gpu (\S+) cpu (\S+)", r[list(r)[-1]]); pg = float(m.group(2)) if m else 0.0
        cells.append((r["machine"], r["model"], r["scheme"], int(r["L"]), r["class"], g, n, pg))
print("cells", len(cells))
rules = {"A current: max(1e-11, 100*PG)": lambda c: max(1e-11, 100 * c[7]),
         "B: max(A, 10*null)": lambda c: max(1e-11, 100 * c[7], 10 * c[6]),
         "B': max(A, 30*null)": lambda c: max(1e-11, 100 * c[7], 30 * c[6])}
for name, fn in rules.items():
    bad = [c for c in cells if c[5] > fn(c)]
    worst = max(c[5] / fn(c) for c in cells)
    print("%-34s above tol: %3d / %d, worst diff/tol %.1f" % (name, len(bad), len(cells), worst))
    for cls in ("scal", "eq", "tau"):
        b = [c for c in bad if c[4] == cls]
        if b: print("     %s: %d, e.g. %s" % (cls, len(b), "; ".join("%s %s L%d %s %.1e>%.1e" % (c[0], c[1], c[3], c[2], c[5], fn(c)) for c in sorted(b, key=lambda c: -c[5] / fn(c))[:3])))
# size scaling of the largest diff per class
print("\nlargest diff per class and size over all models and machines (scheme-independent cells)")
for cls in ("scal", "eq", "tau"):
    for L in (4, 6, 8, 12, 16, 24, 32):
        v = [c[5] for c in cells if c[4] == cls and c[3] == L]; w = [c[6] for c in cells if c[4] == cls and c[3] == L]
        if v: print(cls, "L=%d" % L, "max gpu-cpu %.1e, max null %.1e, n=%d" % (max(v), max(w), len(v)))
