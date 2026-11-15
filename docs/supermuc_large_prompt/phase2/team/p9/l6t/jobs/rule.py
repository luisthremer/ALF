#!/usr/bin/env python3
"""Decision 43 rule. usage: rule.py RESFILE VAR [VAR...]
RES lines: 'RES <sys> k=<k> <var> worst <dev> acc<acc>' (worst rel dev of compare.py, same-seed GPU vs ifx CPU).
L6T: VAR passes if median(dev_VAR) <= 2 median(dev_oneMKL) and max(dev_VAR) <= 3 max(dev_oneMKL) over the seed sets.
Mz=F: every L6F / L12F entry of VAR <= 1e-13. Median of an even count = mean of the two middle values."""
import sys, statistics as st
d = {}
for l in open(sys.argv[1]):
    f = l.split(); s, k, v, w = f[1], int(f[2][2:]), f[3], f[5]
    d[(s, k, v)] = float(w) if w not in ("", "nan") else float("inf")
ks = sorted({k for (s, k, v) in d if s == "L6T"})
print("seed sets:", ks)
def col(s, v): return [d.get((s, k, v), float("inf")) for k in ks]
base = col("L6T", "oneMKL")
mb, xb = st.median(base), max(base)
print("%-8s" % "k" + "".join("%11d" % k for k in ks))
for v in ["oneMKL"] + sys.argv[2:]:
    for s in ("L6T", "L6F", "L12F"):
        if any((s, k, v) in d for k in ks): print("%-8s" % (s[:4] + ":" + v)[:8] + "".join("%11.3e" % x for x in col(s, v)) + "   %s %s" % (s, v))
print("oneMKL L6T: median %.3e  max %.3e" % (mb, xb))
for v in sys.argv[2:]:
    c = col("L6T", v); m, x = st.median(c), max(c)
    okm, okx = m <= 2 * mb, x <= 3 * xb
    mzf = [d[(s, k, v)] for s in ("L6F", "L12F") for k in ks if (s, k, v) in d]
    nexp = sum(1 for s in ("L6F", "L12F") if any((s, k, v) in d for k in ks)) * len(ks)
    okf = len(mzf) == nexp and all(y <= 1e-13 for y in mzf)
    print("RULE %-6s L6T median %.3e (ratio %.2f, limit 2) %s | max %.3e (ratio %.2f, limit 3) %s | Mz=F max %s over %d/%d runs %s | L6T part %s | overall %s" % (
        v, m, m / mb, "ok" if okm else "FAIL", x, x / xb, "ok" if okx else "FAIL", "%.3e" % max(mzf) if mzf else "-", len(mzf), nexp,
        "ok" if okf else "FAIL", "PASS" if okm and okx else "FAIL", "PASS" if okm and okx and okf else "FAIL"))
