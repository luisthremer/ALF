#!/usr/bin/env python3
"""Compare two ALF run directories numerically: every observable file (*_scal, *_eq, *_local, *_localtau),
confout_0 (field configuration) and the acceptance from `info`.
usage: compare.py DIR_A DIR_B [--tol 1e-8]   (exit 1 if max/scaleative deviation > tol)"""
import sys, os, re, glob
num = re.compile(r"[-+]?\d+\.\d*(?:[EeDd][-+]?\d+)?|[-+]?\d+(?:[EeDd][-+]?\d+)|[-+]?\d+")
def floats(path):
    return [float(x.replace("D","E")) for x in num.findall(open(path).read())]
def accept(d):
    for l in open(os.path.join(d,"info")):
        if l.strip().startswith("Acceptance  "): return float(l.split(":")[1])
    return float("nan")
a, b = sys.argv[1], sys.argv[2]
tol = float(sys.argv[sys.argv.index("--tol")+1]) if "--tol" in sys.argv else 1e-8
names = sorted(os.path.basename(p) for p in glob.glob(os.path.join(a,"*")) if not p.endswith(("_info","parameters","seeds","info","run.log","TREE","RUNNING")))
worst = 0.0; bad = []
print(f"{'file':22s} {'n':>7s} {'max|dif|':>11s} {'max/scale':>11s}")
for n in names:
    pa, pb = os.path.join(a,n), os.path.join(b,n)
    if not os.path.exists(pb): bad.append(n+" (missing)"); print(f"{n:22s} MISSING in B"); continue
    x, y = floats(pa), floats(pb)
    if len(x) != len(y): bad.append(n+" (shape)"); print(f"{n:22s} SHAPE {len(x)} vs {len(y)}"); continue
    d = [abs(u-v) for u,v in zip(x,y)]
    scale = max([abs(u) for u in x], default=0) or 1.0   # file-scale norm: avoids blow-up on entries ~ 0
    rel = [dd/scale for dd in d]
    mr = max(rel, default=0); md = max(d, default=0)
    worst = max(worst, mr)
    print(f"{n:22s} {len(x):7d} {md:11.3e} {mr:11.3e}")
aa, ab = accept(a), accept(b)
print(f"acceptance A={aa:.12f} B={ab:.12f}  diff={abs(aa-ab):.2e}")
ok = worst <= tol and not bad and abs(aa-ab) < 1e-12
print(("PASS" if ok else "FAIL"), f"worst rel dev {worst:.3e} (tol {tol:g})", bad if bad else "")
sys.exit(0 if ok else 1)
