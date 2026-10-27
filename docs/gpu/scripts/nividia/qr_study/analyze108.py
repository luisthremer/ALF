#!/usr/bin/env python3
"""Ticket 108: aggregate RUNDIR/<case>_s<k>/<build> over seeds; markdown table per case.
usage: analyze108.py RUNDIR UPBUILD BUILD [BUILD ...]   (venv; pyALF for E, S(pi,pi))
Per variant: Green mean (mean over completed seeds) / Green max (max over seeds) / tau mean / tau max,
ratios to upstream, aborts n/3, E and S(pi,pi) difference to upstream in sigma (pooled over seeds)."""
import glob, os, re, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze import info, obs


def ok(d):
    return os.path.exists(f"{d}/info") and "Precision Green" in open(f"{d}/info", errors="replace").read() \
        and os.path.exists(f"{d}/timing.txt") and "Unstable matrix" not in open(f"{d}/run.log", errors="replace").read() \
        and "NaN" not in open(f"{d}/run.log", errors="replace").read()


def agg(cdirs, b):
    ds = [f"{c}/{b}" for c in cdirs if os.path.exists(f"{c}/{b}/timing.txt")]
    good = [d for d in ds if ok(d)]
    r = dict(n=len(ds), ok=len(good))
    if not good:
        return r
    I = [info(d) for d in good]
    r["gm"] = np.mean([i["green"][0] for i in I]); r["gx"] = max(i["green"][1] for i in I)
    r["tm"] = np.mean([i["tau"][0] for i in I]); r["tx"] = max(i["tau"][1] for i in I)
    r["px"] = max(i["phase"][0] for i in I)
    O = [obs(d) for d in good]
    r["E"] = np.mean([o["E"] for o in O]); r["Ee"] = np.sqrt(sum(o["E_err"] ** 2 for o in O)) / len(O)
    r["S"] = np.mean([o["S"] for o in O]); r["Se"] = np.sqrt(sum(o["S_err"] ** 2 for o in O)) / len(O)
    r["wall"] = np.mean([float(open(f"{d}/timing.txt").read().split("=")[1]) for d in good])
    return r


if __name__ == "__main__":
    run, up, builds = sys.argv[1], sys.argv[2], sys.argv[3:]
    cases = sorted({re.sub(r"_s\d+$", "", os.path.basename(p)) for p in glob.glob(f"{run}/*")})
    print("| case | build | aborts | Green mean | (x up) | Green max | (x up) | tau mean | (x up) | tau max | (x up) | Phase max | E | dE/sig | S(pi,pi) | dS/sig | wall s |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    f = lambda x: f"{x:.2e}"
    for c in cases:
        cd = sorted(glob.glob(f"{run}/{c}_s*"))
        for b in [up] + builds:
            cdb = [x for x in cd if os.path.exists(f"{x}/{b}/timing.txt")]  # compare only on seeds where b ran
            U = agg(cdb, up)
            R = U if b == up else agg(cdb, b)
            if R["n"] == 0:
                continue
            ab = f"{R['n'] - R['ok']}/{R['n']}"
            if R["ok"] == 0:
                print(f"| {c} | {b} | {ab} | | | | | | | | | | | | | | |"); continue
            rt = lambda k: (f"{R[k] / U[k]:.2f}" if U.get("ok") else "-")
            sg = lambda k, ke: (f"{(R[k] - U[k]) / np.hypot(R[ke], U[ke]):+.1f}" if U.get("ok") and b != up else "-")
            print(f"| {c} | {b} | {ab} | {f(R['gm'])} | {rt('gm')} | {f(R['gx'])} | {rt('gx')} | {f(R['tm'])} | {rt('tm')} | {f(R['tx'])} | {rt('tx')} "
                  f"| {f(R['px'])} | {R['E']:.4f}({R['Ee']:.4f}) | {sg('E','Ee')} | {R['S']:.3f}({R['Se']:.3f}) | {sg('S','Se')} | {R['wall']:.0f} |")
