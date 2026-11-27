#!/usr/bin/env python3
"""CPU rounding floor of P2.6 cases: the same CPU binary twice with the same seeds, MKL_CBWR=AUTO vs COMPATIBLE (only the MKL
code path changes), compared by the rule of p25_run.py (scaled diff per file, tol = max(1e-11, 100 x Precision Green mean)).
usage: floor.py EXE START WORK CASES.py label[,label...]"""
import os, sys
exe, start, work, cases, labels = sys.argv[1:6]
P25 = os.path.join(os.path.dirname(os.path.abspath(cases)), "..", "p25_runs")
sys.path.insert(0, P25)
import importlib.util
import p25_run as pr
spec = importlib.util.spec_from_file_location("cases", cases); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
want = labels.split(",")
for label, edits, expect, env, note in m.CASES:
    if label not in want:
        continue
    da, db = os.path.join(work, label, "auto"), os.path.join(work, label, "compat")
    pr.prepare(da, start, edits); pr.prepare(db, start, edits)
    ra, _ = pr.run(exe, da, dict(env, MKL_CBWR="AUTO"), 0, 3600)
    rb, _ = pr.run(exe, db, dict(env, MKL_CBWR="COMPATIBLE"), 0, 3600)
    gp = pr.e2e.green_precision(da)
    tol = max(1e-11, 100 * gp[0]) if gp else 1e-11
    bad, diffs = pr.compare(db, da, False, tol)
    worst = sorted(((v[1], k) for k, v in diffs.items() if v), reverse=True)[:3]
    print("%s rc=%d/%d floor max scaled diff %.2e (tol %.1e) top: %s" % (label, ra, rb, worst[0][0] if worst else 0.0, tol,
          ", ".join("%s %.2e" % (k, v) for v, k in worst)), flush=True)
