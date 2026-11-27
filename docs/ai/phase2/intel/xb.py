#!/usr/bin/env python3
"""Two binaries (or lanes), same seeds, rule of p25_run.py. usage: xb.py EXE_A ENV_A EXE_B ENV_B START WORK CASES.py label[,..]
ENV_x: comma list KEY=VAL or -"""
import os, sys
ea, xa, eb, xb, start, work, cases, labels = sys.argv[1:9]
pe = lambda s: {} if s == "-" else dict(a.split("=", 1) for a in s.split(","))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(cases)), "..", "p25_runs"))
import importlib.util, p25_run as pr
spec = importlib.util.spec_from_file_location("cases", cases); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
for label, edits, expect, env, note in m.CASES:
    if label not in labels.split(","): continue
    da, db = os.path.join(work, label, "a"), os.path.join(work, label, "b")
    pr.prepare(da, start, edits); pr.prepare(db, start, edits)
    ra, _ = pr.run(ea, da, dict(env, **pe(xa)), 0, 7200)
    rb, _ = pr.run(eb, db, dict(env, **pe(xb)), 0, 7200)
    gp = pr.e2e.green_precision(db); tol = max(1e-11, 100 * gp[0]) if gp else 1e-11
    bad, diffs = pr.compare(da, db, False, tol)
    worst = sorted(((v[1], k) for k, v in diffs.items() if v), reverse=True)[:4]
    print("%s rc=%d/%d max scaled diff %.2e (tol %.1e) top: %s" % (label, ra, rb, worst[0][0] if worst else 0.0, tol,
          ", ".join("%s %.2e" % (k, v) for v, k in worst)), flush=True)
