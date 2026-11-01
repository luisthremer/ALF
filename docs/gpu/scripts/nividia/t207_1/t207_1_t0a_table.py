#!/usr/bin/env python3
"""Ticket 207-1 (development only): table of the SELFTEST=4 runs of t207_1_t0a.sh RUN (from the run.log files).  usage: t207_1_t0a_table.py RUN"""
import re, sys, glob, os
R = sys.argv[1]; W = os.path.expanduser("~/alf-gpu/work/t207-1/runs")
num = lambda s: [float(x) for x in re.findall(r"[0-9]\.[0-9]{3}E[-+][0-9]{2}", s)]
rows = []
for d in sorted(glob.glob(f"{W}/t0a{R}_st4*")):
    cfg = os.path.basename(glob.glob(d + "/*")[0]); L = open(f"{d}/{cfg}/run.log").read().splitlines()
    path = "LOG" if "st4log" in d else "def"
    r = dict(cfg=cfg, path=path, passed=any("tau self-test PASSED" in l for l in L))
    for i, l in enumerate(L):
        if "pass over" in l: r["slices"] = re.search(r"over (\d+) slices \((\d+) stab", l).groups(); r["pass"] = num(L[i + 1])
        if "right after CGR2_2" in l: r["after"] = num(l.split("(tol")[0])
        if "natural calls" in l: r["br"] = re.search(r"> D2\(1\): (\d+), else: (\d+)", l).groups()
        if "alone, natural states" in l: r["alone"] = num(l)
        if "D1 scaled x1e3" in l: r["s3"] = num(l.split(":")[-1])
        if "D1 scaled x1e-3" in l: r["sm3"] = num(l.split(":")[-1])
        if "precision reduction vs COMPARE" in l: r["prec"] = num(l.split("tol")[0])
    rows.append(r)
f = lambda v: " ".join("%.2e" % x for x in v)
print(f"{'config':18s} path {'slices/pts':>9s} {'branch D1>D2/else':>8s}  per-slice max G00 GT0 GTT G0T | after-CGR2_2 max | CGR2_2 alone natural max | x1e3 | x1e-3 | prec | verdict")
for r in rows:
    mx = lambda k: "%.2e" % max(r[k]) if k in r else "n/a"
    print(f"{r['cfg']:18s} {r['path']:4s} {'/'.join(r['slices']):>9s} {'/'.join(r['br']):>8s}  {f(r['pass'])} | {mx('after')} | {mx('alone')} | {mx('s3')} | {mx('sm3')} | {mx('prec')} | {'PASSED' if r['passed'] else 'NOT PASSED'}")
print("tolerance 1e-10 (relative) for the slice, after-CGR2_2 and CGR2_2-alone columns, 1e-12 for the precision reduction; x1e3 / x1e-3: CGR2_2 alone with D1 scaled (forced branches)")
