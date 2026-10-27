#!/usr/bin/env python3
"""Ticket 108: for finished runs in RUNDIR compute res.json (E, S(pi,pi)) and delete the bulky data.h5/res/res.pkl
(disk limit). L6 cases are kept (small). usage: post108.py RUNDIR"""
import glob, os, shutil, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze108 import ok
from analyze import obs
for d in sorted(glob.glob(f"{sys.argv[1]}/*/*")):
    if not os.path.exists(f"{d}/timing.txt"):
        continue
    if ok(d) and not os.path.exists(f"{d}/res.json") and os.path.exists(f"{d}/data.h5"):
        obs(d)
    if os.path.exists(f"{d}/res.json") or not ok(d):
        if "/L6_" in d:
            continue
        for f in ("data.h5", "res.pkl"):
            if os.path.exists(f"{d}/{f}"): os.remove(f"{d}/{f}")
        if os.path.isdir(f"{d}/res"): shutil.rmtree(f"{d}/res")
