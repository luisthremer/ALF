#!/usr/bin/env python3
"""Ticket 104: analyse run dirs RUNDIR/<case>/<build>; print markdown rows.
usage: analyze.py RUNDIR   (needs ALF_DIR for pyALF analysis, run in the venv)
"""
import glob
import json
import os
import re
import sys

import h5py
import numpy as np
from py_alf.analysis import analysis
from py_alf.ana import load_res

PAT = {"green": r"Precision Green\s+Mean, Max\s*:\s*(\S+)\s+(\S+)",
       "phase": r"Precision Phase, Max\s*:\s*(\S+)",
       "tau": r"Precision tau\s+Mean, Max\s*:\s*(\S+)\s+(\S+)"}


def info(d):
    t = open(f"{d}/info", errors="replace").read()
    out = {}
    for k, p in PAT.items():
        m = re.findall(p, t, re.M)[-1]
        out[k] = [float(x.replace("D", "E")) for x in (m if isinstance(m, tuple) else (m,))]
    return out


def obs(d):
    if not os.path.exists(f"{d}/res.json"):
        analysis(d)
        r = load_res([d]).iloc[0]
        lat = r["SpinZ_eq_lattice"] if "SpinZ_eq_lattice" in r else None
        K = np.array(r["SpinZ_eqK"])[0, 0]; Ke = np.array(r["SpinZ_eqK_err"])[0, 0]
        # (pi,pi): k-point with largest |k|_1-distance to origin in the staggered sense = maximum of S
        i = int(np.argmax(K))
        res = {"E": float(r["Ener_scal0"]), "E_err": float(r["Ener_scal0_err"]),
               "S": float(K[i]), "S_err": float(Ke[i])}
        json.dump(res, open(f"{d}/res.json", "w"))
    return json.load(open(f"{d}/res.json"))


def first_diff(da, db, tol=1e-8):
    fa, fb = h5py.File(f"{da}/data.h5"), h5py.File(f"{db}/data.h5")
    names = ["Ener_scal/obser", "Kin_scal/obser", "Pot_scal/obser", "SpinZ_eq/obser"]
    n = min(fa[names[0]].shape[0], fb[names[0]].shape[0])
    for b in range(n):
        for nm in names:
            a, c = np.asarray(fa[nm][b]), np.asarray(fb[nm][b])
            if np.max(np.abs(a - c)) > tol * max(np.max(np.abs(a)), 1e-300):  # relative to the largest entry of the bin
                return b + 1
    return "identical"


if __name__ == "__main__":
    run = sys.argv[1]
    rows = []
    for cdir in sorted(glob.glob(f"{run}/*")):
        c = os.path.basename(cdir)
        for stab in ("def", "log"):
            up, pp = f"{cdir}/{stab}-up", f"{cdir}/{stab}-pp"
            if not (os.path.exists(f"{up}/timing.txt") and os.path.exists(f"{pp}/timing.txt")):
                continue
            ok = lambda d: os.path.exists(f"{d}/info") and "Precision Green" in open(f"{d}/info", errors="replace").read()
            if not (ok(up) and ok(pp)):
                st = lambda d: "ok" if ok(d) else ("abort: " + (re.findall(r"(Unstable matrix|NaN|Terminating)", open(f"{d}/run.log", errors="replace").read()) or ["?"])[0])
                rows.append(f"| {c} | {stab} | up: {st(up)}; pp: {st(pp)} | | | | | | | |")
                continue
            iu, ip = info(up), info(pp)
            ou, op = obs(up), obs(pp)
            fd = first_diff(up, pp)
            f = lambda x: f"{x:.2e}"
            rows.append(f"| {c} | {stab} | {f(iu['green'][0])} / {f(ip['green'][0])} | {f(iu['green'][1])} / {f(ip['green'][1])} "
                        f"| {f(iu['tau'][0])} / {f(ip['tau'][0])} | {f(iu['tau'][1])} / {f(ip['tau'][1])} "
                        f"| {f(iu['phase'][0])} / {f(ip['phase'][0])} | {fd} "
                        f"| {ou['E']:.4f}({ou['E_err']:.4f}) / {op['E']:.4f}({op['E_err']:.4f}) "
                        f"| {ou['S']:.3f}({ou['S_err']:.3f}) / {op['S']:.3f}({op['S_err']:.3f}) |")
    print("| case | stab | Green mean up/pp | Green max up/pp | tau mean up/pp | tau max up/pp | Phase max up/pp | first diff bin | E up / pp | S(pi,pi) up / pp |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    print("\n".join(rows))
