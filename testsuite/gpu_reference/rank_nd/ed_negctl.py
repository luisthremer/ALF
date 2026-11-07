#!/usr/bin/env python3
"""Negative controls for the ED references: perturb ONE ingredient of the ED (the 'planted defect') and show that the
validated MC data (summary_<case>.json from ed_validate.py) reject it.  A check that cannot fail is worthless.

usage: ed_negctl.py SUMMARY.json            (run in the pyALF venv, no QMC runs needed)
Output: for each perturbed ED the z-score (MC_extrapolated - ED_perturbed)/err of every scalar observable, next to the
z-score of the unperturbed ED.  A perturbation is 'detected' if some scalar has |z| > 5.
"""
import json
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ed_cases import CASES  # noqa: E402
from ed_validate import wfit  # noqa: E402

VARIANTS = {
    "tV": [("V*1.02", lambda s: {**s, "ham_V": s["ham_V"] * 1.02}),
           ("t*1.02", lambda s: {**s, "ham_T": s["ham_T"] * 1.02}),
           ("beta*1.02", lambda s: {**s, "beta": s["beta"] * 1.02}),
           ("Checkerboard flag flipped (dup. bonds)", lambda s: {**s, "Checkerboard": not s["Checkerboard"]})],
    "Kondo": [("JK*1.02", lambda s: {**s, "Ham_JK": s["Ham_JK"] * 1.02}),
              ("Uc*1.1", lambda s: {**s, "Ham_Uc": s["Ham_Uc"] * 1.1}),
              ("t*1.02", lambda s: {**s, "ham_T": s["ham_T"] * 1.02}),
              ("Uf*1.1 (only Constraint)", lambda s: {**s, "ham_Uf": s["ham_Uf"] * 1.1}),
              ("Checkerboard flag flipped (dup. bonds)", lambda s: {**s, "Checkerboard": not s["Checkerboard"]})],
    "Spin_Peierls": [("J*1.1", lambda s: {**s, "Ham_Jx": s["Ham_Jx"] * 1.1}),
                     ("J*1.3", lambda s: {**s, "Ham_Jx": s["Ham_Jx"] * 1.3}),
                     ("lambda*1.3", lambda s: {**s, "Ham_Lambda": s["Ham_Lambda"] * 1.3}),
                     ("lambda*2", lambda s: {**s, "Ham_Lambda": s["Ham_Lambda"] * 2.0}),
                     ("Omega0*1.3", lambda s: {**s, "Ham_Omega0": s["Ham_Omega0"] * 1.3}),
                     ("U*0.5", lambda s: {**s, "Ham_U": s["Ham_U"] * 0.5}),
                     ("beta*1.2", lambda s: {**s, "beta": s["beta"] * 1.2}),
                     ("no (J/4)Q term / dimer sign test: J*0 (Heisenberg off)", lambda s: {**s, "Ham_Jx": 1e-9})],
}


def mc_extrap(res, nm):
    dts = np.array(res["dtaus"])
    ys = np.array([res["mc"][str(d)][nm] for d in dts])
    es = np.array([res["mc"][str(d)][nm + "_err"] for d in dts])
    if np.all(es < 1e-12):
        return ys[-1], 0.0
    y0, e0, _ = wfit(dts, ys, es)
    return y0, e0


if __name__ == "__main__":
    res = json.load(open(sys.argv[1]))
    case = CASES[res["case"]]
    base = res["sim_dict"]
    ham = case["ham_name"]
    names = [k for k in res["ed"] if k.endswith("_scal0")] + [k.replace("_eqR", "_maxz") for k in res["ed"] if k.endswith("_eqR")]
    lnames = [k for k in res["ed"] if k.endswith("_eqR")]
    zs = {}
    for label, mod in [("unperturbed", lambda s: s)] + VARIANTS[ham]:
        ed = case["ed"](mod(dict(base)))
        row = {}
        for nm in names:
            if nm.endswith("_maxz"):
                continue
            y0, e0 = mc_extrap(res, nm)
            row[nm] = (y0 - ed[nm]) / e0 if e0 > 0 else (0.0 if abs(y0 - ed[nm]) < 1e-8 else float("inf"))
        for nm in lnames:   # lattice observables: largest |z| over the elements with a statistical error
            edv = np.asarray(ed[nm]).ravel()
            dts = np.array(res["dtaus"])
            vals = np.array([np.asarray(res["mc"][str(d)][nm]).ravel() for d in dts])
            errs = np.array([np.asarray(res["mc"][str(d)][nm + "_err"]).ravel() for d in dts])
            zl = [(wfit(dts, vals[:, j], errs[:, j])[0] - edv[j]) / wfit(dts, vals[:, j], errs[:, j])[1]
                  for j in range(edv.size) if np.all(errs[:, j] > 1e-12)]
            row[nm.replace("_eqR", "_maxz")] = max(zl, key=abs) if zl else 0.0
        zs[label] = row
    print(f"case {res['case']}: z = (MC extrapolated to dtau=0 - ED)/err")
    print(f"{'ED variant':42s}" + "".join(f"{n.replace('_scal0',''):>12s}" for n in names) + "   detected")
    for label, row in zs.items():
        det = "" if label == "unperturbed" else ("YES" if max(abs(v) for v in row.values()) > 5 else "no")
        print(f"{label:42s}" + "".join(f"{row[n]:12.1f}" for n in names) + f"   {det}")
