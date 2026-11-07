#!/usr/bin/env python3
"""Validate an ED reference against the CPU build (ALF_cpu.out, GNU, HDF5, noMPI): run, analyse, compare.
Runs ON nividia inside ~/alf-gpu/work/t401ed (pyALF venv, TMPDIR there, OMP_NUM_THREADS=1, <= 3 processes).

usage: ed_validate.py CASE [--dtaus 0.2,0.1,0.05] [--nsweep N] [--nbin N] [--tag T] [--stages prep,run,ana]
CASE is a key of ed_cases.CASES.  Result: runs/<tag>/<CASE>_dt<dtau>/ and runs/<tag>/summary_<CASE>.json
Statistics: pyALF analysis (jackknife over bins, N_skip bins dropped) = the usual ALF analysis, plus lag-1
autocorrelation of the bin series from data.h5.  The ED is compared (a) per dtau with the dtau-dependent MC value
(informational: Trotter error O(dtau^2)) and (b) after weighted extrapolation y0 + a dtau^2 -> dtau = 0 (the test
used by testsuite/test_vs_ed).  Nothing is adjusted to force agreement.
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"          # before numpy is imported: <= 1 core per process (machine rules)
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
import argparse
import json
import secrets
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ed_cases import CASES  # noqa: E402

W = os.path.expanduser("~/alf-gpu/work/t401ed")
BIN = f"{W}/bin/ALF_cpu.out"
os.environ.setdefault("TMPDIR", f"{W}/tmp")


def dname(case, dtau):
    return f"{case['name']}_dt{dtau:g}"


def prep(case, dtau, outroot, nsweep, nbin):
    from py_alf import ALF_source, Simulation
    sd = dict(case["sim_dict"])
    sd["dtau"] = dtau
    sd["NSweep"] = nsweep
    sd["NBin"] = nbin
    src = ALF_source(alf_dir=f"{W}/src")
    sim = Simulation(src, case["ham_name"], sd, sim_root=outroot, sim_dir=dname(case, dtau), machine="GNU", mpi=False)
    sim.run(only_prep=True)
    sf = os.path.join(outroot, dname(case, dtau), "seeds")
    with open(sf, "w") as f:
        for _ in range(8):
            f.write(f"{secrets.randbits(30) + 1:12d}\n")


def run_all(case, dtaus, outroot, nproc=3):
    procs = []
    pending = list(dtaus)
    t0 = time.time()
    while pending or procs:
        while pending and len(procs) < nproc:
            dt = pending.pop(0)
            d = os.path.join(outroot, dname(case, dt))
            p = subprocess.Popen([BIN], cwd=d, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT)
            procs.append((dt, p, time.time()))
        for item in list(procs):
            dt, p, ts = item
            if p.poll() is not None:
                print(f"  dtau={dt:g} rc={p.returncode} wall={time.time() - ts:.1f}s", flush=True)
                procs.remove(item)
        time.sleep(1.0)
    print(f"  all runs done in {time.time() - t0:.0f}s")


def autocorr1(x):
    x = np.asarray(x, float)
    x = x - x.mean()
    return float((x[:-1] * x[1:]).sum() / (x * x).sum()) if (x * x).sum() > 0 else 0.0


def analyse_dir(d, nskip):
    from py_alf.analysis import analysis
    from py_alf.ana import load_res
    import h5py
    analysis(d)
    res = load_res([d]).iloc[0]
    out = {}
    for c in res.index:
        if c.endswith("_scal0") or c.endswith("_eqR") or c.endswith("_eqR_err") or c.endswith("_scal0_err"):
            v = res[c]
            out[c] = np.asarray(v).tolist() if hasattr(v, "shape") else float(v)
    # bin autocorrelation of energy (sign-normalised bins)
    with h5py.File(os.path.join(d, "data.h5")) as f:
        for nm in ("Ener_scal", "Pot_scal", "Kin_scal", "PhiXY_scal", "Part_scal", "Constraint_scal"):
            if nm in f:
                o = f[nm]["obser"][...]
                s = f[nm]["sign"][...]
                series = (o[:, 0, 0] / s)[nskip:]  # obser: (nbin, nobs, re/im) -> real part of first element
                out[f"{nm}_bins"] = int(len(series))
                out[f"{nm}_lag1"] = autocorr1(series)
        first = "Ener_scal" if "Ener_scal" in f else "Pot_scal"
        out["signs_mean"] = float(f[first]["sign"][...][nskip:].mean())
    info = open(os.path.join(d, "info")).read()
    for key in ("Precision Green  Mean, Max", "Acceptance", "CPU Time"):
        for line in info.splitlines():
            if line.strip().startswith(key) and ":" in line:
                out["info_" + key.replace(" ", "_").replace(",", "")] = line.split(":", 1)[1].split()
    return out


def wfit(x, y, dy):
    """weighted fit y = y0 + a x^2 (x = dtau); returns y0, err"""
    A = np.vstack([np.ones_like(x), x ** 2]).T
    Wt = 1.0 / np.asarray(dy) ** 2
    cov = np.linalg.inv(A.T @ (A * Wt[:, None]))
    p = cov @ (A.T @ (Wt * y))
    chi2 = float(((y - A @ p) ** 2 * Wt).sum())
    return float(p[0]), float(np.sqrt(cov[0, 0])), chi2


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("--dtaus", default="0.2,0.1,0.05")
    ap.add_argument("--nsweep", type=int, default=200)
    ap.add_argument("--nbin", type=int, default=100)
    ap.add_argument("--tag", default="val")
    ap.add_argument("--stages", default="prep,run,ana")
    ap.add_argument("--nskip", type=int, default=5)
    ap.add_argument("--set", default="", help="override sim_dict entries, e.g. Ham_Uc=0.0,Ham_JK=0.0 (python literals)")
    ap.add_argument("--nproc", type=int, default=3, help="concurrent ALF processes (machine rule: <= 4 cores in total)")
    a = ap.parse_args()
    import ast
    case = dict(CASES[a.case])
    case["sim_dict"] = dict(case["sim_dict"])
    for kv in filter(None, a.set.split(",")):
        k, v = kv.split("=")
        case["sim_dict"][k] = ast.literal_eval(v)
    case["name"] = a.case + ("_" + a.set.replace("=", "").replace(",", "_") if a.set else "")
    dtaus = [float(x) for x in a.dtaus.split(",")]
    root = f"{W}/runs/{a.tag}"
    os.makedirs(root, exist_ok=True)
    st = a.stages.split(",")
    if "prep" in st:
        for dt in dtaus:
            d = os.path.join(root, dname(case, dt))
            if os.path.exists(d):
                raise SystemExit(f"{d} exists; refusing to overwrite raw data (use another --tag)")
            prep(case, dt, root, a.nsweep, a.nbin)
    if "run" in st:
        run_all(case, dtaus, root, a.nproc)
    if "ana" in st:
        ed = case["ed"](case["sim_dict"])
        res = {"case": a.case, "overrides": a.set, "sim_dict": case["sim_dict"], "dtaus": dtaus, "nsweep": a.nsweep, "nbin": a.nbin,
               "nskip": a.nskip, "mc": {}, "ed": {k: (np.asarray(v).tolist()) for k, v in ed.items()}}
        # pyALF reads N_skip from the parameter file: patch the namelist before analysis
        for dt in dtaus:
            d = os.path.join(root, dname(case, dt))
            ptxt = open(os.path.join(d, "parameters")).read()
            import re
            ptxt = re.sub(r"N_skip = \d+", f"N_skip = {a.nskip}", ptxt)
            open(os.path.join(d, "parameters"), "w").write(ptxt)
            res["mc"][str(dt)] = analyse_dir(d, a.nskip)
        # comparison
        table = []
        names = [k for k in ed if k.endswith("_scal0")]
        for nm in names:
            ys = np.array([res["mc"][str(dt)][nm] for dt in dtaus])
            es = np.array([res["mc"][str(dt)][nm + "_err"] for dt in dtaus])
            if np.all(es < 1e-12):  # exactly conserved quantity (e.g. particle number): no statistical error
                for dt, y in zip(dtaus, ys):
                    table.append((nm, dt, y, 0.0, 0.0 if abs(y - ed[nm]) < 1e-8 else float("inf")))
                continue
            for dt, y, e in zip(dtaus, ys, es):
                table.append((nm, dt, y, e, (y - ed[nm]) / e))
            if len(dtaus) >= 3:
                y0, e0, chi2 = wfit(np.array(dtaus), ys, es)
                table.append((nm, "extrap", y0, e0, (y0 - ed[nm]) / e0))
        print(f"{'obs':12s} {'dtau':>8s} {'MC':>14s} {'err':>10s} {'ED':>14s} {'(MC-ED)/err':>12s}")
        for nm, dt, y, e, z in table:
            print(f"{nm:12s} {dt!s:>8s} {y:14.8f} {e:10.2e} {ed[nm]:14.8f} {z:12.2f}")
        res["table"] = [(nm, dt, float(y), float(e), float(ed[nm]), float(z)) for nm, dt, y, e, z in table]
        # lattice observables: per dtau z-score of every element, then extrapolated
        lat = {}
        for nm in [k for k in ed if k.endswith("_eqR")]:
            edv = np.asarray(ed[nm]).ravel()
            vals = np.array([np.asarray(res["mc"][str(dt)][nm]).ravel() for dt in dtaus])
            errs = np.array([np.asarray(res["mc"][str(dt)][nm + "_err"]).ravel() for dt in dtaus])
            zs = []
            for j in range(edv.size):
                if np.all(errs[:, j] > 1e-12) and len(dtaus) >= 3:
                    y0, e0, _ = wfit(np.array(dtaus), vals[:, j], errs[:, j])
                    zs.append((y0 - edv[j]) / e0)
            zs_per = (vals - edv[None, :]) / np.where(errs > 1e-12, errs, np.nan)
            lat[nm] = {"extrap_z": zs, "per_dtau_z": np.nan_to_num(zs_per).tolist(), "mc_last": vals[-1].tolist(), "ed": edv.tolist()}
            print(f"{nm}: extrapolated z (all elements with err>0): {np.round(zs, 2).tolist()}")
        res["lattice"] = lat
        with open(os.path.join(root, f"summary_{case['name']}.json"), "w") as f:
            json.dump(res, f, indent=1)
        print("wrote", os.path.join(root, f"summary_{case['name']}.json"))
