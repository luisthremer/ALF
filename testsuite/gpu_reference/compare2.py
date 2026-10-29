#!/usr/bin/env python3
"""Statistical comparison of two ALF run directories (CPU vs GPU, compiler vs compiler).

usage: compare2.py RUNDIR_A RUNDIR_B [--alpha 0.01] [--json out.json] [--min-bins 8]

Reads the per-bin data in RUNDIR/data.h5 (HDF5 build, one value per bin, already averaged over
all MPI ranks) through pyALF's bare-bin readers. Replaces compare.py as the observable-level gate.

METHOD
 1. Curated observables (a small, nearly independent set; scalar per bin):
      E, Kin, Pot, Part                               (scalar observables)
      S_zz(k) at (pi,pi), (pi/2,pi/2)                 (SpinZ_eq)
      Den(k) at (pi,pi)                               (Den_eq)
      G(k) at (0,0), (pi/2,pi/2), (pi,0)              (Green_eq)
      G(k=(pi,0), tau=beta/2), S_zz((pi,pi), tau=beta/2)   (Green_tau, SpinZ_tau; needs Ltau=1)
      chi_zz(pi,pi) = int_0^beta dtau S_zz((pi,pi),tau)    (trapezoid over the tau grid)
    k-points are mapped through the lattice stored in data.h5 (py_alf Lattice.k_to_n), never by
    assuming an index order. Points that do not exist on the lattice (e.g. pi/2 on 6x6) are
    skipped and listed. Observables with zero variance in both runs (e.g. Part at half filling)
    are compared exactly (relative 1e-9) and do not count towards the test family.
 2. Errors: the first n_skip bins (parameters file) are dropped; the remaining bins are rebinned
    by factors 1,2,4,... as long as at least --min-bins (8) bins remain. The error at each factor
    is the jackknife error of the ratio <O s>/<s> (sign-aware). The plateau factor f* is the
    smallest factor from which every further doubling raises the error by less than
    1.5/sqrt(2(N-1)) (relative; this is ~1.5 times the statistical noise of an error estimate from
    N bins). If the error is still growing at the largest factor the observable is flagged
    UNDER-RESOLVED (autocorrelation time comparable to the run length: errors are underestimated);
    sigma then is the error at the largest factor. A plateau also cannot be demonstrated, and the
    observable is flagged the same way, if fewer than 3 rebin levels exist (< 4*min_bins bins). tau_int = 0.5 (sigma_f*/sigma_1)^2 (in bins).
 3. Test per observable: t = (a-b)/sqrt(sa^2+sb^2), Welch-Satterthwaite degrees of freedom from
    (N_a(f*)-1, N_b(f*)-1), two-sided Student-t p-value. Family-wise false-alarm rate alpha over the
    m tested observables via Sidak: reject if p < 1-(1-alpha)^(1/m). Sum t^2 vs chi2_m is shown for
    information only.
 4. Informational, never gating: per-array chi2/n of all k (and k,tau) entries at rebin factor 1,
    i.e. the old compare.py statistic.

EXIT CODE  0: PASS, i.e. no observable differs and all errors are resolved (or all means agree to 1e-9 relative,
 e.g. same-seed runs).  1: at least one resolved observable differs (FAIL).
 2: INCONCLUSIVE: no resolved observable differs, but at least one observable is UNDER-RESOLVED (its error
 cannot be trusted; run longer / more bins). Not green.
 3: INPUT ERROR (fail-closed, 2026-10-29): a non-finite value in any bin of either run; no tested observable
 (m = 0); fewer tested observables than --expect-m; an observable named in --expect missing; identical `seeds`
 files without --same-seed (or a missing `seeds` file); identical means without --same-seed. Not green.
"""
import argparse
import json
import math
import os
import sys

import numpy as np
from scipy import stats

from py_alf.ana import Parameters, read_latt, read_scal

K_PI, K_HALF, K_0 = "(pi,pi)", "(pi/2,pi/2)", "(0,0)"


# ------------------------------------------------------------------ data loading
def load_run(rundir):
    """Return dict: name -> (X[bins], S[bins]) per-bin series (skip not yet applied), plus info."""
    if not os.path.exists(os.path.join(rundir, "data.h5")):
        raise SystemExit(f"{rundir}: no data.h5 (HDF5 build required)")
    nskip = Parameters(rundir).N_skip()
    series, notes, arrays = {}, [], {}

    for label, grp in (("E", "Ener_scal"), ("Kin", "Kin_scal"), ("Pot", "Pot_scal"),
                       ("Part", "Part_scal")):
        obs, sign, _ = read_scal(rundir, grp, bare_bins=True)
        series[label] = (obs[:, 0].real, sign)

    def kidx(latt, kind):
        pt = {K_PI: (np.pi, np.pi), K_HALF: (np.pi / 2, np.pi / 2), K_0: (0.0, 0.0),
              "(pi,0)": (np.pi, 0.0)}[kind]
        try:
            return latt.k_to_n(pt)
        except Exception:  # point not on this finite lattice
            return None

    def add(label, key, x, sign):
        series[label] = (x, sign)

    meta = {}
    for grp, kinds in (("SpinZ_eq", (K_PI, K_HALF)), ("Den_eq", (K_PI,)),
                       ("Green_eq", (K_0, K_HALF, "(pi,0)"))):
        obs, _, sign, _, _, _, latt = read_latt(rundir, grp, bare_bins=True)
        tr = obs.trace(axis1=1, axis2=2)  # (bins, nt=1, N): trace over orbitals
        arrays[grp] = (tr[:, 0, :].real, sign)
        meta["N"] = int(latt.N)
        short = {"SpinZ_eq": "Szz", "Den_eq": "Den", "Green_eq": "G"}[grp]
        for kind in kinds:
            n = kidx(latt, kind)
            if n is None:
                notes.append(f"skipped {short}(k={kind}): not on lattice")
                continue
            series[f"{short}(k={kind})"] = (tr[:, 0, n].real, sign)

    have_tau = True
    for grp in ("Green_tau", "SpinZ_tau"):
        try:
            obs, _, sign, _, ntau, dtau, latt = read_latt(rundir, grp, bare_bins=True)
        except (KeyError, OSError):
            have_tau = False
            notes.append(f"{grp} missing (Ltau=0?): tau observables skipped")
            continue
        tr = obs.trace(axis1=1, axis2=2)  # (bins, ntau, N)
        arrays[grp] = (tr.real, sign)
        meta["ntau"], meta["dtau"] = int(ntau), float(dtau)
        if (ntau - 1) % 2:
            notes.append(f"{grp}: odd number of tau intervals, beta/2 not on grid")
            continue
        half = (ntau - 1) // 2
        if grp == "Green_tau":
            n = kidx(latt, "(pi,0)")
            if n is not None:
                series["G(k=(pi,0),tau=beta/2)"] = (tr[:, half, n].real, sign)
        else:
            n = kidx(latt, K_PI)
            if n is not None:
                series["Szz(k=(pi,pi),tau=beta/2)"] = (tr[:, half, n].real, sign)
                series["chi_zz(pi,pi)"] = (np.trapezoid(tr[:, :, n].real, dx=dtau, axis=1), sign)
    meta["beta"] = (meta.get("ntau", 1) - 1) * meta.get("dtau", 0.0)
    return {"series": series, "arrays": arrays, "nskip": nskip, "notes": notes, "meta": meta,
            "dir": rundir}


# ------------------------------------------------------------------ error estimation
def rebin_ratio(X, S, f, nskip):
    """Skip nskip bins, rebin by f (drop leftover bins at the start), return grouped X, S."""
    X, S = X[nskip:], S[nskip:]
    n = len(X) // f
    X, S = X[len(X) - n * f:], S[len(S) - n * f:]
    shp = (n, f) + X.shape[1:]
    return X.reshape(shp).mean(axis=1), S.reshape((n, f)).mean(axis=1)


def jackknife_ratio(X, S):
    """Mean and jackknife error of <X>/<S> over the first axis of X (vectorised in the rest)."""
    n = len(X)
    sx, ss = X.sum(axis=0), S.sum()
    sshape = (n,) + (1,) * (X.ndim - 1)
    J = (sx - X) / (ss - S).reshape(sshape)  # leave-one-out ratio
    mean = sx / ss
    err = np.sqrt((n - 1) * np.mean((J - J.mean(axis=0)) ** 2, axis=0))
    return mean, err


def error_levels(x, s, nskip, min_bins):
    """Error at rebin factors 1,2,4,... with >= min_bins bins. Returns list of (f, N, mean, err)."""
    out, f = [], 1
    while True:
        n = (len(x) - nskip) // f
        if n < min_bins and f > 1:
            break
        if n < 2:
            break
        xr, sr = rebin_ratio(x, s, f, nskip)
        m, e = jackknife_ratio(xr, sr)
        out.append((f, len(xr), float(m), float(e)))
        f *= 2
    return out


def plateau(levels):
    """Choose plateau level. Returns (index, under_resolved, tau_int)."""
    k = len(levels) - 1
    if k < 2:  # fewer than 3 rebin levels: a plateau cannot be demonstrated (needs >= 4*min_bins bins)
        e1 = levels[0][3]
        return k, True, (0.5 * (levels[k][3] / e1) ** 2 if e1 > 0 else float("nan"))

    def grows(j):
        e0, e1, n1 = levels[j][3], levels[j + 1][3], levels[j + 1][1]
        return e1 > e0 * (1 + 1.5 / math.sqrt(2 * (n1 - 1)))

    idx = None
    for i in range(k):  # smallest i from which no further doubling grows the error
        if not any(grows(j) for j in range(i, k)):
            idx = i
            break
    under = idx is None
    if under:
        idx = k
    e1 = levels[0][3]
    tau = 0.5 * (levels[idx][3] / e1) ** 2 if e1 > 0 else float("nan")
    return idx, under, tau


def summarize_obs(x, s, nskip, min_bins):
    lev = error_levels(x, s, nskip, min_bins)
    idx, under, tau = plateau(lev)
    f, n, m, e = lev[idx]
    return {"mean": m, "err": e, "f": f, "N": n, "under": under, "tau_int": tau,
            "levels": [(a, b, d) for a, b, c, d in lev]}


# ------------------------------------------------------------------ comparison
def welch_t(a, b):
    va, vb = a["err"] ** 2, b["err"] ** 2
    s2 = va + vb
    t = (a["mean"] - b["mean"]) / math.sqrt(s2)
    dfa, dfb = max(a["N"] - 1, 1), max(b["N"] - 1, 1)
    nu = s2 ** 2 / (va ** 2 / dfa + vb ** 2 / dfb)
    p = 2 * stats.t.sf(abs(t), nu)
    return t, nu, p


def full_array_chi2(A, B, nskip_a, nskip_b):
    """Informational: per-array mean chi2/n of all entries at rebin 1 (old compare.py statistic)."""
    out = {}
    for grp in sorted(set(A["arrays"]) & set(B["arrays"])):
        xa, sa = A["arrays"][grp]
        xb, sb = B["arrays"][grp]
        if xa.shape[1:] != xb.shape[1:]:
            continue
        ma, ea = jackknife_ratio(xa[nskip_a:], sa[nskip_a:])
        mb, eb = jackknife_ratio(xb[nskip_b:], sb[nskip_b:])
        sig = np.sqrt(ea ** 2 + eb ** 2)
        ok = sig > 0
        z = (ma - mb)[ok] / sig[ok]
        out[grp] = {"n": int(z.size), "chi2_per_n": float(np.mean(z ** 2)) if z.size else None,
                    "max_abs_z": float(np.abs(z).max()) if z.size else None}
    return out


def _is_const(x, rel=1e-9):
    """True if the per-bin series is constant up to round-off (e.g. density at half filling)."""
    return bool(np.ptp(x) <= rel * max(abs(x).max(), 1e-300))


class InputError(Exception):
    pass


def check_inputs(A, B, expect, expect_m, same_seed):
    """Fail-closed input checks (exit 3). Returns the list of observable names to test."""
    for run in (A, B):
        for k, (x, s) in list(run["series"].items()) + list(run["arrays"].items()):
            if not (np.all(np.isfinite(x)) and np.all(np.isfinite(s))):
                raise InputError(f"{run['dir']}: non-finite value in {k}")
    only_a = sorted(set(A["series"]) - set(B["series"]))
    only_b = sorted(set(B["series"]) - set(A["series"]))
    if only_a or only_b:
        raise InputError(f"observables present in one run only: A {only_a} B {only_b}")
    names = [k for k in A["series"] if k in B["series"]]
    missing = [e for e in expect if e not in names]
    if missing:
        raise InputError(f"expected observables missing: {missing}")
    seeds = [os.path.join(r["dir"], "seeds") for r in (A, B)]
    if not same_seed:
        if not all(os.path.isfile(f) for f in seeds):
            raise InputError("a seeds file is missing (needed to prove independent seeds)")
        if open(seeds[0], "rb").read() == open(seeds[1], "rb").read():
            raise InputError("identical seeds files: this test needs independent seeds (use --same-seed otherwise)")
    return names


def compare(A, B, alpha, min_bins, names=None, same_seed=False):
    names = names if names is not None else [k for k in A["series"] if k in B["series"]]
    rows, exact = [], []
    for k in names:
        xa, sa = A["series"][k]
        xb, sb = B["series"][k]
        if _is_const(xa[A["nskip"]:]) and _is_const(xb[B["nskip"]:]):  # constant (symmetry) observable
            ma = xa[A["nskip"]:].mean() / sa[A["nskip"]:].mean()
            mb = xb[B["nskip"]:].mean() / sb[B["nskip"]:].mean()
            exact.append({"name": k, "a": float(ma), "b": float(mb),
                          "ok": bool(abs(ma - mb) <= 1e-9 * max(1.0, abs(ma), abs(mb)))})
            continue
        ra = summarize_obs(xa, sa, A["nskip"], min_bins)
        rb = summarize_obs(xb, sb, B["nskip"], min_bins)
        t, nu, p = welch_t(ra, rb)
        rows.append({"name": k, "a": ra, "b": rb, "t": t, "nu": nu, "p": p})
    m = max(len(rows), 1)
    p_crit = 1 - (1 - alpha) ** (1.0 / m)
    for r in rows:
        r["differs"] = bool(r["p"] < p_crit)
        r["under"] = bool(r["a"]["under"] or r["b"]["under"])
        r["t_crit"] = float(stats.t.isf(p_crit / 2, r["nu"]))
    sumt2 = float(sum(r["t"] ** 2 for r in rows))
    bad_exact = [e for e in exact if not e["ok"]]
    fail = any(r["differs"] and not r["under"] for r in rows) or bool(bad_exact)
    # identical means (same-seed runs with identical trajectories; round-off ~1e-12) need no error estimate
    identical = bool(rows) and all(abs(r["a"]["mean"] - r["b"]["mean"])
                    <= 1e-9 * max(1.0, abs(r["a"]["mean"]), abs(r["b"]["mean"])) for r in rows)
    if identical and not same_seed:
        raise InputError("identical means on independent seeds: the two runs are probably the same run")
    # a PASS needs trustworthy errors: any under-resolved observable makes the result INCONCLUSIVE
    incon = (not fail) and not identical and any(r["differs"] or r["under"] for r in rows)
    return {"rows": rows, "exact": exact, "m": len(rows), "alpha": alpha, "p_crit": p_crit,
            "sum_t2": sumt2, "chi2_p": float(stats.chi2.sf(sumt2, max(len(rows), 1))),
            "verdict": "FAIL" if fail else ("INCONCLUSIVE" if incon else "PASS")}


def report(res, A, B, arrays):
    print(f"A = {A['dir']}\nB = {B['dir']}")
    print(f"family: m={res['m']} observables, alpha={res['alpha']} (Sidak per-test p < {res['p_crit']:.2e})")
    for run in (A, B):
        for n in run["notes"]:
            print(f"note [{os.path.basename(os.path.normpath(run['dir']))}]: {n}")
    hdr = f"{'observable':28s} {'A':>22s} {'B':>22s} {'t':>7s} {'nu':>5s} {'tcrit':>6s} {'p':>9s} {'tauA':>5s} {'tauB':>5s}"
    print(hdr)
    for r in res["rows"]:
        a, b = r["a"], r["b"]
        flag = "DIFF" if r["differs"] else "ok"
        if r["under"]:
            flag += " UNDER-RESOLVED"
        print(f"{r['name']:28s} {a['mean']:12.6g}+-{a['err']:<8.2g} {b['mean']:12.6g}+-{b['err']:<8.2g} "
              f"{r['t']:7.2f} {r['nu']:5.1f} {r['t_crit']:6.2f} {r['p']:9.2e} "
              f"{a['tau_int']:5.1f} {b['tau_int']:5.1f}  {flag}  (f*={a['f']}/{b['f']}, N={a['N']}/{b['N']})")
    for e in res["exact"]:
        print(f"{e['name']:28s} constant: {e['a']:.12g} vs {e['b']:.12g}  {'ok' if e['ok'] else 'DIFF'}")
    print(f"info: sum t^2 = {res['sum_t2']:.1f} vs chi2_{res['m']} (p = {res['chi2_p']:.3f})")
    print("info: full-array chi2/n at rebin 1 (old statistic, not gating):")
    for g, v in arrays.items():
        print(f"   {g:12s} n={v['n']:6d} chi2/n={v['chi2_per_n']:.2f} max|z|={v['max_abs_z']:.2f}")
    print("RESULT:", res["verdict"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--alpha", type=float, default=0.01, help="family-wise false-alarm rate")
    ap.add_argument("--min-bins", type=int, default=8, help="minimum effective bins after rebinning")
    ap.add_argument("--json", help="write full results here")
    ap.add_argument("--expect", default="", help="semicolon-separated observables that must be tested (missing: exit 3)")
    ap.add_argument("--expect-m", type=int, default=1, help="minimum number of tested observables (fewer: exit 3)")
    ap.add_argument("--expect-m-exact", type=int, default=None, help="exact number of tested observables (else exit 3)")
    ap.add_argument("--same-seed", action="store_true", help="the two runs share seeds (allows identical means)")
    args = ap.parse_args(argv)
    expect = [e for e in args.expect.split(";") if e]
    A = B = {"dir": "?"}
    try:
        try:
            A, B = load_run(args.a), load_run(args.b)
        except (SystemExit, KeyError, OSError, ValueError, IndexError) as e:  # unreadable input is an input error
            raise InputError(f"cannot load the runs: {e}")
        names = check_inputs(A, B, expect, args.expect_m, args.same_seed)
        res = compare(A, B, args.alpha, args.min_bins, names, args.same_seed)
        if res["m"] < max(args.expect_m, 1):
            raise InputError(f"only {res['m']} observables tested, expected >= {max(args.expect_m, 1)}")
        if args.expect_m_exact is not None and res["m"] != args.expect_m_exact:
            raise InputError(f"{res['m']} observables tested, expected exactly {args.expect_m_exact}")
        if not all(math.isfinite(v) for r in res["rows"] for v in (r["t"], r["p"], r["a"]["err"], r["b"]["err"])):
            raise InputError("non-finite test statistic")
    except InputError as e:
        print(f"A = {A['dir']}\nB = {B['dir']}")
        print("INPUT ERROR:", e)
        print("RESULT: INPUT_ERROR")
        if args.json:
            with open(args.json, "w", encoding="utf-8") as f:
                json.dump({"verdict": "INPUT_ERROR", "reason": str(e)}, f, indent=1)
        return 3
    arrays = full_array_chi2(A, B, A["nskip"], B["nskip"])
    res["arrays"] = arrays
    report(res, A, B, arrays)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=1, default=float)
    return {"PASS": 0, "FAIL": 1, "INCONCLUSIVE": 2}[res["verdict"]]


if __name__ == "__main__":
    sys.exit(main())
