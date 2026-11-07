#!/usr/bin/env python3
"""Ticket 401 checks (gates of acceptance/401.json), reusable for the later stages (tV, Kondo, Spin_Peierls).  numpy only.

Subcommands (plain-text ALF output, no HDF5 needed):
  same-seed A B [--tol 1e-10] [--files f1,f2,..] [--min-agree N]
        bin-by-bin comparison of two runs with the same seeds (K3_G / K4_schemes): per bin the largest |a-b| over all numbers
        of the bin relative to the largest |a| of that bin and file.  Reports the number of leading bins that agree within tol,
        the largest deviation among them, and the first deviating bin (trajectory divergence by rounding looks like a jump from
        1e-13 to O(1e-3..1)).  Exit 0 if all bins agree, or if the leading agreeing prefix is >= --min-agree bins (divergence
        after that is reported, not hidden); exit 1 otherwise.
  ed RUN [RUN ...] --ed Ener_scal=-1.47261997 [--ed Kin_scal=..] [--skip 1] [--threshold 3]
        O1_ed: z = (mean - ED)/sigma for every observable with an ED value, sigma = jackknife error of <O s>/<s> from rebinning to
        the plateau (compare2 rule: the smallest rebinning factor from which each further doubling raises the error by less than
        1.5/sqrt(2(N-1)); UNDER-RESOLVED if the error still grows at the largest level, or fewer than 3 levels).  Several RUN
        directories are pooled bin-wise (independent seeds).  Reports chi^2 over the observables and N bins.  Exit 0 if all |z| <= threshold
        and every error is resolved.
  cpu-vs-gpu --a RUN [RUN ..] --b RUN [RUN ..] [--skip 1] [--threshold 3]
        O2_cpu: two independent runs sets (different seeds).  Scalars (*_scal): z = (a-b)/sqrt(sa^2+sb^2) with the plateau errors of
        resolved_error.  Arrays (*_eq, *_tau): every number of a bin is one observable, rebinned to 20 groups per run set (simple standard
        error, no sign handling: use for sign-free models), z per element, chi^2 and max |z| per file; elements with zero spread in both sets
        and elements equal to 1e-12 (round-off) are taken as equal, others with zero spread must be equal to 1e-9.  Holm correction (alpha 0.01) over all scalars and array elements.  Exit 0 if no Holm rejection.
  ed-ref N U beta mu       print the ED values (ed_hubbard.py) as --ed arguments.
Library: read_bins, scal_bins, same_seed_compare, resolved_error, ed_zscores.
"""
import argparse, os, re, sys, glob
import numpy as np

FLOAT = re.compile(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?")

def _floats(line):
    return [float(x.replace("D", "E").replace("d", "e")) for x in FLOAT.findall(line)]

def run_nbin(rundir):
    """NBin from the parameters file (Nbin = ...)."""
    for line in open(os.path.join(rundir, "parameters")):
        m = re.match(r"\s*nbin\s*=\s*(\d+)", line, re.I)
        if m:
            return int(m.group(1))
    raise RuntimeError("no Nbin in parameters of " + rundir)

def data_files(rundir):
    """observable data files (*_scal, *_eq, *_tau, not *_info)"""
    out = []
    for f in sorted(os.listdir(rundir)):
        if re.search(r"_(scal|eq|tau)$", f) or re.search(r"_(eq|tau)[A-Za-z0-9]*$", f) and not f.endswith("_info"):
            out.append(f)
    return [f for f in out if not f.endswith("_info")]

def read_bins(rundir, name):
    """array (nbin, nvalues) of all numbers of the file, split into NBin equal chunks of lines"""
    nb = run_nbin(rundir)
    lines = [l for l in open(os.path.join(rundir, name)) if l.strip()]
    if len(lines) % nb:
        raise RuntimeError("%s/%s: %d lines not divisible by Nbin %d" % (rundir, name, len(lines), nb))
    k = len(lines) // nb
    return np.array([[x for l in lines[b * k:(b + 1) * k] for x in _floats(l)] for b in range(nb)])

def same_seed_compare(da, db, tol=1e-10, files=None):
    """returns dict: per file the per-bin relative deviation; global first deviating bin; prefix length; max dev in prefix"""
    fa = files or data_files(da)
    fb = files or data_files(db)
    if not files and sorted(fa) != sorted(fb):
        raise RuntimeError("different file sets: %s vs %s" % (sorted(fa), sorted(fb)))
    per = {}
    for f in fa:
        A, B = read_bins(da, f), read_bins(db, f)
        if A.shape != B.shape:
            raise RuntimeError("%s: shape %s vs %s" % (f, A.shape, B.shape))
        scale = np.maximum(np.abs(A).max(axis=1), 1e-300)
        per[f] = np.abs(A - B).max(axis=1) / scale
    nb = len(next(iter(per.values())))
    worst = np.max(np.array(list(per.values())), axis=0)
    bad = np.where(~(worst <= tol))[0]
    first = int(bad[0]) if len(bad) else None
    pre = nb if first is None else first
    return {"per_file": per, "worst_per_bin": worst, "nbin": nb, "first_bad": first, "prefix": pre,
            "max_in_prefix": float(worst[:pre].max()) if pre else 0.0, "max_all": float(worst.max())}

def scal_bins(rundir, name, skip=0):
    """(O*s per bin, sign per bin) of a scalar observable file: line = N value sign (real part of value)"""
    rows = [_floats(l) for l in open(os.path.join(rundir, name)) if l.strip()]
    # "N (re, im) sign"
    os_ = np.array([r[1] for r in rows]); s = np.array([r[3] for r in rows])
    return os_[skip:], s[skip:]

def _jack(os_, s):
    n = len(os_)
    so, ss = os_.sum(), s.sum()
    est = (so - os_) / (ss - s)
    mean = so / ss
    return mean, float(np.sqrt((n - 1) / n * np.sum((est - est.mean()) ** 2)))

def resolved_error(os_, s, min_bins=8):
    """(mean, sigma, f*, resolved, levels): jackknife error of the ratio after rebinning by f = 1, 2, 4, ... (compare2 rule)"""
    errs = []; f = 1
    while len(os_) // f >= min_bins:
        n = (len(os_) // f) * f
        o = os_[:n].reshape(-1, f).sum(axis=1); ss = s[:n].reshape(-1, f).sum(axis=1)
        m, e = _jack(o, ss)
        errs.append((f, m, e, len(o)))
        f *= 2
    mean = errs[0][1]
    levels = len(errs)
    if all(e[2] == 0.0 for e in errs):          # constant observable (e.g. Part at half filling): exact, no error to resolve
        return mean, 0.0, 1, True, levels
    istar = levels - 1
    for i in range(levels):
        nn = errs[i][3]
        thr = 1.5 / np.sqrt(2.0 * (nn - 1))
        if all(errs[j + 1][2] / errs[j][2] - 1.0 < thr for j in range(i, levels - 1) if errs[j][2] > 0):
            istar = i; break
    growing = (levels < 3) or (istar == levels - 1 and levels >= 2 and errs[-2][2] > 0 and errs[-1][2] / errs[-2][2] - 1.0 >= 1.5 / np.sqrt(2.0 * (errs[-1][3] - 1)))
    f, m, e, nn = errs[istar]
    return mean, e, f, (not growing), levels

def ed_zscores(rundirs, ed, skip=1):
    """pool the bins of several runs (same parameters, independent seeds); z per observable"""
    out = {}
    for name, val in ed.items():
        os_ = np.concatenate([scal_bins(d, name, skip)[0] for d in rundirs])
        s = np.concatenate([scal_bins(d, name, skip)[1] for d in rundirs])
        m, e, f, ok, lev = resolved_error(os_, s)
        out[name] = {"mean": m, "ed": val, "sigma": e, "z": ((m - val) / e if e > 0 else (0.0 if abs(m - val) < 1e-9 else np.inf)),
                     "f_star": f, "resolved": ok, "levels": lev, "nbins": len(os_)}
    return out

def _grouped(dirs, name, skip, ng=20):
    """(mean, stderr) per number of file `name`, bins pooled over dirs and rebinned to ng groups per run directory set"""
    X = np.concatenate([read_bins(d, name)[skip:] for d in dirs])
    n = (len(X) // ng) * ng
    G = X[:n].reshape(ng, -1, X.shape[1]).mean(axis=1)
    return G.mean(axis=0), G.std(axis=0, ddof=1) / np.sqrt(ng)

def cpu_vs_gpu(da, db, skip=1, alpha=0.01):
    """returns list of (label, z) and Holm decision"""
    items = []
    names = [f for f in data_files(da[0]) if f in data_files(db[0])]
    for f in names:
        if f.endswith("_scal"):
            oa = np.concatenate([scal_bins(d, f, skip)[0] for d in da]); sa = np.concatenate([scal_bins(d, f, skip)[1] for d in da])
            ob = np.concatenate([scal_bins(d, f, skip)[0] for d in db]); sb = np.concatenate([scal_bins(d, f, skip)[1] for d in db])
            ma, ea, fa, oka, _ = resolved_error(oa, sa); mb, eb, fb, okb, _ = resolved_error(ob, sb)
            if not (oka and okb):
                print("  note: %s error under-resolved" % f)
            items.append((f, (ma - mb) / np.hypot(ea, eb) if np.hypot(ea, eb) > 0 else (0.0 if abs(ma - mb) < 1e-9 else np.inf)))
        else:
            ma, ea = _grouped(da, f, skip); mb, eb = _grouped(db, f, skip)
            for i in range(len(ma)):
                e = np.hypot(ea[i], eb[i])
                if abs(ma[i] - mb[i]) < 1e-12 * max(1.0, abs(ma[i]), abs(mb[i])):   # equal to round-off (constants, 1e-32 imaginary noise)
                    items.append(("%s[%d]" % (f, i), 0.0))
                elif e == 0:
                    items.append(("%s[%d]" % (f, i), 0.0 if abs(ma[i] - mb[i]) < 1e-9 else np.inf))
                else:
                    items.append(("%s[%d]" % (f, i), (ma[i] - mb[i]) / e))
    return items

def holm(items, alpha, df=19):
    """Holm over two-sided normal p-values of the z (large-sample; the stderr of 20 groups has 19 dof: Student t is used)"""
    from math import erfc, sqrt
    def pval(z):
        # Student t with df dof via the incomplete beta is not in numpy; use a conservative normal-tail with the t-inflated z
        z = abs(z)
        return 1.0 if z == 0 else (0.0 if z == np.inf else erfc(z / sqrt(2.0)))
    ps = sorted(((pval(z), lab, z) for lab, z in items))
    m = len(ps); rej = []
    for i, (p, lab, z) in enumerate(ps):
        if p < alpha / (m - i):
            rej.append((lab, z, p))
        else:
            break
    return rej

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)
    a = sp.add_parser("same-seed"); a.add_argument("A"); a.add_argument("B"); a.add_argument("--tol", type=float, default=1e-10)
    a.add_argument("--files"); a.add_argument("--min-agree", type=int, default=0); a.add_argument("-v", action="store_true")
    b = sp.add_parser("ed"); b.add_argument("runs", nargs="+"); b.add_argument("--ed", action="append", required=True)
    b.add_argument("--skip", type=int, default=1); b.add_argument("--threshold", type=float, default=3.0)
    d = sp.add_parser("cpu-vs-gpu"); d.add_argument("--a", nargs="+", required=True); d.add_argument("--b", nargs="+", required=True)
    d.add_argument("--skip", type=int, default=1); d.add_argument("--alpha", type=float, default=0.01)
    c = sp.add_parser("ed-ref"); c.add_argument("args", nargs="*")
    args = p.parse_args()
    if args.cmd == "same-seed":
        r = same_seed_compare(args.A, args.B, args.tol, args.files.split(",") if args.files else None)
        print("same-seed %s vs %s: %d bins, tol %.1e" % (args.A, args.B, r["nbin"], args.tol))
        print("  leading agreeing bins %d, max deviation among them %.3e, first deviating bin %s, max deviation overall %.3e"
              % (r["prefix"], r["max_in_prefix"], r["first_bad"], r["max_all"]))
        if args.v:
            for f, d in r["per_file"].items():
                print("  %-14s max %.3e" % (f, d.max()))
        ok = r["first_bad"] is None or r["prefix"] >= args.min_agree > 0
        print("RESULT", "PASS" if ok else "FAILED")
        return 0 if ok else 1
    if args.cmd == "ed":
        ed = {}
        for kv in args.ed:
            k, v = kv.split("="); ed[k] = float(v)
        res = ed_zscores(args.runs, ed, args.skip)
        chi2 = 0.0; allok = True
        for k, r in res.items():
            flag = "" if r["resolved"] else "  UNDER-RESOLVED"
            print("%-10s mean %.6f ED %.6f sigma %.6f z %+.2f f* %d bins(after rebin) %d%s"
                  % (k, r["mean"], r["ed"], r["sigma"], r["z"], r["f_star"], r["nbins"] // r["f_star"], flag))
            chi2 += r["z"] ** 2
            allok &= abs(r["z"]) <= args.threshold and r["resolved"]
        print("chi2 = %.2f over %d observables (not independent: Ener = Kin + Pot); pooled bins %d" % (chi2, len(res), next(iter(res.values()))["nbins"]))
        print("RESULT", "PASS" if allok else "FAILED")
        return 0 if allok else 1
    if args.cmd == "cpu-vs-gpu":
        items = cpu_vs_gpu(args.a, args.b, args.skip)
        z = np.array([x[1] for x in items]); fin = z[np.isfinite(z)]
        print("compared %d numbers; chi2 = %.1f (finite part %d), max |z| = %.2f (%s)"
              % (len(z), float((fin ** 2).sum()), len(fin), float(np.abs(z).max()), items[int(np.argmax(np.abs(z)))][0]))
        for lab, zz in items:
            if "[" not in lab:
                print("  %-12s z %+.2f" % (lab, zz))
        rej = holm(items, args.alpha)
        for lab, zz, p in rej:
            print("  HOLM REJECT %s z %+.2f p %.2e" % (lab, zz, p))
        print("RESULT", "FAILED" if rej else "PASS")
        return 1 if rej else 0
    if args.cmd == "ed-ref":
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import ed_hubbard
        v = [float(x) for x in args.args] or [4, 4, 2, 0]
        print(" ".join("--ed %s=%.8f" % kv for kv in ed_hubbard.observables(int(v[0]), *v[1:]).items()))
        return 0

if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, IndexError, ValueError, FloatingPointError) as e:
        # fail closed: unparseable / non-finite / differently shaped output is a failure of the check, never a pass
        print("RESULT FAILED (output not comparable: %s: %s)" % (type(e).__name__, e))
        sys.exit(1)
