#!/usr/bin/env python3
"""Statistical comparison of two ALF run directories (independent chains, different seeds).

Usage: python3 -I stat_compare.py DIR_A DIR_B      (A = reference, e.g. standard; B = candidate)
       python3 -I stat_compare.py --selftest

Reads the analysed outputs (ana.out, run by stat_runs.sh):
  *_scalJ   observable rows "idx value error" (jackknife, sign-reweighted); the LAST row of each file is the
            average sign/phase. Every observable row of every *_scalJ present in both directories is a scalar
            (Ener_scal = energy, Pot_scal = double-occupancy proxy, Kin_scal, Part_scal); the sign is compared once
            (Ener_scal's) and checked for consistency with the others.
  *_eqJK / *_eqJR   k- / r-resolved equal-time correlations; entries "i j Re errRe Im errIm" under each
            coordinate line; Re and Im are separate entries (SpinZ_eq, Green_eq).
z = (a-b)/sqrt(sa^2+sb^2).
PASS: every scalar |z| < 3; for each k-point array (*_eqJK): fraction(|z|>2) <= 10 % and max|z| < 4.5.
(*_eqJR, the Fourier transform of the same data, is printed for information and does not decide.)
Entries that are numerically zero in both runs (|value| < 1e-10 and error < 1e-10, e.g. Im parts of
symmetric correlations) are skipped from the array statistics (their count is reported).
Also reports 'Precision Green/Phase' lines of both info files; flags (does not fail) if B's max > 10x A's.
Exit 0 PASS, 1 FAIL, 2 unusable input.
"""
import glob
import math
import os
import re
import sys

Z_SCALAR = 3.0
Z_ARR_FRAC = 2.0
FRAC_MAX = 0.10
Z_ARR_MAX = 4.5
ZERO = 1e-10
PREC_FLAG = 10.0


def floats(line):
    return [float(x.replace("D", "E")) for x in line.split()]


def parse_scalJ(path):
    """-> list of (idx, value, error); also number of bins line is ignored."""
    rows = []
    for line in open(path):
        t = line.split()
        if len(t) == 3 and t[0].isdigit():
            rows.append((int(t[0]), float(t[1]), float(t[2])))
    return rows


def parse_eqJ(path):
    """-> dict (coord tuple, i, j, 'Re'|'Im') -> (value, error)."""
    out = {}
    coord = None
    for line in open(path):
        t = line.split()
        if len(t) == 2:
            coord = (float(t[0]), float(t[1]))
        elif len(t) == 6 and coord is not None:
            i, j = int(t[0]), int(t[1])
            re_v, re_e, im_v, im_e = map(float, t[2:])
            out[(coord, i, j, "Re")] = (re_v, re_e)
            out[(coord, i, j, "Im")] = (im_v, im_e)
    return out


def zscore(a, ea, b, eb):
    """z = (a-b)/sqrt(ea^2+eb^2). A quantity that is constant to roundoff in both runs (errors and
    difference below 1e-9 relative, e.g. the particle number at half filling) has z = 0; a difference
    larger than that with zero/roundoff errors gives |z| = inf."""
    s = math.hypot(ea, eb)
    d = a - b
    tol = 1e-9 * max(1.0, abs(a), abs(b))
    if s <= tol:
        return 0.0 if abs(d) <= tol else math.copysign(math.inf, d)
    return d / s


def parse_precision(path):
    """-> {'Green mean','Green max','Phase max',...} from the info file."""
    res = {}
    for line in open(path):
        m = re.search(r"Precision\s+Green\s+Mean,\s*Max\s*:\s*(\S+)\s+(\S+)", line)
        if m:
            res["Green mean"], res["Green max"] = float(m.group(1)), float(m.group(2))
        m = re.search(r"Precision\s+Phase,\s*Max\s*:\s*(\S+)", line)
        if m:
            res["Phase max"] = float(m.group(1))
        m = re.search(r"Precision\s+tau\s+Mean,\s*Max\s*:\s*(\S+)\s+(\S+)", line, re.I)
        if m:
            res["tau mean"], res["tau max"] = float(m.group(1)), float(m.group(2))
    return res


def compare(da, db, out=print):
    ok = True
    # scalars
    names_a = {os.path.basename(p) for p in glob.glob(os.path.join(da, "*_scalJ"))}
    names_b = {os.path.basename(p) for p in glob.glob(os.path.join(db, "*_scalJ"))}
    common = sorted(names_a & names_b)
    if not common:
        out("no common *_scalJ files (run ana.out first)")
        return 2
    if names_a ^ names_b:
        out("warning: scalar files only in one directory: %s" % sorted(names_a ^ names_b))
    out("%-22s %14s %10s %14s %10s %9s  %s" % ("scalar", "A", "errA", "B", "errB", "z", ""))
    signs = []
    for f in common:
        ra, rb = parse_scalJ(os.path.join(da, f)), parse_scalJ(os.path.join(db, f))
        if len(ra) != len(rb) or len(ra) < 2:
            out("%s: row count mismatch or <2 rows (%d vs %d)" % (f, len(ra), len(rb)))
            return 2
        base = f[:-len("_scalJ")]
        for k in range(len(ra) - 1):
            z = zscore(ra[k][1], ra[k][2], rb[k][1], rb[k][2])
            bad = not abs(z) < Z_SCALAR
            ok &= not bad
            out("%-22s %14.6g %10.3g %14.6g %10.3g %9.2f  %s" % (
                "%s[%d]" % (base, ra[k][0]), ra[k][1], ra[k][2], rb[k][1], rb[k][2], z, "FAIL" if bad else "ok"))
        signs.append((base, ra[-1], rb[-1]))
    # sign / phase: Ener_scal's, others must agree with it
    ref = [s for s in signs if s[0] == "Ener_scal"] or signs[:1]
    base, sa, sb = ref[0]
    z = zscore(sa[1], sa[2], sb[1], sb[2])
    bad = not abs(z) < Z_SCALAR
    ok &= not bad
    out("%-22s %14.6g %10.3g %14.6g %10.3g %9.2f  %s" % (
        "avg sign/phase", sa[1], sa[2], sb[1], sb[2], z, "FAIL" if bad else "ok"))
    for b2, a2, c2 in signs:
        if abs(a2[1] - sa[1]) > 1e-10 or abs(c2[1] - sb[1]) > 1e-10:
            out("warning: sign in %s_scalJ differs from %s_scalJ" % (b2, base))

    # arrays
    arr_files = []
    for pat in ("*_eqJK", "*_eqJR"):
        na = {os.path.basename(p) for p in glob.glob(os.path.join(da, pat))}
        nb = {os.path.basename(p) for p in glob.glob(os.path.join(db, pat))}
        if na ^ nb:
            out("warning: array files only in one directory: %s" % sorted(na ^ nb))
        arr_files += sorted(na & nb)
    if not arr_files:
        out("no common *_eqJK/*_eqJR files")
        return 2
    out("")
    out("%-12s %7s %8s %8s %7s %6s | %s" % ("array", "entries", "skipped0", "frac>2", "max|z|", "mean z",
                                          "informational: symmetry-unique entries, frac>2"))
    for f in arr_files:
        A, B = parse_eqJ(os.path.join(da, f)), parse_eqJ(os.path.join(db, f))
        if set(A) != set(B):
            out("%s: entry sets differ" % f)
            return 2
        zs, skipped = [], 0
        for key in A:
            (a, ea), (b, eb) = A[key], B[key]
            if max(abs(a), abs(b)) < ZERO and max(ea, eb) < ZERO:
                skipped += 1
                continue
            zs.append(zscore(a, ea, b, eb))
        n = len(zs)
        if n == 0:
            out("%s: no usable entries" % f)
            return 2
        frac = sum(1 for z in zs if abs(z) > Z_ARR_FRAC) / n
        mx = max(abs(z) for z in zs)
        # lattice-symmetry copies of one entry have identical z (up to sign): count each |z| once
        uz = sorted({round(abs(z), 5) for z in zs})
        ufrac = sum(1 for z in uz if z > Z_ARR_FRAC) / len(uz)
        # only the k-point arrays (*_eqJK) decide the verdict; *_eqJR is the Fourier transform (info only)
        decides = f.endswith("JK")
        bad = decides and not (frac <= FRAC_MAX and mx < Z_ARR_MAX)
        ok &= not bad
        out("%-12s %7d %8d %7.1f%% %7.2f %6.2f | %d unique, %.1f%%  %s" % (
            f, n, skipped, 100 * frac, mx, sum(zs) / n, len(uz), 100 * ufrac,
            "FAIL" if bad else ("ok" if decides else "(info)")))

    # precision diagnostics
    out("")
    try:
        pa, pb = parse_precision(os.path.join(da, "info")), parse_precision(os.path.join(db, "info"))
    except OSError as e:
        out("info not readable: %s" % e)
        pa = pb = {}
    out("%-14s %12s %12s %8s  %s" % ("precision", "A", "B", "B/A", ""))
    for k in sorted(set(pa) & set(pb)):
        r = pb[k] / pa[k] if pa[k] > 0 else (math.inf if pb[k] > 0 else 1.0)
        flag = "FLAG (B max > %gx A)" % PREC_FLAG if k.endswith("max") and r > PREC_FLAG else ""
        out("%-14s %12.3e %12.3e %8.2f  %s" % (k, pa[k], pb[k], r, flag))
    out("")
    out("VERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


def selftest():
    assert zscore(1.0, 0.1, 1.0, 0.1) == 0.0
    assert abs(zscore(2.0, 0.3, 1.0, 0.4) - 2.0) < 1e-12
    assert zscore(1.0, 0.0, 1.0, 0.0) == 0.0 and zscore(1.0, 0.0, 2.0, 0.0) == -math.inf
    assert zscore(36.0, 3e-13, 36.0, 1e-12) == 0.0 and zscore(36.0, 3e-13, 36.001, 1e-12) == -math.inf
    import tempfile
    t = tempfile.mkdtemp()
    p = os.path.join(t, "x_scalJ")
    open(p, "w").write(" Effective number of bins, and bins:  5  6\n\n 1 -2.5E+01 1.0E-01\n\n 2 1.0E+00 0.0E+00\n")
    assert parse_scalJ(p) == [(1, -25.0, 0.1), (2, 1.0, 0.0)]
    p = os.path.join(t, "y_eqJK")
    open(p, "w").write(" 0.5 -0.5\n 1 1 1.0E+00 1.0E-01 2.0E-17 1.0E-17\n")
    d = parse_eqJ(p)
    assert d[((0.5, -0.5), 1, 1, "Re")] == (1.0, 0.1)
    print("selftest ok")
    return 0


def main(argv):
    if argv[1:] == ["--selftest"]:
        return selftest()
    if len(argv) != 3:
        sys.stderr.write(__doc__)
        return 2
    for d in argv[1:]:
        if not os.path.isdir(d):
            sys.stderr.write("not a directory: %s\n" % d)
            return 2
    return compare(argv[1], argv[2])


if __name__ == "__main__":
    sys.exit(main(sys.argv))
