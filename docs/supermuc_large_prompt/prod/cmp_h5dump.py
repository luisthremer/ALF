#!/usr/bin/env python3
"""Compare two ALF HDF5 files (data.h5 or confout_*.h5) numerically through `h5dump` text (no h5py needed).

usage: cmp_h5dump.py A.h5 B.h5 [--tol 1e-10] [--h5dump PATH]     A = reference, B = restarted / recovered run
       (a file ending in .txt is taken as an existing `h5dump -m %.17g` text dump)

Bin datasets (maximum first dimension H5S_UNLIMITED, i.e. <group>/obser, /sign, /back; first dimension = bin index) are compared bin by bin.
  * "duplicated bin": two consecutive bins identical in ALL bin datasets (what a kill between the data.h5 write and the confout write
    causes after the restart). The later copy is dropped before comparing; the positions are reported.
  * raw number of bins of B must equal that of A (the restart is meant to deliver the same number of bins, one of them possibly a copy).
  * after dropping duplicates, the first min(nA, nB') bins are compared.
Other datasets (confout: configuration, seed; lattice datasets) must agree in all values. Attributes (e.g. the parameter groups store NBin)
are listed when they differ, but do not decide the verdict.
Deviation: rel = |a-b| / max(|a|,|b|) elementwise (0 if both 0); scaled = |a-b| / max|a| of the dataset. Verdict PASS if max scaled <= tol
(and rel is reported). With bit-identical runs both are 0.
Machine-readable last lines:  DUPS=<n>  NBINS_A=<n> NBINS_B=<n>  MAXREL=<x> MAXSCALED=<x>  RESULT=PASS|FAIL
exit status 0 = PASS, 1 = FAIL, 2 = usage / dump error."""
import os, re, subprocess, sys


def dump(path, h5dump):
    if path.endswith(".txt"):
        return open(path).read()
    p = subprocess.Popen([h5dump, "-m", "%.17g", path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    out, err = p.communicate()
    if p.returncode != 0:
        sys.stderr.write("h5dump failed on %s: %s\n" % (path, err.strip()[:300]))
        sys.exit(2)
    return out


def tok(t):
    t = t.strip()
    try:
        return float(t)
    except ValueError:
        return t.strip('"')


def parse(text):
    """-> (datasets, attrs); datasets[path] = {"dims": [...], "unlimited": bool, "vals": {first_index: [values]}}"""
    ds, attrs = {}, {}
    stack = []            # (kind, name)
    cur = None            # current dataset/attribute record
    indata = False
    for raw in text.splitlines():
        line = raw.strip()
        m = re.match(r'^(GROUP|DATASET|ATTRIBUTE)\s+"([^"]*)"\s*\{$', line)
        if m:
            kind, name = m.group(1), m.group(2)
            stack.append((kind, name))
            if kind in ("DATASET", "ATTRIBUTE"):
                groups = [n for k, n in stack if k == "GROUP" and n != "/"]
                dsname = [n for k, n in stack if k == "DATASET"]
                path = "/" + "/".join(groups + dsname[:1])
                if kind == "ATTRIBUTE":
                    path = path + "@" + name
                cur = {"dims": [], "unlimited": False, "vals": {}, "path": path, "kind": kind}
                (ds if kind == "DATASET" else attrs)[path] = cur
            continue
        if cur is not None and line.startswith("DATASPACE"):
            m = re.search(r"SIMPLE\s*\{\s*\(([^)]*)\)\s*/\s*\(([^)]*)\)", line)
            if m:
                cur["dims"] = [int(x) for x in m.group(1).split(",")]
                cur["unlimited"] = "UNLIMITED" in m.group(2).split(",")[0]
            continue
        if line == "DATA {":
            indata = True
            continue
        if indata:
            if line == "}":
                indata = False
                continue
            m = re.match(r"^\(([0-9,\s]*)\):\s*(.*)$", line)
            if m and cur is not None:
                idx = m.group(1).split(",")[0].strip() or "0"
                vals = [tok(t) for t in m.group(2).split(",") if t.strip() != ""]
                cur["vals"].setdefault(int(idx), []).extend(vals)
            continue
        if line == "}":
            if stack:
                k, _ = stack.pop()
                if k in ("DATASET", "ATTRIBUTE"):
                    cur = None
                    # an attribute inside a dataset: restore the dataset as current
                    for kk, nn in reversed(stack):
                        if kk == "DATASET":
                            groups = [n for k2, n in stack if k2 == "GROUP" and n != "/"]
                            cur = ds.get("/" + "/".join(groups + [nn]))
                            break
            continue
        if line.endswith("{"):
            stack.append(("OTHER", line))
    return ds, attrs


def dedup(bins):
    """bins: list of tuples (one per bin, containing the values of all bin datasets) -> (kept indices, duplicate positions)"""
    keep, dups = [], []
    for i, b in enumerate(bins):
        if i > 0 and b == bins[i - 1]:
            dups.append(i)
        else:
            keep.append(i)
    return keep, dups


def main():
    args = sys.argv[1:]
    tol, h5d = 1e-10, os.environ.get("H5DUMP", "h5dump")
    if "--tol" in args:
        i = args.index("--tol"); tol = float(args[i + 1]); del args[i:i + 2]
    if "--h5dump" in args:
        i = args.index("--h5dump"); h5d = args[i + 1]; del args[i:i + 2]
    if len(args) != 2:
        print(__doc__)
        return 2
    A, B = args
    da, aa = parse(dump(A, h5d))
    db, ab = parse(dump(B, h5d))
    fails = []
    if sorted(da) != sorted(db):
        fails.append("dataset lists differ: only in A %s, only in B %s" % (sorted(set(da) - set(db)), sorted(set(db) - set(da))))
    common = sorted(set(da) & set(db))
    binned = [p for p in common if da[p]["unlimited"] and db[p]["unlimited"]]
    static = [p for p in common if p not in binned]

    def nb(d, p):
        return len(d[p]["vals"])
    nA = min([nb(da, p) for p in binned] or [0]); nA_max = max([nb(da, p) for p in binned] or [0])
    nB = min([nb(db, p) for p in binned] or [0]); nB_max = max([nb(db, p) for p in binned] or [0])
    if nA != nA_max:
        fails.append("A: bin datasets disagree in length (%d..%d)" % (nA, nA_max))
    if nB != nB_max:
        fails.append("B: bin datasets disagree in length (%d..%d)" % (nB, nB_max))
    dupsA = dedup([tuple(tuple(da[p]["vals"][i]) for p in binned) for i in range(nA)])[1] if binned else []
    keepB, dupsB = dedup([tuple(tuple(db[p]["vals"][i]) for p in binned) for i in range(nB)]) if binned else ([], [])
    print("A: %s  bins=%d  identical consecutive bins at %s" % (A, nA, dupsA or "none"))
    print("B: %s  bins=%d  identical consecutive bins at %s" % (B, nB, dupsB or "none"))
    if binned and nA != nB:
        fails.append("raw number of bins differs: A %d, B %d" % (nA, nB))
    if dupsA:
        fails.append("reference A contains duplicated bins at %s" % dupsA)
    if len(dupsB) > 1:
        fails.append("B contains %d duplicated bins (at most one expected per kill)" % len(dupsB))
    keepA = list(range(nA))
    ncmp = min(len(keepA), len(keepB))
    if binned:
        print("B after dropping %d duplicate(s): %d unique bins; comparing the first %d bins of A and B" % (len(dupsB), len(keepB), ncmp))

    maxrel, maxscaled, rows = 0.0, 0.0, []

    def cmpvals(a, b):
        """-> (maxrel, maxscaled_abs_diff, maxabs_a, problem)"""
        if len(a) != len(b):
            return 0, 0, 0, "length %d vs %d" % (len(a), len(b))
        mr = md = ma = 0.0
        for x, y in zip(a, b):
            if isinstance(x, float) and isinstance(y, float):
                d = abs(x - y)
                s = max(abs(x), abs(y))
                mr = max(mr, d / s if s > 0 else 0.0)
                md, ma = max(md, d), max(ma, abs(x))
            elif x != y:
                return mr, md, ma, "%r != %r" % (x, y)
        return mr, md, ma, None

    for p in binned:
        mr = md = ma = 0.0
        prob = None
        for k in range(ncmp):
            r, d, a, pr = cmpvals(da[p]["vals"][keepA[k]], db[p]["vals"][keepB[k]])
            mr, md, ma = max(mr, r), max(md, d), max(ma, a)
            prob = prob or ("bin %d: %s" % (k, pr) if pr else None)
        sc = md / ma if ma > 0 else 0.0
        rows.append((p, "bins", ncmp, mr, sc, prob))
    for p in static:
        va = [v for k in sorted(da[p]["vals"]) for v in da[p]["vals"][k]]
        vb = [v for k in sorted(db[p]["vals"]) for v in db[p]["vals"][k]]
        mr, md, ma, prob = cmpvals(va, vb)
        sc = md / ma if ma > 0 else 0.0
        if da[p]["dims"] != db[p]["dims"]:
            prob = prob or "dims %s vs %s" % (da[p]["dims"], db[p]["dims"])
        rows.append((p, "static", len(va), mr, sc, prob))
    print("\n| dataset | kind | n | max rel dev | max scaled dev | note |\n|---|---|---|---|---|---|")
    for p, kind, n, mr, sc, prob in rows:
        print("| %s | %s | %d | %.3e | %.3e | %s |" % (p, kind, n, mr, sc, prob or ""))
        maxrel, maxscaled = max(maxrel, mr), max(maxscaled, sc)
        if prob:
            fails.append("%s: %s" % (p, prob))
        elif sc > tol:
            fails.append("%s: scaled deviation %.3e > tol %.1e" % (p, sc, tol))
    adiff = [k for k in sorted(set(aa) | set(ab)) if aa.get(k, {}).get("vals") != ab.get(k, {}).get("vals")]
    if adiff:
        print("\nattributes that differ (informational, not part of the verdict): %s" % ", ".join(adiff[:12]) + (" ..." if len(adiff) > 12 else ""))
    print("")
    for f in fails:
        print("FAIL: " + f)
    ok = not fails
    print("DUPS=%d" % len(dupsB))
    print("NBINS_A=%d NBINS_B=%d" % (nA, nB))
    print("MAXREL=%.3e MAXSCALED=%.3e" % (maxrel, maxscaled))
    print("RESULT=%s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
