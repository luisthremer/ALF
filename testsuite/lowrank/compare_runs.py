#!/usr/bin/env python3
"""Compare the outputs of two ALF runs (confout_0 and the observable bin files).

Usage: compare_runs.py [--info-only] [--rtol R] DIR_A DIR_B SCHEME_A SCHEME_B

Checks (exit 0 iff all pass, 1 otherwise, 2 on unusable input):
  * info: for every non-"standard" scheme the run's info file must contain the line
    "Sequential scheme: <NAME>" (NAME = scheme in upper case), so a run that silently fell back
    to the standard update fails. For "standard" nothing is checked.
  * confout_0 (RNG state and field configuration): identical field by field (all tokens equal).
  * every other output file except info*, parameters, seeds, *_info, stdout.log and trace_updates.txt
    (the trace is compared by trace_diff.py) (bin files *_scal, *_eq, *_tau):
    same set of files in both runs, same shape, and every number within relative tolerance
    |a-b| <= rtol*max(|a|,|b|) + 1e-12*max|value in that file| (default rtol 1e-8; the second term
    absorbs round-off noise on values that are zero up to ~1e-16, e.g. imaginary parts). The absolute
    floor overrides rtol for small entries: an entry below 1e-4 of the file maximum can differ by up
    to 100% at rtol 1e-8 without failing.
--info-only restricts the checks to the info line (used after a legitimate round-off divergence of
the trajectories, where confout and bins differ by construction).
"""
import math
import os
import re
import sys

NUM = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?")


def numbers(path):
    out = []
    with open(path) as f:
        for line in f:
            out.append([float(t.replace("D", "e").replace("d", "e")) for t in NUM.findall(line)])
    return out


def check_info(d, scheme):
    if scheme.lower() == "standard":
        return True, "standard (no info line to check)"
    want = "Sequential scheme: " + scheme.upper()
    path = os.path.join(d, "info")
    if not os.path.isfile(path):
        return False, "no info file in " + d
    with open(path) as f:
        if any(line.strip() == want for line in f):
            return True, "'" + want + "' found"
    return False, "line '%s' missing in %s/info (silent fallback to standard?)" % (want, d)


def compare_confout(da, db):
    pa, pb = os.path.join(da, "confout_0"), os.path.join(db, "confout_0")
    if not (os.path.isfile(pa) and os.path.isfile(pb)):
        return False, "confout_0 missing"
    with open(pa) as fa, open(pb) as fb:
        ta, tb = fa.read().split(), fb.read().split()
    if len(ta) != len(tb):
        return False, "confout_0 has different sizes (%d vs %d tokens)" % (len(ta), len(tb))
    bad = [i for i, (x, y) in enumerate(zip(ta, tb)) if x != y]
    if bad:
        return False, "confout_0 differs in %d of %d fields (first at field %d: %s vs %s)" \
            % (len(bad), len(ta), bad[0], ta[bad[0]], tb[bad[0]])
    return True, "confout_0 identical (%d fields)" % len(ta)


def compare_bins(da, db, rtol):
    # stdout.log (run log) and trace_updates.txt (compared by trace_diff.py) are not bin files
    skip = lambda n: n in ("parameters", "seeds", "confout_0", "stdout.log", "trace_updates.txt") \
        or n.startswith("info") or n.endswith("_info")
    fa = sorted(n for n in os.listdir(da) if not skip(n) and os.path.isfile(os.path.join(da, n)))
    fb = sorted(n for n in os.listdir(db) if not skip(n) and os.path.isfile(os.path.join(db, n)))
    if fa != fb:
        return False, "different output files: only A %s, only B %s" % (sorted(set(fa) - set(fb)),
                                                                          sorted(set(fb) - set(fa)))
    if not any(n.endswith("_scal") for n in fa):
        return False, "no *_scal bin files found"
    worst, worst_file, ok = 0.0, "", True
    msgs = []
    for n in fa:
        a, b = numbers(os.path.join(da, n)), numbers(os.path.join(db, n))
        if [len(r) for r in a] != [len(r) for r in b]:
            return False, "%s: different line/field structure" % n
        flat_a = [v for r in a for v in r]
        flat_b = [v for r in b for v in r]
        if not all(math.isfinite(v) for v in flat_a + flat_b):
            return False, "%s: non-finite value" % n
        scale = max([abs(v) for v in flat_a + flat_b] + [0.0])
        for x, y in zip(flat_a, flat_b):
            err = abs(x - y)
            lim = rtol * max(abs(x), abs(y)) + 1e-12 * scale
            rel = err / max(abs(x), abs(y), 1e-300) if err > 0 else 0.0
            if err > lim:
                ok = False
                msgs.append("%s: %.17g vs %.17g (rel %.2e)" % (n, x, y, rel))
                break
            if err > 0 and err <= lim and rel > worst:
                worst, worst_file = rel, n
    if not ok:
        return False, "bins differ beyond rtol=%g: %s" % (rtol, "; ".join(msgs[:5]))
    return True, "%d bin files within rtol=%g (max relative difference %.2e%s)" \
        % (len(fa), rtol, worst, " in " + worst_file if worst_file else "")


def main(argv):
    args = argv[1:]
    info_only, rtol = False, 1e-8
    while args and args[0].startswith("--"):
        o = args.pop(0)
        if o == "--info-only":
            info_only = True
        elif o == "--rtol":
            rtol = float(args.pop(0))
        else:
            sys.exit(__doc__)
    if len(args) != 4:
        sys.exit(__doc__)
    da, db, sa, sb = args
    for d in (da, db):
        if not os.path.isdir(d):
            print("compare_runs.py: no such run directory " + d)
            return 2
    ok = True
    for tag, d, s in (("A", da, sa), ("B", db, sb)):
        good, msg = check_info(d, s)
        ok &= good
        print("info %s: %s %s" % (tag, "ok  " if good else "FAIL", msg))
    if not info_only:
        for fn in (compare_confout, lambda x, y: compare_bins(x, y, rtol)):
            good, msg = fn(da, db)
            ok &= good
            print("%s %s" % ("ok  " if good else "FAIL", msg))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
