#!/usr/bin/env python3
"""Compare two Metropolis-decision traces (trace_updates.txt) of ALF runs.

Usage: trace_diff.py TRACE_A TRACE_B
       trace_diff.py --selftest

Each trace line (written by Prog/Trace_updates_mod.F90, ALF built with -DALF_TRACE_UPDATES) is
    nt n_op Re(R) Im(R) Weight ranf toggle
For runs with identical seeds and parameters the RNG is consumed in the same order, so line i of both
traces is the same proposal (same nt, n_op, same ranf) as long as all earlier decisions agree.

The traces are walked in step. Reported: the number of proposals compared and the maximum relative
ratio difference |R_a-R_b|/max(|R_a|,|R_b|) over the common prefix up to and including the first
differing toggle (after a differing toggle the configurations differ and the ratios are not comparable).
At the first differing toggle the decision is classified by the margin |Weight-ranf|/Weight of both runs:
both <= 1e-10 -> ROUNDOFF (the two runs sit on opposite sides of the Metropolis threshold within
floating-point noise; the trajectories legitimately diverge afterwards), otherwise BUG.
Structural problems before any differing toggle (different nt/n_op or ranf at line i, different number of
lines, empty traces) are BUG as well.

The last stdout line is "RESULT: IDENTICAL ...", "RESULT: ROUNDOFF k margin ..." or "RESULT: BUG ...".
Exit code: 0 for IDENTICAL and ROUNDOFF, 1 for BUG, 2 for unusable input (missing or empty trace).
"""
import math
import os
import sys
import tempfile

MARGIN_TOL = 1e-10


def read_trace(path):
    rows = []
    with open(path) as f:
        for ln, line in enumerate(f, 1):
            p = line.split()
            if not p:
                continue
            if len(p) != 7:
                sys.exit("trace_diff.py: %s:%d: expected 7 fields, got %d" % (path, ln, len(p)))
            rows.append((int(p[0]), int(p[1]), complex(float(p[2]), float(p[3])),
                         float(p[4]), float(p[5]), int(p[6])))
    return rows


def margin(w, r):
    return abs(w - r) / w if w > 0.0 else math.inf


def rel_diff(a, b):
    d = max(abs(a), abs(b))
    return 0.0 if d == 0.0 else abs(a - b) / d


def compare(ta, tb):
    """Returns (status, message, exit_code); status in IDENTICAL, ROUNDOFF, BUG, ERROR."""
    if not ta or not tb:
        return "ERROR", "ERROR empty trace (A: %d lines, B: %d lines): was the binary built with -DALF_TRACE_UPDATES?" \
            % (len(ta), len(tb)), 2
    n = min(len(ta), len(tb))
    maxrel = 0.0
    for i in range(n):
        a, b = ta[i], tb[i]
        k = i + 1
        if a[0] != b[0] or a[1] != b[1] or a[4] != b[4]:
            print("proposals compared: %d, max relative ratio difference: %.3e" % (i, maxrel))
            return "BUG", ("BUG proposal %d differs in (nt, n_op, ranf): A=(%d,%d,%.17g) B=(%d,%d,%.17g), "
                           "RNG order or proposal sequence out of step"
                           % (k, a[0], a[1], a[4], b[0], b[1], b[4])), 1
        maxrel = max(maxrel, rel_diff(a[2], b[2]))
        if a[5] != b[5]:
            ma, mb = margin(a[3], a[4]), margin(b[3], b[4])
            print("proposals compared: %d, max relative ratio difference: %.3e" % (k, maxrel))
            print("first differing toggle at proposal %d (nt=%d n_op=%d): A toggle=%d Weight=%.17g ranf=%.17g, "
                  "B toggle=%d Weight=%.17g ranf=%.17g" % (k, a[0], a[1], a[5], a[3], a[4], b[5], b[3], b[4]))
            print("margin |Weight-ranf|/Weight: A %.3e, B %.3e" % (ma, mb))
            if max(ma, mb) <= MARGIN_TOL:
                return "ROUNDOFF", "ROUNDOFF %d %.3e" % (k, max(ma, mb)), 0
            return "BUG", "BUG at proposal %d (margin %.3e exceeds %.0e)" % (k, max(ma, mb), MARGIN_TOL), 1
    print("proposals compared: %d, max relative ratio difference: %.3e" % (n, maxrel))
    if len(ta) != len(tb):
        return "BUG", ("BUG traces have different lengths (%d vs %d) but identical decisions over the common prefix"
                      % (len(ta), len(tb))), 1
    return "IDENTICAL", "IDENTICAL decisions over all %d proposals" % n, 0


def selftest():
    def row(nt, n_op, r, w, ran, tog):
        return "%d %d %.16e %.16e %.16e %.16e %d\n" % (nt, n_op, r.real, r.imag, w, ran, tog)

    base = [row(1, 1, 0.5 + 0j, 0.5, 0.3, 1), row(1, 2, 0.2 + 0j, 0.2, 0.7, 0), row(1, 3, 1.5 + 0j, 1.5, 0.9, 1)]
    cases = []
    # large margin: B rejects a proposal that A accepts with Weight 0.5 vs ranf 0.3
    big = list(base)
    big[0] = row(1, 1, 0.5 + 0j, 0.5, 0.3, 0)
    cases.append(("large margin", big, "BUG", 1))
    # margin 1e-13: A accepts, B rejects, Weight within 1e-13 of ranf
    ra = list(base)
    ra[0] = row(1, 1, 0.5 + 0j, 0.3 * (1 + 1e-13), 0.3, 1)
    rb = list(base)
    rb[0] = row(1, 1, 0.5 + 0j, 0.3 * (1 - 1e-13), 0.3, 0)
    cases.append(("margin 1e-13", (ra, rb), "ROUNDOFF", 0))
    cases.append(("identical", base, "IDENTICAL", 0))
    short = base[:2]
    cases.append(("different length", (base, short), "BUG", 1))
    ranf_bad = list(base)
    ranf_bad[1] = row(1, 2, 0.2 + 0j, 0.2, 0.71, 0)
    cases.append(("ranf out of step", ranf_bad, "BUG", 1))
    ok = True
    with tempfile.TemporaryDirectory() as d:
        for name, spec, want, want_code in cases:
            a, b = spec if isinstance(spec, tuple) else (base, spec)
            pa, pb = os.path.join(d, "a"), os.path.join(d, "b")
            open(pa, "w").writelines(a)
            open(pb, "w").writelines(b)
            status, msg, code = compare(read_trace(pa), read_trace(pb))
            good = status == want and code == want_code
            ok &= good
            print("selftest %-18s -> %-9s exit %d  %s" % (name, status, code, "ok" if good else "FAIL"))
    return 0 if ok else 1


def main(argv):
    if len(argv) == 2 and argv[1] == "--selftest":
        return selftest()
    if len(argv) != 3:
        sys.exit(__doc__)
    for p in argv[1:]:
        if not os.path.isfile(p):
            print("trace_diff.py: missing trace file " + p)
            print("RESULT: ERROR missing trace file")
            return 2
    status, msg, code = compare(read_trace(argv[1]), read_trace(argv[2]))
    print("RESULT: " + msg)
    return code


if __name__ == "__main__":
    sys.exit(main(sys.argv))
