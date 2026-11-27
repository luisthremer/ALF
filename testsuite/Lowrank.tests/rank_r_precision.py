#!/usr/bin/env python3
"""Precision of the rank-r accept of the delayed and submatrix updates (KNOWN_ISSUES.md, G11).

Spin_Peierls 12x12 has non-diagonal rank-2 vertices (type-4 J terms). The delayed and submatrix accepts
applied them in the Woodbury form (G(:,P) Delta Mat^-1) V^T: Precision Green mean 2.5-2.8e-10 against
9.0e-12 of the standard update (Upgrade2's grouping (G(:,P) X) (V Tri)^T). Runs standard, delayed and
submatrix (N_delay = 4) on the same input and seeds. Passes if every run ends with return code 0 and the
Precision Green mean of delayed and submatrix is at most RATIO times that of standard.

Usage: rank_r_precision.py --exe <ALF>/Prog/ALF.out --start <ALF>/Scripts_and_Parameters_files/Start
                           --work DIR [--mpiexec "mpirun ..."]
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_schemes import SCHEMES, execute, make_text, prepare, with_scheme  # noqa: E402

BASE = {"ham_name": "Spin_Peierls", "Model": "Spin_Peierls", "L1": 12, "L2": 12, "Beta": 2.0, "ham_U": 1.5,
        "ham_Lambda": 0.1, "ham_omega0": 0.5, "N_SUN": 2, "N_FL": 1, "NSweep": 4, "NBin": 1, "Ltau": 0}
N_DELAY = 4
RATIO = 2.0


def precision_mean(d):
    info = os.path.join(d, "info")
    m = re.search(r"Precision Green\s+Mean, Max :\s*(\S+)\s+(\S+)", open(info).read() if os.path.exists(info) else "")
    return float(m.group(1).replace("D", "E")) if m else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--start", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--mpiexec", default=None)
    args = ap.parse_args()
    text, missing = make_text(open(os.path.join(args.start, "parameters")).read(), BASE)
    assert not missing, missing
    mean = {}
    bad = 0
    for scheme in SCHEMES:
        d = os.path.join(args.work, scheme)
        prepare(d, [with_scheme(text, scheme, 0 if scheme == "standard" else N_DELAY)],
                os.path.join(args.start, "seeds"), 1, "serial")
        rc = execute(args.exe, d, args.mpiexec, 1)
        mean[scheme] = precision_mean(d)
        if rc != 0:
            print("FAIL %-10s rc=%s" % (scheme, rc))
            bad += 1
    for scheme in SCHEMES:
        ratio = mean[scheme] / mean["standard"] if mean["standard"] > 0 else float("inf")
        ok = ratio <= RATIO
        print("%-4s %-10s Precision Green mean = %.3e  ratio to standard = %.2f (<= %.1f)"
              % ("OK" if ok else "FAIL", scheme, mean[scheme], ratio, RATIO))
        bad += not ok
    print("PASS: rank-r precision" if not bad else "FAILED: %d" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
