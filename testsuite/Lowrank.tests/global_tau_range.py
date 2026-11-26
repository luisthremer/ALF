#!/usr/bin/env python3
"""Global_tau_moves with a sequential range that does not cover every vertex (KNOWN_ISSUES.md, H5).

WRAPGRUP / WRAPGRDO must move G through the vertices outside Nt_sequential_start..Nt_sequential_end also
when no global tau move is done (N_Global_tau = 0). Before the fix G was not wrapped through them and
"Precision Green" reached ~1 with return code 0. Hubbard 4x4 (Mz, U = 1, 16 vertices per slice):

  partial: Nt_sequential_start = 3, Nt_sequential_end = 16, N_Global_tau = 0
  empty:   empty range (template: 0 .. -1), N_Global_tau = 0

each with Sequential_scheme = standard, delayed and submatrix. Passes if every run ends with return code 0
and Precision Green Max < 1e-8 (start = 1, end = 16 gives ~1e-14).

Usage: global_tau_range.py --exe <ALF>/Prog/ALF.out --start <ALF>/Scripts_and_Parameters_files/Start
                           --work DIR [--mpiexec "mpirun ..."]
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_schemes import SCHEMES, execute, make_text, prepare, with_scheme  # noqa: E402

CASES = {
    "partial": {"Nt_sequential_start": 3, "Nt_sequential_end": 16},
    "empty": {},
}
BASE = {"L1": 4, "L2": 4, "Beta": 2.0, "NSweep": 5, "NBin": 1, "Ltau": 0, "ham_U": 1.0,
        "Global_tau_moves": True, "N_Global_tau": 0}
TOL = 1e-8


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--start", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--mpiexec", default=None)
    args = ap.parse_args()
    template = open(os.path.join(args.start, "parameters")).read()
    bad = 0
    for case, pars in CASES.items():
        text, missing = make_text(template, dict(BASE, **pars))
        assert not missing, missing
        for scheme in SCHEMES:
            d = os.path.join(args.work, "%s-%s" % (case, scheme))
            prepare(d, [with_scheme(text, scheme, 0)], os.path.join(args.start, "seeds"), 1, "serial")
            rc = execute(args.exe, d, args.mpiexec, 1)
            info = os.path.join(d, "info")
            m = re.search(r"Precision Green\s+Mean, Max :\s*(\S+)\s+(\S+)",
                          open(info).read() if os.path.exists(info) else "")
            xmax = float(m.group(2).replace("D", "E")) if m else float("nan")
            ok = rc == 0 and xmax < TOL
            print("%-4s %-8s %-10s rc=%s Precision Green Max = %s" % ("OK" if ok else "FAIL", case, scheme, rc, xmax))
            bad += not ok
    print("PASS: global tau range" if not bad else "FAILED: %d run(s)" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
