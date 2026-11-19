#!/usr/bin/env python3
"""Compare the low-rank sequential schemes with the standard update on fixed seeds.

Written with AI assistance on master_ai; not part of the ALF baseline.

One ALF executable built with the LOWRANK switch runs every serial scenario of
testsuite/test_branch_parameters.json three times, with Sequential_scheme =
"standard", "delayed" and "submatrix" (N_delay = 0, automatic) and the same
seeds. The low-rank schemes compute the same determinant ratios and Green
functions in a different summation order and draw the same random numbers, so
the Markov chains coincide and the bins agree up to roundoff. Every number in
every bin file (and confout_0) is compared; the largest relative difference
max|a - b| / max(|a|, |b|, 1) per scheme must stay below --tol.

A differing number of tokens (for example an accept/reject decision that
flipped at a near-tie) is reported as a failure, since the chains then differ.

Usage:
    compare_schemes.py --exe <ALF>/Prog/ALF.out
                       --start <ALF>/Scripts_and_Parameters_files/Start
                       [--work DIR] [--only NAME ...] [--max-bins N] [--tol T]
"""
import argparse
import json
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from baseline_guard import make_parameters, run  # noqa: E402

SCHEMES = ("standard", "delayed", "submatrix")
IGNORED = {"info", "info_git", "parameters", "seeds", "stdout"}
NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?")


def numbers(path):
    with open(path) as f:
        text = f.read()
    return [float(t.replace("D", "E").replace("d", "e")) for t in NUMBER.findall(text)]


def max_rel_diff(dir_a, dir_b):
    """Largest relative difference over all comparable files, or None if the shapes differ."""
    worst = 0.0
    for name in sorted(os.listdir(dir_a)):
        if name in IGNORED or name.endswith("_info"):
            continue
        a = numbers(os.path.join(dir_a, name))
        b = numbers(os.path.join(dir_b, name))
        if len(a) != len(b):
            return None, name
        for x, y in zip(a, b):
            d = abs(x - y) / max(abs(x), abs(y), 1.0)
            if d > worst:
                worst = d
    return worst, None


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, help="ALF.out built with the LOWRANK switch")
    ap.add_argument("--start", required=True, help="directory with template 'parameters' and 'seeds'")
    ap.add_argument("--scenarios", default=os.path.join(here, "..", "test_branch_parameters.json"))
    ap.add_argument("--work", default="lowrank_runs")
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--max-bins", type=int, default=None)
    ap.add_argument("--tol", type=float, default=1e-8)
    args = ap.parse_args()

    with open(os.path.join(args.start, "parameters")) as f:
        template = f.read()
    seeds = os.path.join(args.start, "seeds")
    with open(args.scenarios) as f:
        scenarios = json.load(f)

    failed = 0
    for name, (_, pars) in scenarios.items():
        if args.only and name not in args.only:
            continue
        if isinstance(pars, list):
            continue
        if pars.get("Langevin") or pars.get("HMC") or pars.get("Sequential") is False:
            continue  # no sequential updates: nothing to compare
        base, missing = make_parameters(template, pars, args.max_bins)
        if missing:
            continue
        dirs = {}
        rcs = {}
        for scheme in SCHEMES:
            text = base + '\n&VAR_LOWRANK\nSequential_scheme = "%s"\nN_delay = 0\n/\n' % scheme
            dirs[scheme] = os.path.join(args.work, name, scheme)
            rcs[scheme] = run(os.path.abspath(args.exe), dirs[scheme], text, seeds)
        if rcs["standard"] != 0:
            print("SKIP  %s (standard run stops with exit code %d)" % (name, rcs["standard"]))
            continue
        for scheme in SCHEMES[1:]:
            if rcs[scheme] != 0:
                print("FAIL  %-26s %-9s exit code %d" % (name, scheme, rcs[scheme]))
                failed += 1
                continue
            d, bad = max_rel_diff(dirs["standard"], dirs[scheme])
            if d is None:
                print("FAIL  %-26s %-9s %s has a different number of entries" % (name, scheme, bad))
                failed += 1
            else:
                ok = d <= args.tol
                failed += 0 if ok else 1
                print("%s  %-26s %-9s max rel diff %.2e" % ("PASS" if ok else "FAIL", name, scheme, d))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
