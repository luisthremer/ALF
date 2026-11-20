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
                       [--work DIR] [--only NAME ...] [--max-bins N] [--tol T] [--n-delay 0|-1]

--n-delay -1 runs the low-rank schemes with the timed automatic N_delay (the candidates switch between
slices and the choice follows the machine timing; the chain must still agree up to roundoff). The
"Auto N_delay:" decision line of each run is printed when timing took place (>= 3 sweeps, at least two candidates).
"""
import argparse
import json
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from baseline_guard import make_parameters, run  # noqa: E402

SCHEMES = ("standard", "delayed", "submatrix")

# Scenarios beyond test_branch_parameters.json: the low-rank schemes together with Global_tau
# moves, sequential MALA, global-tau MALA, continuous fields and the projector code
_HUB = {"Model": "Hubbard", "Lattice_type": "N_leg_ladder", "L1": 6, "L2": 2, "Beta": 4.0,
        "Ham_U": 4.0, "Ham_T": 1.0, "Nsweep": 20, "NBin": 5, "Ltau": 1, "Mz": True}
EXTRA = {
    # sequential vertices 1..12 (all sites) followed by the Global_tau moves on every slice
    "Hubbard_Global_tau": dict(_HUB, Global_tau_moves=True, N_Global_tau=2,
                               Nt_sequential_start=1, Nt_sequential_end=12),
    "Hubbard_Global_tau_MALA_seq": dict(_HUB, Continuous=True, Global_tau_MALA_moves=True,
                                        Nt_sequential_start=1, Nt_sequential_end=12),
    "Hubbard_Continuous": dict(_HUB, Continuous=True),
    "Hubbard_Sequential_MALA": dict(_HUB, Continuous=True, Sequential_MALA=True),
    "Hubbard_Global_tau_MALA": dict(_HUB, Continuous=True, Global_tau_MALA_moves=True),
    "Hubbard_Projector": dict(_HUB, Projector=True, Theta=4.0),
}
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
    ap.add_argument("--n-delay", type=int, default=0, choices=(0, -1),
                    help="N_delay of the low-rank runs: 0 static rule, -1 timed automatic")
    args = ap.parse_args()

    with open(os.path.join(args.start, "parameters")) as f:
        template = f.read()
    seeds = os.path.join(args.start, "seeds")
    with open(args.scenarios) as f:
        scenarios = json.load(f)

    for name, pars in EXTRA.items():
        scenarios[name] = ("Hubbard", pars)

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
            print("SKIP  %s (keys not in template: %s)" % (name, ", ".join(missing)))
            continue
        dirs = {}
        rcs = {}
        for scheme in SCHEMES:
            text = base + '\n&VAR_LOWRANK\nSequential_scheme = "%s"\nN_delay = %d\n/\n' % (scheme, args.n_delay)
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
            decision = ""
            if args.n_delay == -1:
                with open(os.path.join(dirs[scheme], "info")) as f:
                    lines = [l.strip() for l in f if l.startswith("Auto N_delay")]
                decision = lines[0] if lines else "no decision line (fewer than 3 sweeps, or a single candidate)"
            d, bad = max_rel_diff(dirs["standard"], dirs[scheme])
            if d is None:
                print("FAIL  %-26s %-9s %s has a different number of entries" % (name, scheme, bad))
                failed += 1
            else:
                ok = d <= args.tol
                failed += 0 if ok else 1
                print("%s  %-26s %-9s max rel diff %.2e" % ("PASS" if ok else "FAIL", name, scheme, d))
                if decision:
                    print("      %s" % decision)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
