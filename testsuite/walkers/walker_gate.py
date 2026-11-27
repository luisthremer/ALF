#!/usr/bin/env python3
"""Gate of task T12 (several walkers per process): N walkers in one process run the same chains as N
separate runs with the same seeds, chain by chain.

For each scenario of testsuite/test_branch_parameters.json named below (2 bins), with C = 4 chains:
  W   one noMPI process, N_walkers = 4, full seeds file                -> confout_0..3
  S_c one noMPI process per chain c, seeds file = line C - c only      -> confout_0 == W confout_c
  M   (with --mpi-exe) 4 MPI ranks, N_walkers = 1                      -> confout_c == W confout_c
  M2  (with --mpi-exe) 2 MPI ranks, N_walkers = 2                      -> confout_c == W confout_c
  R   restart: 1 bin, confout -> confin, 1 more bin, with N_walkers = 4 and per chain -> same confout_c
confout_c holds the field configuration and the random number generator state, so "byte-identical" means
the chain is the same. Acceptance lines in info are compared too; observable bins (W against M) only
within roundoff (sums over walkers and over ranks are added in a different order), reported, not gated.
Control: W confout_1 must differ from S_0 (the comparison sees which chain is which).

usage: walker_gate.py <repo> <noMPI ALF.out> <workdir> [--mpi-exe ALF.out --mpiexec "mpirun --oversubscribe"]
"""
import argparse, filecmp, os, re, shutil, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../Lowrank.tests"))
from compare_schemes import make_text, prepare, execute, restart_files, set_parameter  # noqa: E402
import json  # noqa: E402

SCENARIOS = ["Hubbard_N_Leg_Ladder", "Hubbard_PAM_1", "Kondo", "LRC", "tV", "Z2_Matter", "Spin_Peierls",
             "Hubbard_Global_moves", "Hubbard_Bilayer"]
C = 4


def with_walkers(text, n):
    return text.replace("&VAR_QMC", "&VAR_QMC\nN_walkers = %d" % n, 1)


def acceptance(path):
    return [l for l in open(path, errors="replace") if re.search(r"(?i)accept", l)]


def numbers(path):
    out = []
    for l in open(path, errors="replace"):
        for t in l.split():
            try:
                out.append(float(t.replace("D", "E")))
            except ValueError:
                pass
    return out


def max_rel(a_dir, b_dir):
    worst = 0.0
    for f in os.listdir(a_dir):
        if not re.search(r"_(scal|eq|tau)$", f):
            continue
        a, b = numbers(os.path.join(a_dir, f)), numbers(os.path.join(b_dir, f))
        if len(a) != len(b):
            return float("inf")
        for x, y in zip(a, b):
            worst = max(worst, abs(x - y) / max(abs(x), abs(y), 1e-300) if x != y else 0.0)
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("repo"); ap.add_argument("exe"); ap.add_argument("work")
    ap.add_argument("--mpi-exe"); ap.add_argument("--mpiexec", default="mpirun --oversubscribe")
    ap.add_argument("--scenarios", default=",".join(SCENARIOS))
    a = ap.parse_args()
    start = os.path.join(a.repo, "Scripts_and_Parameters_files/Start")
    template = open(os.path.join(start, "parameters")).read()
    seed_lines = open(os.path.join(start, "seeds")).read().split()[:C]
    seeds = os.path.join(a.work, "seeds_%d" % C)
    os.makedirs(a.work, exist_ok=True)
    open(seeds, "w").write("\n".join(seed_lines) + "\n")
    d = json.load(open(os.path.join(a.repo, "testsuite/test_branch_parameters.json")))
    bad = 0
    for name in a.scenarios.split(","):
        model, pars = d[name]
        text, missing = make_text(template, pars, 2, None)
        if name == "Hubbard_Bilayer":  # also the low-rank scheme (fixed N_delay: the timed choice is not reproducible)
            text = text.replace("&VAR_QMC", '&VAR_QMC\nSequential_scheme = "submatrix"\nN_delay = 3', 1)
        if missing:
            print("SKIP %-22s missing %s" % (name, missing)); continue
        base = os.path.join(a.work, name)
        fails = []

        def run(lab, txt, sd, exe, ranks=1, launcher=None, restart_from=None):
            dd = os.path.join(base, lab)
            prepare(dd, [txt], sd, 1, "serial")
            if restart_from:
                for f in os.listdir(restart_from):
                    if f.startswith("confout"):
                        shutil.copy(os.path.join(restart_from, f), dd)
                restart_files(dd)
            rc = execute(exe, dd, launcher, ranks)
            if rc:
                fails.append("%s rc=%d" % (lab, rc))
            return dd

        def same(x, y, what):
            if not (os.path.exists(x) and os.path.exists(y) and filecmp.cmp(x, y, shallow=False)):
                fails.append(what)

        W = run("W", with_walkers(text, C), seeds, a.exe)
        for c in range(C):
            sd = os.path.join(base, "seed_c%d" % c)
            open(sd, "w").write(seed_lines[C - 1 - c] + "\n")  # chain c reads line C - c
            S = run("S%d" % c, text, sd, a.exe)
            same(os.path.join(W, "confout_%d" % c), os.path.join(S, "confout_0"), "S%d" % c)
        if filecmp.cmp(os.path.join(W, "confout_1"), os.path.join(base, "S0", "confout_0"), shallow=False):
            fails.append("control: W1 == S0")
        # restart (1 bin, confout -> confin, 1 bin) against the restart of the separate run of each chain
        one = set_parameter(text, "NBin", 1)[0]
        R1 = run("R1", with_walkers(one, C), seeds, a.exe)
        R2 = run("R2", with_walkers(one, C), seeds, a.exe, restart_from=R1)
        for c in range(C):
            sd = os.path.join(base, "seed_c%d" % c)
            T1 = run("T1_%d" % c, one, sd, a.exe)
            T2 = run("T2_%d" % c, one, sd, a.exe, restart_from=T1)
            same(os.path.join(R2, "confout_%d" % c), os.path.join(T2, "confout_0"), "restart c%d" % c)
        extra = ""
        if a.mpi_exe:
            M = run("M", text, seeds, a.mpi_exe, C, a.mpiexec)
            M2 = run("M2", with_walkers(text, 2), seeds, a.mpi_exe, 2, a.mpiexec)
            for c in range(C):
                same(os.path.join(W, "confout_%d" % c), os.path.join(M, "confout_%d" % c), "M c%d" % c)
                same(os.path.join(W, "confout_%d" % c), os.path.join(M2, "confout_%d" % c), "M2 c%d" % c)
            if acceptance(os.path.join(W, "info")) != acceptance(os.path.join(M, "info")):
                fails.append("acceptance W/M")
            extra = "obs max rel W/M %.1e" % max_rel(W, M)
        print("%-4s %-22s %s %s" % ("OK" if not fails else "FAIL", name, " ".join(fails), extra))
        bad += bool(fails)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
