#!/usr/bin/env python3
"""Job 3 (P2.8.7): long CPU runs against exact diagonalization.  Cases = the 4-site Hubbard chain of testsuite/test_vs_ed/test_specs.yaml
(finite T and projector, ED values from there) plus SU(2) form.  Per seed: dtau scan, fit E(dtau) = y0 + a dtau^2 (weights 1/err^2, as
analysis.py), y0 compared with ED.  Counts |y0-ED| < 1 sigma, 2 sigma, 3 sigma over the seeds (expected 68.3 / 95.4 / 99.7 %).
usage: ed_stats.py EXE WORK OUT NSEEDS JOBS"""
import os, re, sys, math, shutil, subprocess
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "docs", "ai", "phase2", "p25_runs")); import p25_run as p25
START = os.path.join(REPO, "Scripts_and_Parameters_files", "Start")
BASE = {"Lattice_type": '"N_leg_ladder"', "L1": "4", "L2": "1", "Checkerboard": ".F.", "Symm": ".T.", "ham_T": "1.d0", "ham_U": "4.d0",
        "ham_Tperp": "0.d0", "Ltau": "0", "NSweep": "1000", "NBin": "300"}
DT = [(0.05, 25), (0.1, 15), (0.14285714285714285, 10), (0.2, 5)]
ALL_DT = DT; DT = DT[:int(os.environ.get("ED_NDT", "4"))]      # ED_NDT=3: fit without dtau = 0.2 (dtau^4 bias check)
CASES = {  # name: (edits, ED energy)
 "finT-mz":   (dict(BASE, Beta="2.d0", Mz=".T."), -1.47261997),
 "finT-su2":  (dict(BASE, Beta="2.d0", Mz=".F."), -1.47261997),
 "T0-mz":     (dict(BASE, Beta="1.d0", Projector=".T.", Theta="5.d0", Mz=".T.", NBin="200"), -2.10274848),
}
def one(args):
    exe, d, ed, dtau, nw, seed = args
    if os.path.exists(d + "/.done"): return
    p25.prepare(d, START, dict(ed, Dtau=repr(dtau), Nwrap=str(nw)))
    open(d + "/seeds", "w").write("%d\n" % seed)
    p25.run(exe, d, {}, 0, 7200); open(d + "/.done", "w").write("1")
def ener(d):
    v = []
    for l in open(d + "/Ener_scal"):
        m = re.match(r"\s*\d+\s+\(\s*(\S+),\s*\S+\)\s+(\S+)", l)
        if m: v.append(float(m.group(1).replace("D", "E")) / float(m.group(2).replace("D", "E")))   # sign-normalised bin value
    v = v[1:]; m = sum(v) / len(v); return m, math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1) / len(v))      # n_skip = 1
def main():
    exe, W, out, ns, jobs = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
    tasks = [(exe, "%s/%s/s%d/dt%g" % (W, c, s, dt), ed, dt, nw, 100003 + 7919 * s) for c, (ed, _) in CASES.items() for s in range(ns) for dt, nw in ALL_DT]
    with ThreadPoolExecutor(jobs) as p: list(p.map(one, tasks))
    os.makedirs(out, exist_ok=True)
    rows = ["case\tseeds\tED\tmean_y0\tstd_of_y0\tmean_sigma\thit_1sigma\thit_2sigma\thit_3sigma\tpull_mean\tpull_std"]
    det = ["case\tseed\ty0\tsigma\tpull"]
    for c, (_, ref) in CASES.items():
        pulls, y0s, sg = [], [], []
        for s in range(ns):
            E = [ener("%s/%s/s%d/dt%g" % (W, c, s, dt)) for dt, _ in DT]
            x = [dt for dt, _ in DT]; y = [e[0] for e in E]; w = [1 / e[1] ** 2 for e in E]; z = [t * t for t in x]
            S, Sz, Szz = sum(w), sum(a * b for a, b in zip(w, z)), sum(a * b * b for a, b in zip(w, z))   # normal equations of y0 + a z
            Sy, Szy = sum(a * b for a, b in zip(w, y)), sum(a * b * c for a, b, c in zip(w, z, y)); det_ = S * Szz - Sz * Sz
            p0 = (Szz * Sy - Sz * Szy) / det_; sig = math.sqrt(Szz / det_); pulls.append((p0 - ref) / sig); y0s.append(p0); sg.append(sig)
            det.append("%s\t%d\t%.6f\t%.6f\t%.3f" % (c, s, p0, sig, pulls[-1]))
        pu = [abs(q) for q in pulls]; n = len(pu)
        mean = lambda l: sum(l) / len(l); std = lambda l: math.sqrt(sum((q - mean(l)) ** 2 for q in l) / (len(l) - 1))
        rows.append("%s\t%d\t%.8f\t%.6f\t%.6f\t%.6f\t%d/%d\t%d/%d\t%d/%d\t%.2f\t%.2f" % (c, n, ref, mean(y0s), std(y0s), mean(sg),
                    sum(q < 1 for q in pu), n, sum(q < 2 for q in pu), n, sum(q < 3 for q in pu), n, mean(pulls), std(pulls)))
    open(out + "/ed_stats%s.tsv" % os.environ.get("ED_SUF", ""), "w").write("\n".join(rows) + "\n"); open(out + "/ed_stats_seeds%s.tsv" % os.environ.get("ED_SUF", ""), "w").write("\n".join(det) + "\n"); print("\n".join(rows))
main()
