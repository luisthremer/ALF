#!/usr/bin/env python3
"""Throughput before (shared-core pinning, production settings) and after tuning, beta = 4 and beta = L.
usage: compare_tuned.py OLD_BETA4.dat OLD_BETAL.dat NEW_CSV [NEW_CSV ...] -- NEW_BETAL_CSV [...] OUT.png
OLD_*: antrag_best.dat of the reports (columns L thr_cpu thr_gpu thr_gpub); NEW: analyse.py csv of the tuned runs (disjoint cores)."""
import csv, collections, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    import scienceplots  # noqa
    plt.style.use(["science", "no-latex"])
except Exception:
    pass
a = sys.argv[1:]
old4, oldL = a[0], a[1]
i = a.index("--")
new4, newL, out = a[2:i], a[i + 1:-1], a[-1]


def load_old(f):
    return {int(r["L"]): r for r in csv.DictReader(open(f), delimiter=" ")}


def load_new(files):
    best = collections.defaultdict(float)
    for f in files:
        for r in csv.DictReader(open(f)):
            if r["rc"] == "0" and r.get("throughput_node_h"):
                k = (r["side"], int(r["L"])); best[k] = max(best[k], float(r["throughput_node_h"]))
    return best


fig, ax = plt.subplots(2, 2, figsize=(11, 7.5), sharex="col", gridspec_kw=dict(height_ratios=[3, 1.6]))
for col, (title, oldf, newf) in enumerate([(r"$\beta = 4$", old4, new4), (r"$\beta = L$", oldL, newL)]):
    old, new = load_old(oldf), load_new(newf)
    Ls = sorted(old)
    a0, a1 = ax[0][col], ax[1][col]
    cpu = [(L, float(old[L]["thr_cpu"])) for L in Ls if float(old[L]["thr_cpu"]) == float(old[L]["thr_cpu"])]
    a0.plot(*zip(*cpu), "o-", color="C0", label="CPU node, ALF CPU")
    for side, name, c, m in [("gpu", "ALF GPU", "C2", "s"), ("gpub", "ALF GPU Intel-opt.", "C1", "^")]:
        o = [(L, float(old[L]["thr_" + side])) for L in Ls]
        a0.plot(*zip(*o), m + "--", color=c, mfc="none", alpha=.8, label=f"{name}, before tuning")
        n = [(L, new[(side, L)]) for L in Ls if (side, L) in new]
        if n:
            a0.plot(*zip(*n), m + "-", color=c, label=f"{name}, tuned")
    rr = [(L, new[("gpub", L)] / new[("gpu", L)]) for L in Ls if ("gpub", L) in new and ("gpu", L) in new]
    if rr:
        a1.plot(*zip(*rr), "^-", color="C1", label="tuned")
    ro = [(L, float(old[L]["thr_gpub"]) / float(old[L]["thr_gpu"])) for L in Ls]
    a1.plot(*zip(*ro), "^--", color="C1", mfc="none", alpha=.7, label="before tuning")
    a0.set_yscale("log"); a0.set_title(title); a0.grid(alpha=.3)
    a1.axhline(1, color="gray", lw=.8); a1.set_xlabel(r"system size $L$ ($N=L^2$)"); a1.grid(alpha=.3)
    if col == 0:
        a0.set_ylabel("work units per node-hour"); a1.set_ylabel("Intel-opt. / ALF GPU")
        a0.legend(fontsize=7); a1.legend(fontsize=7)
    if not any((s, L) in new for s in ("gpu", "gpub") for L in Ls):
        a0.text(.5, .1, "no tuned runs yet", transform=a0.transAxes, ha="center", color="red")
    a1.set_ylim(0.8, None)
fig.tight_layout(); fig.savefig(out, dpi=140)
