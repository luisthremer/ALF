#!/usr/bin/env python3
"""Log-log plot of the node throughput against L (N = L^2) with the asymptotic large-matrix limit: every curve ~ A L^-p (p = 6 at beta = 4,
7 at beta = L: cost per sweep ~ N^3 Nt). CPU: amplitude from all measured points; GPU codes: amplitude from the tail L >= TAIL.
The ratio of the amplitudes is the limiting speedup. usage: loglog_limit.py best4.dat best_L.dat OUT.png"""
import csv, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
try:
    import scienceplots  # noqa
    plt.style.use(["science", "no-latex"])
except Exception:
    pass
TAIL = 44
def load(f): return {int(r["L"]): r for r in csv.DictReader(open(f), delimiter=" ")}
fig, ax = plt.subplots(1, 2, figsize=(11.5, 5))
Lx = np.logspace(np.log10(16), np.log10(160), 100)
for a, f, title, p in [(ax[0], sys.argv[1], r"$\beta=4$", 6), (ax[1], sys.argv[2], r"$\beta=L$", 7)]:
    d = load(f); Ls = np.array(sorted(d), float)
    amp = {}; slope = {}
    for key, nm, col, m in [("thr_cpu", "CPU node", "C0", "o"), ("thr_gpu", "ALF GPU", "C2", "s"), ("thr_gpub", "ALF GPU Intel-opt.", "C1", "^")]:
        pts = [(L, float(d[int(L)][key])) for L in Ls if float(d[int(L)][key]) == float(d[int(L)][key]) and (key == "thr_cpu" or L >= TAIL)]
        allp = [(L, float(d[int(L)][key])) for L in Ls if float(d[int(L)][key]) == float(d[int(L)][key])]
        l, y = np.array(pts).T
        amp[key] = np.exp(np.mean(np.log(y) + p * np.log(l)))
        slope[key] = -np.polyfit(np.log(l), np.log(y), 1)[0]
        a.loglog(*np.array(allp).T, m, color=col, label=nm if key == "thr_cpu" else nm + " (tuned)")
        a.loglog(Lx, amp[key] * Lx ** -p, "--", color=col, alpha=.8, lw=1)
    lim = amp["thr_gpub"] / amp["thr_cpu"]; limb = amp["thr_gpu"] / amp["thr_cpu"]
    print(f"{title}: p={p}; limiting speedup Intel-opt. {lim:.1f}x, ALF GPU {limb:.1f}x; tail slopes (free fit, L>={TAIL}): GPU {slope['thr_gpu']:.2f}, Intel {slope['thr_gpub']:.2f}, CPU (all) {slope['thr_cpu']:.2f}")
    a.text(.97, .95, f"limit (dashed, $\\propto L^{{-{p}}}$):\nIntel-opt./CPU = {lim:.0f}$\\times$\nALF GPU/CPU = {limb:.0f}$\\times$", transform=a.transAxes, ha="right", va="top", fontsize=8)
    a.set_xlabel(r"system size $L$ ($N=L^2$)"); a.set_title(title); a.grid(alpha=.3, which="both")
    a.xaxis.set_minor_formatter(NullFormatter()); a.set_xticks([16, 24, 32, 48, 64, 96, 128, 160]); a.set_xticklabels([16, 24, 32, 48, 64, 96, 128, 160])
    sec = a.secondary_xaxis("top", functions=(lambda x: x ** 2, lambda n: np.sqrt(n))); sec.set_xlabel("matrix size $N$"); sec.xaxis.set_minor_formatter(NullFormatter()); sec.set_xticks([256, 1024, 4096, 16384]); sec.set_xticklabels([256, 1024, 4096, 16384])
ax[0].set_ylabel("work units per node-hour"); ax[0].legend(fontsize=7, loc="lower left")
fig.tight_layout(); fig.savefig(sys.argv[3], dpi=140)
