#!/usr/bin/env python3
"""CPU node throughput fitted with a power law A * L^-p (p fixed: cost per sweep ~ N^3 = L^6 at beta = 4, ~ L^7 at beta = L because
Nt = 10 L), extended to L = 56, with the tuned GPU curves and the speedups (measured: plain, extrapolated: in brackets).
usage: cpu_fit.py antrag_best.dat P OUT_PREFIX     -> OUT_PREFIX.pdf and OUT_PREFIX_table.tex (macro \\FitTable)"""
import csv, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    import scienceplots  # noqa
    plt.style.use(["science", "no-latex"])
except Exception:
    pass
d = {int(r["L"]): r for r in csv.DictReader(open(sys.argv[1]), delimiter=" ")}
p = int(sys.argv[2]); out = sys.argv[3]
Ls = sorted(d)
cl = [L for L in Ls if float(d[L]["thr_cpu"]) == float(d[L]["thr_cpu"])]
c = np.array([float(d[L]["thr_cpu"]) for L in cl]); l = np.array(cl, float)
A = np.exp(np.mean(np.log(c) + p * np.log(l)))
pf = -np.polyfit(np.log(l), np.log(c), 1)[0]
fig, a = plt.subplots(figsize=(6.4, 4.6))
Lx = np.linspace(16, 56, 100)
a.plot(l, c, "o", color="C0", label="CPU node (measured)")
a.plot(Lx, A * Lx ** -p, "--", color="C0", label=rf"CPU fit $\propto L^{{-{p}}}$ (free fit: $L^{{-{pf:.1f}}}$)")
a.plot(Ls, [float(d[L]["thr_gpu"]) for L in Ls], "s-", color="C2", label="ALF GPU")
a.plot(Ls, [float(d[L]["thr_gpub"]) for L in Ls], "^-", color="C1", label="ALF GPU Intel-opt.")
rows = []
for L in Ls:
    g = float(d[L]["thr_gpub"]); gb = float(d[L]["thr_gpu"])
    meas = L in cl
    cpu = float(d[L]["thr_cpu"]) if meas else A * L ** -p
    fac, facb = g / cpu, gb / cpu
    a.annotate((f"{fac:.1f}$\\times$" if meas else f"({fac:.0f}$\\times$)"), (L, g), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=7, color="C1")
    rows.append((L, cpu, meas, g, fac, facb))
a.set_yscale("log"); a.set_xlabel(r"system size $L$ ($N=L^2$)"); a.set_ylabel("work units per node-hour"); a.grid(alpha=.3); a.legend(fontsize=7)
fig.tight_layout(); fig.savefig(out + ".pdf"); fig.savefig(out + ".png", dpi=140)
t = ["\\begin{tabular}{rrrrrr}", "$L$ & CPU node & Intel-opt. GPU & speedup & ALF GPU speedup & CPU\\\\", "\\hline"]
for L, cpu, meas, g, fac, facb in rows:
    f1 = f"{fac:.1f}" if meas else f"({fac:.0f})"; f2 = f"{facb:.1f}" if meas else f"({facb:.0f})"
    t.append(f"{L} & {cpu:,.0f} & {g:,.0f} & {f1}$\\times$ & {f2}$\\times$ & {'measured' if meas else 'fit'}\\\\".replace(",", "\\,"))
t.append("\\end{tabular}")
open(out + "_table.tex", "w").write("\\newcommand{\\FitTable}{%\n" + "\n".join(t) + "}\n")
