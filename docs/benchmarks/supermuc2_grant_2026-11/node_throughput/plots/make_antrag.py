#!/usr/bin/env python3
"""Numbers and tables for the proposal document (report/antrag.tex) from results.csv.

usage: make_antrag.py RESULTS.csv > antrag_data.tex
Defines \\Lmax (largest L with CPU and both GPU results), \\FacIntel, \\FacBase, \\Eratio (at \\Lmax),
\\Complete (yes/no) and the table \\AntragTable (best layout per side and L, commit column, energy per sweep).
"""
import os
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from plot_all import POWER, SIDES, best, load  # noqa: E402

NAME = {"cpu": "CPU node (2 $\\times$ Xeon 8480L)", "gpu": "GPU node, standard GPU code", "gpub": "GPU node, Intel-optimised code"}
COMMIT = {"cpu": "410c7b04", "gpu": "410c7b04", "gpub": "d1060305"}
LAYOUT = {"cpu": lambda r: f"{r['ranks']} MPI proc.\\ $\\times$ {r['threads']} thread{'s' if r['threads'] > 1 else ''}",
          "gpu": lambda r: f"{r['w']} MPI proc./GPU tile, {r['threads']} thread{'s' if r['threads'] > 1 else ''}",
          "gpub": lambda r: f"{r['w']} MPI proc./GPU tile, {r['threads']} thread{'s' if r['threads'] > 1 else ''}"}
PLANNED = set(int(x) for x in os.environ.get("PLANNED_L", "16,24,32,36").split(","))


def num(x, nd=0):
    return f"{x:,.{nd}f}".replace(",", "\\,")


rows = load(sys.argv[1])
b = best(rows)
Ls = sorted({r["L"] for r in rows})
import os
_n = int(os.environ.get("NLAY", "5"))
NLAY = {"cpu": _n, "gpu": _n, "gpub": _n}   # configurations per side and L in the test
REQ = {k: set(v.split(",")) for k, v in (x.split("=") for x in os.environ.get("REQ_TAGS", "").split(";") if x)}   # e.g. cpu=t1,t2


NLAY_L = {int(k): int(v) for k, v in (x.split("=") for x in os.environ.get("NLAY_L", "").split(";") if x)}   # e.g. 40=1


def done(s_, L):
    tags = {r["tag"] for r in rows if r["side"] == s_ and r["L"] == L}
    return len(tags) >= NLAY_L.get(L, NLAY[s_]) and REQ.get(s_, set()) <= tags
full = [L for L in Ls if all((s, L) in b and done(s, L) for s in SIDES)]
out = []
if full:
    L = full[-1]
    c, g, gb = b[("cpu", L)], b[("gpu", L)], b[("gpub", L)]
    fi, fb = gb["throughput_node_h"] / c["throughput_node_h"], g["throughput_node_h"] / c["throughput_node_h"]
    er = (POWER["cpu"] / c["throughput_node_h"]) / (POWER["gpub"] / gb["throughput_node_h"])
    out += [f"\\newcommand{{\\Lmax}}{{{L}}}", f"\\newcommand{{\\FacIntel}}{{{fi:.1f}}}",
            f"\\newcommand{{\\FacBase}}{{{fb:.1f}}}", f"\\newcommand{{\\Eratio}}{{{er:.1f}}}"]
else:
    out += ["\\newcommand{\\Lmax}{--}", "\\newcommand{\\FacIntel}{--}", "\\newcommand{\\FacBase}{--}", "\\newcommand{\\Eratio}{--}"]
facs = [(L, b[("gpub", L)]["throughput_node_h"] / b[("cpu", L)]["throughput_node_h"]) for L in Ls if L in full]
out.append("\\newcommand{\\FacList}{" + ", ".join(f"{f:.1f}$\\times$ ($L={L}$)" for L, f in facs) + "}")
complete = set(full) >= PLANNED
partial = any(not done(s_, L) for s_ in SIDES for L in Ls if (s_, L) in b)   # a side without any run at that L (CPU above L = 40) is not partial
out.append("\\newcommand{\\StarNote}{" + ("$^*$\\,best of the configurations finished so far; no speedup given." if partial else "") + "}")
out.append(f"\\newcommand{{\\Complete}}{{{'yes' if complete else 'no'}}}")

SHORT = {"cpu": "CPU", "gpu": "GPU", "gpub": "GPU"}
CODE = {"cpu": "ALF CPU", "gpu": "ALF GPU", "gpub": "ALF GPU Intel-opt."}
CONF = {"cpu": lambda r: f"{r['ranks']} proc.\\ $\\times$ {r['threads']} thr.",
        "gpu": lambda r: f"{r['w']} proc./tile, {r['threads']} thr.",
        "gpub": lambda r: f"{r['w']} proc./tile, {r['threads']} thr."}
t = ["\\begin{tabular}{rlllrrr}",
     "$L$ & node & code & configuration & work units / node-h & speedup & kJ / unit\\\\", "\\hline"]
for i, L in enumerate(Ls):
    if i:
        t.append("\\hline")
    c = b.get(("cpu", L))
    for s_ in SIDES:
        r = b.get((s_, L))
        head = f"{L} & {SHORT[s_]} & {CODE[s_]}"
        if not r:
            t.append(f"{head} & \\multicolumn{{4}}{{l}}{{\\emph{{{'not measured' if s_ == 'cpu' else 'running'}}}}}\\\\")
            continue
        thr = r["throughput_node_h"]
        fac = f"{thr / c['throughput_node_h']:.1f}" if c else "--"
        e = POWER[s_] * 3600 / thr / 1000
        star = "" if done(s_, L) else "$^*$"
        if c and not done("cpu", L):
            fac = "--"
        t.append(f"{head} & {CONF[s_](r)}{star} & {num(thr)} & {fac} & {e:.3f}\\\\")
t.append("\\end{tabular}")
out.append("\\newcommand{\\AntragTable}{%\n" + "\n".join(t) + "}")
# power-law fit of the CPU node, throughput = A L^-p (p = CPU_EXP: 6 at beta = 4, 7 at beta = L, cost per sweep ~ N^3 Nt), amplitude by least squares in log space
import math
CP = float(os.environ.get("CPU_EXP", "6"))
cpuL = [L for L in Ls if ("cpu", L) in b]
AFIT = math.exp(sum(math.log(b[("cpu", L)]["throughput_node_h"]) + CP * math.log(L) for L in cpuL) / len(cpuL)) if cpuL else float("nan")
fit = lambda L: AFIT * L ** -CP
# data file for the pgfplots figures (best layout per side; nan = not finished)
with open(__file__.rsplit("/", 2)[0] + "/" + os.environ.get("REPORTDIR", "report") + "/antrag_best.dat", "w") as f:
    f.write("L thr_cpu thr_gpu thr_gpub e_cpu e_gpu e_gpub fac_gpub thr_cpufit\n")
    for L in Ls:
        v = []
        for s_ in SIDES:
            r = b.get((s_, L)); v.append(r["throughput_node_h"] if r else float("nan"))
        e = [POWER[s_] * 3600 / x / 1000 if x == x else float("nan") for s_, x in zip(SIDES, v)]
        fac = v[2] / v[0] if v[0] == v[0] and v[2] == v[2] else float("nan")
        f.write(f"{L} " + " ".join(f"{x:.6g}" for x in v + e + [fac, fit(L)]) + "\n")
with open(__file__.rsplit("/", 2)[0] + "/" + os.environ.get("REPORTDIR", "report") + "/antrag_fac.dat", "w") as f:   # speedup labels, only where both exist
    f.write("L thr_gpub fac\n")
    for L in Ls:
        c, g = b.get(("cpu", L)), b.get(("gpub", L))
        if c and g and done("cpu", L) and done("gpub", L):
            f.write(f"{L} {g['throughput_node_h']:.6g} {g['throughput_node_h'] / c['throughput_node_h']:.4g}\n")
# extrapolated speedups (no CPU measurement): GPU / fitted CPU
with open(__file__.rsplit("/", 2)[0] + "/" + os.environ.get("REPORTDIR", "report") + "/antrag_facx.dat", "w") as f:
    f.write("L thr_gpub fac\n")
    ex = [(L, b[("gpub", L)]["throughput_node_h"] / fit(L), b[("gpu", L)]["throughput_node_h"] / fit(L)) for L in Ls if ("cpu", L) not in b and ("gpub", L) in b and ("gpu", L) in b]
    for L, fi, fb in ex:
        f.write(f"{L} {b[('gpub', L)]['throughput_node_h']:.6g} {fi:.4g}\n")
if ex:
    out.append(f"\\newcommand{{\\FacExtrL}}{{{ex[-1][0]}}}\\newcommand{{\\FacExtr}}{{{ex[-1][1]:.0f}}}\\newcommand{{\\FacExtrBase}}{{{ex[-1][2]:.0f}}}\\newcommand{{\\CpuExp}}{{{int(CP)}}}")
    out.append("\\newcommand{\\FacExtrList}{" + ", ".join(f"{fi:.0f}$\\times$ ($L={L}$)" for L, fi, fb in ex) + "}")
else:
    out.append("\\newcommand{\\FacExtrL}{--}\\newcommand{\\FacExtr}{--}\\newcommand{\\FacExtrBase}{--}\\newcommand{\\CpuExp}{6}\\newcommand{\\FacExtrList}{--}")
print("\n".join(out))
