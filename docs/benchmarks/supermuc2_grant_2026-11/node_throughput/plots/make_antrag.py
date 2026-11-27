#!/usr/bin/env python3
"""Numbers and tables for the proposal document (report/antrag.tex) from results.csv.

usage: make_antrag.py RESULTS.csv > antrag_data.tex
Defines \\Lmax (largest L with CPU and both GPU results), \\FacIntel, \\FacBase, \\Eratio (at \\Lmax),
\\Complete (yes/no) and the table \\AntragTable (best layout per side and L, commit column, energy per sweep).
"""
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from plot_all import POWER, SIDES, best, load  # noqa: E402

NAME = {"cpu": "CPU node (2 $\\times$ Xeon 8480L)", "gpu": "GPU node, standard GPU code", "gpub": "GPU node, Intel-optimised code"}
COMMIT = {"cpu": "410c7b04", "gpu": "410c7b04", "gpub": "d1060305"}
LAYOUT = {"cpu": lambda r: f"{r['ranks']} MPI proc.\\ $\\times$ {r['threads']} thread{'s' if r['threads'] > 1 else ''}",
          "gpu": lambda r: f"{r['w']} MPI proc./GPU tile, {r['threads']} thread{'s' if r['threads'] > 1 else ''}",
          "gpub": lambda r: f"{r['w']} MPI proc./GPU tile, {r['threads']} thread{'s' if r['threads'] > 1 else ''}"}
PLANNED = {16, 24, 32, 36}


def num(x, nd=0):
    return f"{x:,.{nd}f}".replace(",", "\\,")


rows = load(sys.argv[1])
b = best(rows)
Ls = sorted({r["L"] for r in rows})
import os
_n = int(os.environ.get("NLAY", "5"))
NLAY = {"cpu": _n, "gpu": _n, "gpub": _n}   # configurations per side and L in the test
done = lambda s_, L: len({r["tag"] for r in rows if r["side"] == s_ and r["L"] == L}) >= NLAY[s_]
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
complete = set(full) >= PLANNED
partial = any(not done(s_, L) for s_ in SIDES for L in Ls)
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
            t.append(f"{head} & \\multicolumn{{4}}{{l}}{{\\emph{{running}}}}\\\\")
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
# data file for the pgfplots figures (best layout per side; nan = not finished)
with open(__file__.rsplit("/", 2)[0] + "/" + os.environ.get("REPORTDIR", "report") + "/antrag_best.dat", "w") as f:
    f.write("L thr_cpu thr_gpu thr_gpub e_cpu e_gpu e_gpub fac_gpub\n")
    for L in Ls:
        v = []
        for s_ in SIDES:
            r = b.get((s_, L)); v.append(r["throughput_node_h"] if r else float("nan"))
        e = [POWER[s_] * 3600 / x / 1000 if x == x else float("nan") for s_, x in zip(SIDES, v)]
        fac = v[2] / v[0] if v[0] == v[0] and v[2] == v[2] else float("nan")
        f.write(f"{L} " + " ".join(f"{x:.6g}" for x in v + e + [fac]) + "\n")
with open(__file__.rsplit("/", 2)[0] + "/" + os.environ.get("REPORTDIR", "report") + "/antrag_fac.dat", "w") as f:   # speedup labels, only where both exist
    f.write("L thr_gpub fac\n")
    for L in Ls:
        c, g = b.get(("cpu", L)), b.get(("gpub", L))
        if c and g and done("cpu", L) and done("gpub", L):
            f.write(f"{L} {g['throughput_node_h']:.6g} {g['throughput_node_h'] / c['throughput_node_h']:.4g}\n")
print("\n".join(out))
