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
full = [L for L in Ls if all((s, L) in b for s in SIDES)]
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
out.append(f"\\newcommand{{\\Complete}}{{{'yes' if complete else 'no'}}}")

t = ["\\begin{tabular}{@{}rllrrr@{}}\\toprule",
     "$L$ & node, code (commit) & best configuration & work units per node-hour & factor vs CPU & energy per unit (kJ)\\\\\\midrule"]
for i, L in enumerate(Ls):
    if i:
        t.append("\\addlinespace")
    c = b.get(("cpu", L))
    for s in SIDES:
        r = b.get((s, L))
        lab = f"{NAME[s]} (\\texttt{{{COMMIT[s]}}})"
        if not r:
            t.append(f"{L} & {lab} & \\multicolumn{{4}}{{l}}{{\\emph{{run not finished}}}}\\\\")
            continue
        thr, err = r["throughput_node_h"], r["throughput_err"]
        fac = f"{thr / c['throughput_node_h']:.1f}" if c else "--"
        e = POWER[s] * 3600 / thr / 1000
        t.append(f"{L} & {lab} & {LAYOUT[s](r)} & ${num(thr)} \\pm {num(err)}$ & {fac} & {e:.3f}\\\\")
t.append("\\bottomrule\\end{tabular}")
out.append("\\newcommand{\\AntragTable}{%\n" + "\n".join(t) + "}")
print("\n".join(out))
