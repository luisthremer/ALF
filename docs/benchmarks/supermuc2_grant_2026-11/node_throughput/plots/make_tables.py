#!/usr/bin/env python3
"""LaTeX tables from results.csv (same averaging over repeats as plot_all.py).

usage: make_tables.py RESULTS.csv > tables.tex
Writes: one table per L with every layout, and the summary table (best layout per side, GPU/CPU factor).
"""
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from plot_all import SIDES, best, load  # noqa: E402

NAME = {"cpu": "CPU", "gpu": "GPU 410c7b04", "gpub": "GPU Intel-opt."}


def pm(v, e, nd=0):
    return f"${v:,.{nd}f} \\pm {e:,.{nd}f}$".replace(",", "\\,")


def layout(r):
    if r["side"] == "cpu":
        return f"{r['ranks']} ranks $\\times$ {r['threads']} thr."
    return f"{r['w']} walker/tile, {r['threads']} thr."


rows = load(sys.argv[1])
Ls = sorted({r["L"] for r in rows})
out = []
b = best(rows)
out.append("\\begin{table}[h]\\centering\\small")
out.append("\\caption{Best layout per side: node throughput (sweeps per node-hour) and GPU/CPU factor.}\\label{tab:summary}")
out.append("\\resizebox{\\textwidth}{!}{\\begin{tabular}{@{}rlllrr@{}}\\toprule")
out.append("$L$ & CPU (best) & GPU 410c7b04 (best) & GPU Intel-opt.\\ (best) & factor 410c7b04 & factor Intel-opt.\\\\\\midrule")
for L in Ls:
    c, g, gb = b.get(("cpu", L)), b.get(("gpu", L)), b.get(("gpub", L))
    cell = lambda r: f"\\texttt{{{r['tag']}}}: {pm(r['throughput_node_h'], r['throughput_err'])}" if r else "--"
    fac = lambda r: f"{r['throughput_node_h'] / c['throughput_node_h']:.2f}" if (r and c) else "--"
    out.append(f"{L} & {cell(c)} & {cell(g)} & {cell(gb)} & {fac(g)} & {fac(gb)}\\\\")
out.append("\\bottomrule\\end{tabular}}\\end{table}")
for L in Ls:
    sel = sorted([r for r in rows if r["L"] == L], key=lambda r: (list(SIDES).index(r["side"]), r["tag"]))
    out.append("\\begin{table}[h]\\centering\\small")
    out.append(f"\\caption{{All layouts at $L={L}$ ($\\beta=4$). Time per sweep of one walker (mean of 10 timed sweeps, "
               f"$\\pm$ standard error; averaged over repeats) and node throughput.}}")
    out.append("\\begin{tabular}{@{}lllrr@{}}\\toprule")
    out.append("side & tag & layout & time per sweep (s) & sweeps per node-hour\\\\\\midrule")
    for r in sel:
        out.append(f"{NAME[r['side']]} & \\texttt{{{r['tag']}}} & {layout(r)} & {pm(r['t_sweep_s'], r['t_sweep_err_s'], 3)} & "
                   f"{pm(r['throughput_node_h'], r['throughput_err'])}\\\\")
    out.append("\\bottomrule\\end{tabular}\\end{table}")
print("\n".join(out))
