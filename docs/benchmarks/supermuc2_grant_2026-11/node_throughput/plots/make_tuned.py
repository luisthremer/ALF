#!/usr/bin/env python3
"""LaTeX table of the tuned settings: usage make_tuned.py tuned.json [final_betaL.json] > tuned.tex   (macro \\TunedTable)"""
import json, sys
NAME = {"gpu": "ALF GPU", "gpub": "ALF GPU Intel-opt."}
d = json.load(open(sys.argv[1]))
fb = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else {}
def lay(x): w, t = x.split("_"); return f"{w} proc./tile $\\times$ {t} thr."
def env(e): return (e.replace("NDELAY=", "$N_\\text{delay}=$").replace("ALFSH_", "").replace("_", "\\_").replace(",", ", ")) or "default"
rows = ["\\begin{tabular}{rllll}", "$L$ & code & layout $\\beta=4$ & settings $\\beta=4$ & layout $\\beta=L$\\\\", "\\hline"]
Ls = sorted({int(k.split("_")[1]) for k in d})
for L in Ls:
    for s in ("gpub", "gpu"):
        t = d.get(f"{s}_{L}")
        if not t:
            continue
        b = fb.get(f"{s}_{L}")
        rows.append(f"{L} & {NAME[s]} & {lay(t['layout'])} & {env(t['env'])} & {lay(b['layout']) if b else '--'}\\\\")
    rows.append("\\hline")
rows.append("\\end{tabular}")
print("\\newcommand{\\TunedTable}{%\n" + "\n".join(rows) + "}")
