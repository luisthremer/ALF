#!/usr/bin/env python3
"""Plots for the proposal from results.csv (analyse.py) - SciencePlots style.

usage: plot_all.py RESULTS.csv OUTDIR [--fake]
  --fake  stamps every figure with "FAKE DATA" (results_fake.csv from make_fake_data.py)

Figures (PDF + PNG):
  1 throughput_vs_L      best layout per side: sweeps per node-hour vs L, GPU/CPU factor annotated
  2 layout_scan          throughput of every layout, one panel per L
  3 walker_time          time per sweep of one walker vs layout (throughput vs single-chain speed)
  4 energy_per_sweep     estimated energy per sweep: GPU node vs CPU-only node (TDP sums)
  5 gpu_usage            mean tile utilisation and node GPU power vs walkers per tile
"""
import csv
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import scienceplots  # noqa: E402,F401

plt.style.use(["science", "no-latex"])

# Node power for the energy estimate (vendor TDP, to be checked against the data sheets before use in the proposal):
P_CPU_NODE = 2 * 350.0              # W, 2 x Xeon Platinum 8480L, notional CPU-only node
P_GPU_NODE = 4 * 600.0 + P_CPU_NODE  # W, + 4 x Data Center GPU Max 1550
SIDES = {"cpu": ("CPU node, 112 cores", "C0", "o"),
         "gpu": ("GPU node, 410c7b04", "C1", "s"),
         "gpub": ("GPU node, Intel-optimised", "C2", "^")}
POWER = {"cpu": P_CPU_NODE, "gpu": P_GPU_NODE, "gpub": P_GPU_NODE}


def load(path):
    """Read results.csv; repeats of one layout (same side, L, tag) are averaged, errors combined as sqrt(sum e^2)/n."""
    grp = defaultdict(list)
    for r in csv.DictReader(open(path)):
        if not r["throughput_node_h"] or r["rc"] not in ("0", 0):
            continue
        for k in ("L", "ranks", "threads"):
            r[k] = int(r[k])
        for k in ("t_sweep_s", "t_sweep_err_s", "throughput_node_h", "throughput_err"):
            r[k] = float(r[k])
        grp[(r["side"], r["L"], r["tag"])].append(r)
    rows = []
    for rs in grp.values():
        r = dict(rs[0]); n = len(rs)
        for k, e in (("t_sweep_s", "t_sweep_err_s"), ("throughput_node_h", "throughput_err")):
            r[k] = sum(x[k] for x in rs) / n
            r[e] = sum(x[e] ** 2 for x in rs) ** 0.5 / n
        r["n_rep"] = n
        rows.append(r)
    return rows


def stamp(fig, fake):
    if fake:
        fig.text(0.5, 0.5, "FAKE DATA", fontsize=28, color="red", alpha=0.25, ha="center", va="center", rotation=25)


def save(fig, out, name, fake):
    stamp(fig, fake)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=300)
    plt.close(fig)


def best(rows):
    b = {}
    for r in rows:
        k = (r["side"], r["L"])
        if k not in b or r["throughput_node_h"] > b[k]["throughput_node_h"]:
            b[k] = r
    return b


def fig_throughput(rows, out, fake):
    b = best(rows)
    Ls = sorted({r["L"] for r in rows})
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    for s, (lab, col, mk) in SIDES.items():
        pts = [b[(s, L)] for L in Ls if (s, L) in b]
        if pts:
            ax.errorbar([p["L"] for p in pts], [p["throughput_node_h"] for p in pts],
                        yerr=[p["throughput_err"] for p in pts], color=col, marker=mk, ms=4, capsize=2, label=lab)
    for L in Ls:
        if ("cpu", L) in b:
            g = max((b[(s, L)]["throughput_node_h"] for s in ("gpu", "gpub") if (s, L) in b), default=None)
            if g:
                ax.annotate(f"{g / b[('cpu', L)]['throughput_node_h']:.1f}x", (L, g), textcoords="offset points",
                            xytext=(0, 5), ha="center", fontsize=6)
    ax.set_yscale("log"); ax.set_xticks(Ls)
    ax.set_xlabel("$L$ ($\\beta = 4$)"); ax.set_ylabel("sweeps per node-hour")
    ax.legend(fontsize=6)
    save(fig, out, "throughput_vs_L", fake)


def fig_layout(rows, out, fake):
    Ls = sorted({r["L"] for r in rows})
    fig, axs = plt.subplots(1, len(Ls), figsize=(1.9 * len(Ls), 2.4), sharey=False)
    axs = axs if len(Ls) > 1 else [axs]
    for ax, L in zip(axs, Ls):
        sel = [r for r in rows if r["L"] == L]
        sel.sort(key=lambda r: (list(SIDES).index(r["side"]), r["tag"]))
        xs = range(len(sel))
        ax.bar(xs, [r["throughput_node_h"] for r in sel], yerr=[r["throughput_err"] for r in sel],
               color=[SIDES[r["side"]][1] for r in sel], capsize=1.5)
        ax.set_xticks(list(xs)); ax.set_xticklabels([r["tag"] for r in sel], rotation=90, fontsize=5)
        ax.set_title(f"$L={L}$", fontsize=7); ax.set_yscale("log")
    axs[0].set_ylabel("sweeps per node-hour")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c, _ in SIDES.values()]
    fig.legend(handles, [l for l, _, _ in SIDES.values()], loc="upper center", ncol=3, fontsize=6, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save(fig, out, "layout_scan", fake)


def fig_walker(rows, out, fake):
    Ls = sorted({r["L"] for r in rows})
    fig, axs = plt.subplots(1, 2, figsize=(5.0, 2.4))
    for L, ls in zip(Ls, ["-", "--", "-.", ":"]):
        for s, (lab, col, mk) in SIDES.items():
            sel = [r for r in rows if r["side"] == s and r["L"] == L and not (s != "cpu" and r["tag"] == "w4t1")]
            if not sel:
                continue
            x = [r["threads"] if s == "cpu" else int(r["w"]) for r in sel]
            o = sorted(range(len(x)), key=x.__getitem__)
            ax = axs[0] if s == "cpu" else axs[1]
            ax.errorbar([x[i] for i in o], [sel[i]["t_sweep_s"] for i in o], yerr=[sel[i]["t_sweep_err_s"] for i in o],
                        color=col, ls=ls, marker=mk, ms=3, capsize=1.5, label=f"{lab.split(',')[-1].strip()}, $L={L}$")
    axs[0].set_xlabel("OpenMP threads per rank"); axs[1].set_xlabel("walkers per tile")
    for ax in axs:
        ax.set_yscale("log"); ax.set_ylabel("time per sweep of one walker (s)")
    axs[1].set_xscale("log", base=2)
    axs[1].legend(fontsize=4.5, ncol=2)
    axs[0].legend(fontsize=4.5)
    fig.tight_layout()
    save(fig, out, "walker_time", fake)


def fig_energy(rows, out, fake):
    b = best(rows)
    Ls = sorted({r["L"] for r in rows})
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    for s, (lab, col, mk) in SIDES.items():
        pts = [b[(s, L)] for L in Ls if (s, L) in b]
        if pts:
            e = [POWER[s] * 3600 / p["throughput_node_h"] / 1000 for p in pts]   # kJ per sweep
            ee = [x * p["throughput_err"] / p["throughput_node_h"] for x, p in zip(e, pts)]
            ax.errorbar([p["L"] for p in pts], e, yerr=ee, color=col, marker=mk, ms=4, capsize=2, label=lab)
    ax.set_yscale("log"); ax.set_xticks(Ls)
    ax.set_xlabel("$L$ ($\\beta = 4$)"); ax.set_ylabel("energy per sweep (kJ, TDP estimate)")
    ax.legend(fontsize=6)
    ax.text(0.02, 0.02, f"CPU node {P_CPU_NODE:.0f} W, GPU node {P_GPU_NODE:.0f} W", transform=ax.transAxes, fontsize=5)
    save(fig, out, "energy_per_sweep", fake)


def fig_gpu(rows, out, fake):
    sel = [r for r in rows if r["side"] != "cpu" and r["gpu_util_pct"] and r["tag"] != "w4t1"]
    if not sel:
        return
    fig, axs = plt.subplots(1, 2, figsize=(5.0, 2.3))
    grp = defaultdict(list)
    for r in sel:
        grp[(r["side"], r["L"])].append(r)
    for (s, L), rs in sorted(grp.items()):
        rs.sort(key=lambda r: int(r["w"]))
        x = [int(r["w"]) for r in rs]
        ls = {16: "-", 24: "--", 32: "-.", 36: ":"}.get(L, "-")
        axs[0].plot(x, [float(r["gpu_util_pct"]) for r in rs], color=SIDES[s][1], ls=ls, marker=SIDES[s][2], ms=3,
                    label=f"{SIDES[s][0].split(',')[-1].strip()}, $L={L}$")
        axs[1].plot(x, [float(r["gpu_power_w"]) for r in rs], color=SIDES[s][1], ls=ls, marker=SIDES[s][2], ms=3)
    for ax in axs:
        ax.set_xscale("log", base=2); ax.set_xlabel("walkers per tile")
    axs[0].set_ylabel("mean GPU utilisation (%)"); axs[1].set_ylabel("GPU power, 4 cards (W)")
    axs[1].axhline(4 * 600, color="grey", lw=0.6, ls=":"); axs[1].text(1, 4 * 600, "4 x 600 W TDP", fontsize=5, va="bottom")
    axs[0].legend(fontsize=4.5, ncol=2)
    fig.tight_layout()
    save(fig, out, "gpu_usage", fake)


def table(rows, out):
    b = best(rows)
    Ls = sorted({r["L"] for r in rows})
    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write("| L | best CPU (layout, sweeps/node-h) | best GPU 410c7b04 | best GPU Intel-opt. | GPU/CPU (410c7b04) | GPU/CPU (Intel-opt.) |\n")
        f.write("|---|---|---|---|---|---|\n")
        for L in Ls:
            c, g, gb = b.get(("cpu", L)), b.get(("gpu", L)), b.get(("gpub", L))
            fmt = lambda r: f"{r['tag']}: {r['throughput_node_h']:.0f} ± {r['throughput_err']:.0f}" if r else "-"
            rat = lambda r: f"{r['throughput_node_h'] / c['throughput_node_h']:.1f}" if (r and c) else "-"
            f.write(f"| {L} | {fmt(c)} | {fmt(g)} | {fmt(gb)} | {rat(g)} | {rat(gb)} |\n")


def main():
    src, out = sys.argv[1], sys.argv[2]
    fake = "--fake" in sys.argv
    os.makedirs(out, exist_ok=True)
    rows = load(src)
    fig_throughput(rows, out, fake); fig_layout(rows, out, fake); fig_walker(rows, out, fake)
    fig_energy(rows, out, fake); fig_gpu(rows, out, fake); table(rows, out)
    print(f"wrote 5 figures + summary.md to {out}")


if __name__ == "__main__":
    main()
