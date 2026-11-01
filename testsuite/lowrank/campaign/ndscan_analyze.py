#!/usr/bin/env python3
"""N_delay scan analysis: sweep time vs N_delay per parameter point, optimum, automatic choice, plots.

Usage: python3 -I testsuite/lowrank/campaign/ndscan_analyze.py COLLECTED.json [COLLECTED.json ...] --out DIR [--groups FILE]

Input: collect.py output of campaigns whose runs vary Sequential_scheme and N_delay per point (smuc_ndscan). Every number comes from
the job logs: per run, the sweep time is the mean over the bins >= 2 of the T1 timers' `sweep` (mean over ranks) divided by the sweeps of
that bin (bin 1 holds allocation/first-touch/warm-up); runs without timers fall back to collect.py's sec_per_sweep (flagged).
Per (point, layout, scheme, N_delay requested): median over the repetitions, min and max (range).
Optimum: the fixed N_delay with the smallest median; "within noise" are all fixed values with median <= best * (1 + tol),
tol = max(2 %, 2 x the larger relative range of the two). Automatic (N_delay 0) and automatic+timing (-1): the value(s) used, the
median time and the ratio to the best fixed value.
Output in DIR: ndscan.md (tables), ndscan_vs_ndelay.png/.svg (relative time vs N_delay, one panel per parameter group),
ndscan_abs.png (absolute s/sweep incl. the standard scheme as a horizontal line), ndscan_rows.tsv (per-run values).
Groups: a JSON file {"group title": ["point", ...], ...}; default: groups by point-name prefix as in smuc_ndscan.
"""
import argparse
import json
import math
import os
import statistics
import sys

DEFAULT_GROUPS = {
    "System size (Hubbard U=4, Beta=2, 2x24)": ["hubL20", "hubL30", "hubL40", "hubL50"],
    "Time slices (Hubbard L=40, U=4)": ["hubL40b1", "hubL40", "hubL40b4"],
    "Interaction U (Hubbard L=40, Beta=2)": ["hubL40U1", "hubL40", "hubL40U8"],
    "Vertex rank / model (Beta=2)": ["hubL20", "tvL20", "kondoL20", "hubL40", "tvL40", "tvL30"],
    "Threads per rank (Hubbard L=40)": ["hubL40@1x48", "hubL40", "hubL40@4x12"],
}


def run_time(r):
    """(seconds per sweep, source) of one collected row."""
    t = r.get("timers") or {}
    per_bin = t.get("timers") or {}
    sweeps = t.get("sweeps") or {}
    vals = []
    for b, cats in per_bin.items():
        try:
            bi = int(b)
        except ValueError:
            continue
        n = sweeps.get(b) or sweeps.get(str(bi))
        if bi >= 2 and n and "sweep" in cats:
            vals.append(cats["sweep"]["mean"] / float(n))
    if vals:
        return sum(vals) / len(vals), "timers"
    return float(r["sec_per_sweep"]), "cpu_time"


def key_point(r):
    lay = r["layout"]
    return r["point"] if lay == "2x24" else "%s@%s" % (r["point"], lay)


def load(paths):
    rows = []
    for p in paths:
        d = json.load(open(p))
        for r in d["rows"]:
            if not str(r.get("status", "OK")).startswith("OK"):
                continue
            t, src = run_time(r)
            ndr = r.get("n_delay_requested")
            ndr = None if ndr in (None, "", "n/a") else int(ndr)
            ndu = r.get("n_delay_used")
            ndu = None if ndu in (None, "", "n/a", "na") else int(ndu)
            rows.append(dict(point=key_point(r), scheme=r["scheme"], ndr=ndr, ndu=ndu, rep=r["repetition"], t=t, src=src,
                             L=r["L"], ndim=r.get("ndim"), acc=r.get("acceptance"), campaign=d.get("campaign_id"),
                             run_id=r["run_id"], n0=r.get("n_delay_n0")))
    return rows


def summarize(rows):
    g = {}
    for r in rows:
        k = (r["point"], r["scheme"], r["ndr"] if r["scheme"] != "standard" else None)
        g.setdefault(k, []).append(r)
    s = {}
    for k, rs in g.items():
        ts = [r["t"] for r in rs]
        med = statistics.median(ts)
        s[k] = dict(n=len(ts), med=med, lo=min(ts), hi=max(ts), rel=(max(ts) - min(ts)) / med if med > 0 else 0.0,
                    used=sorted({r["ndu"] for r in rs if r["ndu"] is not None}), rows=rs)
    return s


def analyse(s):
    points = sorted({k[0] for k in s})
    res = []
    for p in points:
        std = s.get((p, "standard", None))
        for sc in ("delayed", "submatrix"):
            fixed = {k[2]: v for k, v in s.items() if k[0] == p and k[1] == sc and k[2] is not None and k[2] >= 1}
            if not fixed:
                continue
            best_nd = min(fixed, key=lambda n: fixed[n]["med"])
            best = fixed[best_nd]
            within = []
            for n, v in sorted(fixed.items()):
                tol = max(0.02, 2 * max(best["rel"], v["rel"]))
                if v["med"] <= best["med"] * (1 + tol):
                    within.append(n)
            entry = dict(point=p, scheme=sc, best_nd=best_nd, best=best, within=within, fixed=fixed, std=std)
            for req, name in ((0, "auto"), (-1, "tune")):
                v = s.get((p, sc, req))
                if v:
                    entry[name] = dict(v, ratio=v["med"] / best["med"])
            res.append(entry)
    return res


def fmt_t(v):
    return "%.4g" % v


def write_md(res, s, out, rows):
    L = []
    L.append("# N_delay scan: sweep time vs N_delay\n")
    L.append("Generated by `testsuite/lowrank/campaign/ndscan_analyze.py` from collect.py output (%s). Sweep time = timers, bins >= 2, "
             "mean over ranks; median over repetitions [min, max]. Within noise: median <= best x (1 + max(2 %%, 2 x relative range)).\n"
             % ", ".join(sorted({r["campaign"] for r in rows})))
    nsrc = sum(1 for r in rows if r["src"] != "timers")
    if nsrc:
        L.append("**%d runs without timers: CPU-time fallback (includes bin 1).**\n" % nsrc)
    L.append("## Optimum and automatic choice\n")
    L.append("| point | scheme | best fixed N_delay | best s/sweep [range] | within noise | auto: used, ratio to best | tune: used, ratio | standard s/sweep | speedup best vs std |")
    L.append("|---|---|---:|---|---|---|---|---:|---:|")
    for e in res:
        b = e["best"]
        a = e.get("auto"); tn = e.get("tune")
        at = "%s, %.3f" % ("/".join(map(str, a["used"])), a["ratio"]) if a else "-"
        tt = "%s, %.3f" % ("/".join(map(str, tn["used"])), tn["ratio"]) if tn else "-"
        stdm = e["std"]["med"] if e["std"] else None
        L.append("| %s | %s | %d | %s [%s, %s] (n=%d) | %s | %s | %s | %s | %s |" % (
            e["point"], e["scheme"], e["best_nd"], fmt_t(b["med"]), fmt_t(b["lo"]), fmt_t(b["hi"]), b["n"],
            ", ".join(map(str, e["within"])), at, tt, fmt_t(stdm) if stdm else "-",
            "%.2f" % (stdm / b["med"]) if stdm else "-"))
    L.append("\n## All values (median s/sweep, n repetitions)\n")
    nds = sorted({k[2] for k in s if k[2] is not None})
    L.append("| point | scheme | " + " | ".join(("auto" if n == 0 else "tune" if n == -1 else str(n)) for n in nds) + " |")
    L.append("|---|---|" + "---:|" * len(nds))
    for e in res:
        cells = []
        for n in nds:
            v = s.get((e["point"], e["scheme"], n))
            cells.append("%s (%d)" % (fmt_t(v["med"]), v["n"]) if v else "-")
        L.append("| %s | %s | %s |" % (e["point"], e["scheme"], " | ".join(cells)))
    open(os.path.join(out, "ndscan.md"), "w").write("\n".join(L) + "\n")


def plot(res, groups, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    byp = {}
    for e in res:
        byp.setdefault(e["point"], {})[e["scheme"]] = e
    gl = [(t, [p for p in ps if p in byp]) for t, ps in groups.items()]
    gl = [(t, ps) for t, ps in gl if ps]
    if not gl:
        print("ndscan_analyze: no group has data, no plot", file=sys.stderr)
        return
    for kind in ("rel", "abs"):
        n = len(gl)
        cols = min(3, n)
        rws = int(math.ceil(n / float(cols)))
        fig, axes = plt.subplots(rws, cols, figsize=(5.2 * cols, 4.0 * rws), squeeze=False)
        for i, (title, ps) in enumerate(gl):
            ax = axes[i // cols][i % cols]
            colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
            for j, p in enumerate(ps):
                c = colors[j % len(colors)]
                for sc, ls, mk in (("delayed", "-", "o"), ("submatrix", "--", "s")):
                    e = byp[p].get(sc)
                    if not e:
                        continue
                    xs = sorted(e["fixed"])
                    norm = e["best"]["med"] if kind == "rel" else 1.0
                    ys = [e["fixed"][x]["med"] / norm for x in xs]
                    lo = [e["fixed"][x]["med"] / norm - e["fixed"][x]["lo"] / norm for x in xs]
                    hi = [e["fixed"][x]["hi"] / norm - e["fixed"][x]["med"] / norm for x in xs]
                    ax.errorbar(xs, ys, yerr=[lo, hi], color=c, ls=ls, marker=mk, ms=4, capsize=2,
                                label="%s %s" % (p, sc[:3]))
                    a = e.get("auto")
                    if a and a["used"]:
                        ax.plot([a["used"][0]], [a["med"] / norm], marker="*", ms=13, color=c, mec="k", ls="none")
                    if kind == "abs" and sc == "delayed" and e["std"]:
                        ax.axhline(e["std"]["med"], color=c, lw=0.8, ls=":")
            ax.set_xscale("log", base=2)
            ax.set_xlabel("N_delay")
            ax.set_ylabel("sweep time / best fixed" if kind == "rel" else "s per sweep")
            if kind == "abs":
                ax.set_yscale("log")
            ax.set_title(title, fontsize=9)
            ax.grid(True, which="both", alpha=0.3)
            ax.legend(fontsize=6, ncol=2)
        for k in range(len(gl), rws * cols):
            axes[k // cols][k % cols].axis("off")
        fig.suptitle("N_delay scan (SuperMUC-NG)\nsolid delayed, dashed submatrix, star = automatic N_delay"
                     + ("" if kind == "rel" else ", dotted = standard scheme"), fontsize=10)
        fig.tight_layout()
        base = os.path.join(out, "ndscan_vs_ndelay" if kind == "rel" else "ndscan_abs")
        fig.savefig(base + ".png", dpi=130)
        if kind == "rel":
            fig.savefig(base + ".svg")
        plt.close(fig)


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("collected", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--groups")
    a = ap.parse_args(argv[1:])
    os.makedirs(a.out, exist_ok=True)
    rows = load(a.collected)
    if not rows:
        print("ndscan_analyze: no OK rows", file=sys.stderr)
        return 1
    with open(os.path.join(a.out, "ndscan_rows.tsv"), "w") as fh:
        fh.write("run_id\tpoint\tscheme\tn_delay_req\tn_delay_used\trep\tsec_per_sweep\tsource\n")
        for r in rows:
            fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%.6g\t%s\n" % (r["run_id"], r["point"], r["scheme"], r["ndr"], r["ndu"], r["rep"], r["t"], r["src"]))
    s = summarize(rows)
    res = analyse(s)
    write_md(res, s, a.out, rows)
    groups = json.load(open(a.groups)) if a.groups else DEFAULT_GROUPS
    plot(res, groups, a.out)
    print("ndscan_analyze: %d runs, %d (point, scheme) curves -> %s" % (len(rows), len(res), a.out))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
