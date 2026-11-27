#!/usr/bin/env python3
"""Automatic search for the best GPU-node configuration per code and system size (SuperMUC-NG Phase 2), runs on the login node.

usage: tune.py search [--Ls 16,24,...] [--sides gpub,gpu]   beta = 4: layouts -> one parameter at a time -> combination -> final runs
       tune.py final  --beta L [--Ls ...]                    tuned settings of the beta = 4 search at another beta, with a layout fallback
       add --dry-run to print the planned runs only.
Procedure per (side, L), every run is one full node (jobs run through multi2.sh, one srun step per run):
  A  layouts (processes per tile x threads, disjoint core slices), defaults for everything else
  B  from the best layout, one parameter group at a time: N_delay, ASYNCBLK, BLKFUSE, QRNB (the last three Intel-optimised code only)
  C  winners of B (gain >= GAIN) combined in one run, plus a repeat of the baseline to measure the noise
  F  the two best configurations of A..C again with 10 timed sweeps; the faster one is the tuned setting
Results: $A/auto/tuned.json (side_L -> layout, env, throughput, alternatives), runs in $A/runs_auto (beta = 4), $A/runs_auto_betaL.
Metric: work units per node-hour (analyse.py). Noise is about 1-3 %, hence the GAIN threshold.
"""
import argparse
import csv
import io
import json
import os
import subprocess
import sys
import threading
import time

A = "/hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27"
J = A + "/jobs/nt"
LAYOUTS = ["2_7", "4_1", "4_2", "4_3", "6_2", "7_2", "8_1"]          # processes per tile _ threads (1_14 was never best)
GROUPS = {   # group -> [(variant, env settings)], first group value is the default (not run)
    "nd": [(f"nd{n}", f"NDELAY={n}") for n in (128, 256, 512, 1024)],
    "asy": [("asy0", "ALFSH_ASYNCBLK=0")],
    "bf": [("bf3", "ALFSH_BLKFUSE=3"), ("bf8", "ALFSH_BLKFUSE=8")],
    "qr": [("qr8", "ALFSH_QRNB=8"), ("qr32", "ALFSH_QRNB=32"), ("qr64", "ALFSH_QRNB=64")],
}
INTEL_ONLY = {"asy", "bf", "qr"}
GAIN = 0.03
BASE_ENV = "PINMODE=excl"
MAXPAIRS = 6    # (side, L) pairs searched at the same time

log_lock = threading.Lock()
sem = threading.Semaphore(MAXPAIRS)


def log(*a):
    with log_lock:
        print(time.strftime("%H:%M:%S"), *a, flush=True)


def nbin(L, beta):
    if beta == "L":
        return 3 if L >= 56 else 4
    return 4 if L >= 48 else 6


def walltime(L, beta):
    if beta == "L":
        return "02:00:00" if L >= 52 else "01:00:00"
    return "01:30:00" if L >= 48 else ("00:50:00" if L >= 36 else "00:30:00")


def run_items(items, name, L, beta, rundir, dry):
    """items: [(side, layout, nbin, variant, env)]; submits one job, waits, returns rows of analyse.py for this job"""
    args = [f"{s}:{L}:{lay}:{nb}:{var}:{BASE_ENV}{',' + env if env else ''}" for s, lay, nb, var, env in items]
    if dry:
        log("DRY", name, len(items), "nodes:", " ".join(args))
        return []
    cmd = ["sbatch", "--parsable", "--export=BETA,RUNDIR,EXE_GPU,BINB", "-N", str(len(items)), "-J", name, "-p", "general",
           "-t", walltime(L, beta), f"{J}/multi2.sh"] + args
    env = dict(os.environ, BETA=beta, RUNDIR=rundir, EXE_GPU=f"{A}/build_pp_alfgpu_detfix/Prog/ALF.out", BINB=f"{A}/bin_batched_65800b82/pp")
    jid = subprocess.run(cmd, env=env, cwd=J, capture_output=True, text=True, check=True).stdout.strip().split(";")[0]
    log("submitted", jid, name, len(items), "runs")
    while subprocess.run(["squeue", "-j", jid, "-h", "-o", "%T"], capture_output=True, text=True).stdout.strip():
        time.sleep(60)
    return [r for r in analyse(rundir) if r["job"] == jid]


def analyse(rundir):
    out = subprocess.run(["python3", f"{J}/analyse.py", f"{A}/{rundir}"], capture_output=True, text=True).stdout
    lines = [ln for ln in out.splitlines() if ln.startswith(("side,", "gpu", "gpub"))]
    return list(csv.DictReader(io.StringIO("\n".join(lines)))) if lines else []


def thr(r):
    return float(r["throughput_node_h"]) if r["rc"] == "0" and r["throughput_node_h"] else 0.0


def search_pair(side, L, beta, rundir, results, dry):
    with sem:
        nb = nbin(L, beta)
        N = L * L
        tagsel = lambda rows, tag: max([thr(r) for r in rows if r["tag"] == tag] or [0.0])
        # A: layouts
        rowsA = run_items([(side, lay, nb, "A", "") for lay in LAYOUTS], f"au_A_{side}_{L}", L, beta, rundir, dry)
        if dry:
            rowsA = [dict(tag=f"w{lay.split('_')[0]}t{lay.split('_')[1]}_A", rc="0", throughput_node_h="1", job="0") for lay in LAYOUTS]
        rankA = sorted(rowsA, key=thr, reverse=True)
        best = rankA[0]
        lay = best["tag"].rsplit("_", 1)[0]
        lay_arg = f"{best['w']}_{best['threads']}" if best.get("w") else LAYOUTS[0]
        base = thr(best)
        log(side, L, "A best", best["tag"], round(base), "ranking", [(r["tag"].split("_")[0], round(thr(r))) for r in rankA[:4]])
        # B: one group at a time
        items, groupof = [], {}
        for g, vs in GROUPS.items():
            if g in INTEL_ONLY and side != "gpub":
                continue
            for var, env in vs:
                if g == "nd" and int(var[2:]) > N:
                    continue
                items.append((side, lay_arg, nb, var, env)); groupof[var] = (g, env)
        rowsB = run_items(items, f"au_B_{side}_{L}", L, beta, rundir, dry) if items else []
        win = {}
        for r in rowsB:
            var = r["tag"].split("_", 1)[1]
            g, env = groupof[var]
            if thr(r) > base * (1 + GAIN) and thr(r) > win.get(g, (0, ""))[0]:
                win[g] = (thr(r), env, var)
        log(side, L, "B winners", {g: (v[2], round(v[0])) for g, v in win.items()})
        # C: combination + baseline repeat
        cand = [(base, lay_arg, "", "A")] + [(v[0], lay_arg, v[1], v[2]) for v in win.values()]
        itemsC = [(side, lay_arg, nb, "rep", "")]
        if len(win) > 1:
            itemsC.append((side, lay_arg, nb, "combo", ",".join(v[1] for v in win.values())))
        rowsC = run_items(itemsC, f"au_C_{side}_{L}", L, beta, rundir, dry)
        for r in rowsC:
            if r["tag"].endswith("_combo"):
                cand.append((thr(r), lay_arg, ",".join(v[1] for v in win.values()), "combo"))
        cand.sort(key=lambda c: c[0], reverse=True)
        # also keep the runner-up layouts for the fallback (memory) at other beta
        fall = [f"{r['w']}_{r['threads']}" for r in rankA[1:4] if r.get("w")]
        # F: the two best configurations with 10 timed sweeps
        fin = [(side, c[1], 11 if beta == "4" else nb, f"F{i}", c[2]) for i, c in enumerate(cand[:2])]
        rowsF = run_items(fin, f"au_F_{side}_{L}", L, beta, rundir, dry)
        bestF = max(rowsF, key=thr) if rowsF else None
        if bestF is None and not dry:
            bestF = best
        i = int(bestF["tag"].rsplit("_F", 1)[1]) if bestF is not None and "_F" in bestF["tag"] else 0
        choice = cand[min(i, len(cand) - 1)]
        results[f"{side}_{L}"] = dict(layout=choice[1], env=choice[2], variant=choice[3], throughput=thr(bestF) if bestF else 0,
                                      fallback_layouts=fall, base_env=BASE_ENV, beta=beta)
        log(side, L, "TUNED", results[f"{side}_{L}"])


def final_pair(side, L, beta, rundir, tuned, results, dry):
    with sem:
        t = tuned.get(f"{side}_{L}")
        if not t:
            log("no tuned setting for", side, L); return
        nb = nbin(L, beta)
        for lay in [t["layout"]] + t.get("fallback_layouts", []):
            rows = run_items([(side, lay, nb, "F", t["env"])], f"au_Fb_{side}_{L}", L, beta, rundir, dry)
            if dry or (rows and thr(rows[0]) > 0):
                results[f"{side}_{L}"] = dict(t, layout=lay, throughput=thr(rows[0]) if rows else 0, beta=beta)
                log(side, L, beta, "FINAL", results[f"{side}_{L}"]); return
            log(side, L, beta, "layout", lay, "failed, next")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["search", "final"])
    ap.add_argument("--Ls", default="16,24,32,36,40,44,48,52,56")
    ap.add_argument("--sides", default="gpub,gpu")
    ap.add_argument("--beta", default="4")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    if not a.dry_run:
        os.makedirs(f"{A}/auto", exist_ok=True)
    jf = f"{A}/auto/tuned.json"
    tuned = json.load(open(jf)) if os.path.exists(jf) else {}
    rundir = "runs_auto" if a.beta == "4" else f"runs_auto_beta{a.beta}"
    results, threads = {}, []
    for side in a.sides.split(","):
        for L in [int(x) for x in a.Ls.split(",")]:
            tgt = (lambda s=side, L=L: search_pair(s, L, a.beta, rundir, results, a.dry_run)) if a.mode == "search" else \
                  (lambda s=side, L=L: final_pair(s, L, a.beta, rundir, tuned, results, a.dry_run))
            th = threading.Thread(target=tgt); th.start(); threads.append(th)
    for th in threads:
        th.join()
    if not a.dry_run:
        if a.mode == "search":
            tuned.update(results); json.dump(tuned, open(jf, "w"), indent=1)
        else:
            json.dump(results, open(f"{A}/auto/final_beta{a.beta}.json", "w"), indent=1)
    log("done", len(results), "settings")


if __name__ == "__main__":
    main()
