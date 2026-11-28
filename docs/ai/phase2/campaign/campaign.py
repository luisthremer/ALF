#!/usr/bin/env python3
"""Data campaign harness (jobs 1 and 2 of brief-sim-campaign).  No Prog/ code.  Reuses p25_run.py (prepare/run) and gpu_e2e.py (numdiff).

  tol    GPU path vs CPU path of the same binary (ALF_GPU_DISABLE=1), max scaled difference |x-y|/max(1,|y|) per observable class
         (scal = *_scal, eq = *_eq, tau = *_tau, other), against size, per model and scheme; null pairs on the CPU path:
         null_nw = CPU Nwrap 10 vs 5 (std scheme), null_sc = CPU scheme vs CPU standard.  Resumable (dir/.rc marks a finished run).
  bench  seconds per sweep = (t(NSweep=6) - t(NSweep=2)) / 4, Ltau = 0, per model / size / scheme, GPU path or CPU path.
usage: campaign.py tol|bench --exe EXE --machine NAME --work DIR --out DIR [--cpu] [--jobs N] [--sizes ...]"""
import argparse, os, re, sys, time
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "testsuite", "GPU.tests")); sys.path.insert(0, os.path.join(REPO, "docs", "ai", "phase2", "p25_runs"))
import gpu_e2e as e2e, p25_run as p25

KONDO = {"ham_name": '"Kondo"', "Model": '"Kondo"', "Lattice_type": '"Bilayer_square"', "ham_Uf": "1.d0"}
TV = {"ham_name": '"tV"', "Model": '"tV"', "N_SUN": "1", "ham_V": "-1.d0", "ham_T2": "0.d0", "ham_Tperp": "0.d0", "Beta": "2.d0", "Ltau": "0"}
SP = {"ham_name": '"Spin_Peierls"', "Model": '"Spin_Peierls"', "Beta": "2.d0", "ham_U": "1.5d0", "ham_Lambda": "0.1d0",
      "ham_omega0": "0.5d0", "N_SUN": "2", "N_FL": "1"}
MODELS = {"hub-mz": {}, "hub-su2": {"Mz": ".F."}, "kondo": KONDO, "tV": TV, "sp": SP}
SCHEMES = {"std": {}, "dly16": {"Sequential_scheme": '"delayed"', "N_delay": "16"}, "sub16": {"Sequential_scheme": '"submatrix"', "N_delay": "16"}}
TOL_SIZES = {"hub-mz": [8, 12, 16, 24, 32], "hub-su2": [8, 12, 16, 24, 32], "kondo": [4, 6, 8, 12], "tV": [8, 12, 16], "sp": [8, 12, 16]}
BENCH_SIZES = {"hub-mz": [8, 16, 24], "kondo": [6, 8, 12], "tV": [8, 12, 16], "sp": [8, 12, 16]}

def edits(model, L, scheme, extra=None):
    e = {"L1": str(L), "L2": str(L)}; e.update(MODELS[model]); e.update(SCHEMES[scheme]); e.update(extra or {}); return e

def do(exe, d, start, ed, env, timeout):
    """run once (resumable); returns (rc, seconds)"""
    mark = os.path.join(d, ".rc")
    if os.path.exists(mark):
        rc, t = open(mark).read().split(); return int(rc), float(t)
    p25.prepare(d, start, ed)
    rc, t = p25.run(exe, d, env, 0, timeout)
    open(mark, "w").write("%d %.3f" % (rc, t)); return rc, t

def state(d):
    l = p25.gpu_lines(d)
    if any("ACTIVE (mirror" in x for x in l): return "ACTIVE"
    r = next((x.split("reason:", 1)[1].strip() for x in l if "INACTIVE" in x and "reason:" in x), "")
    return "INACTIVE:" + r if r or any("INACTIVE" in x for x in l) else "NOINFO"

def cls(f):
    return "scal" if f.endswith("_scal") else "eq" if f.endswith("_eq") else "tau" if f.endswith("_tau") else None

def fdiffs(a, b):
    """{file: max scaled diff or nan}, confout identical?"""
    import filecmp
    out, conf = {}, True
    for f in e2e.files(a):
        fa, fb = os.path.join(a, f), os.path.join(b, f)
        if f.endswith(".h5") or not os.path.exists(fb) or (cls(f) is None and not f.startswith("confout")): continue
        if f.startswith("confout"):
            conf = conf and filecmp.cmp(fa, fb, shallow=False); continue
        if filecmp.cmp(fa, fb, shallow=False): out[f] = 0.0; continue
        dd = e2e.numdiff(fa, fb); out[f] = float("nan") if dd is None else dd[1]
    return out, conf

def cmax(df):
    r = {}
    for f, v in df.items(): r[cls(f)] = max(r.get(cls(f), 0.0), v) if v == v else float("nan")
    return r

def tol(a):
    jobs = []; cells = []
    models = a.models.split(",")
    for m in models:
        for L in [int(x) for x in a.sizes.split(",")] if a.sizes else TOL_SIZES[m]:
            if a.maxL and L > a.maxL: continue
            base = os.path.join(a.work, "%s-L%d" % (m, L)); cells.append((m, L, base))
            jobs.append((os.path.join(base, "cpu-nw5"), edits(m, L, "std", {"Nwrap": "5"}), {"ALF_GPU_DISABLE": "1"}))
            for s in SCHEMES:
                jobs.append((os.path.join(base, "gpu-" + s), edits(m, L, s), {}))
                jobs.append((os.path.join(base, "cpu-" + s), edits(m, L, s), {"ALF_GPU_DISABLE": "1"}))
    ex = dict(a.extra)
    for _, e, _ in jobs: e.update({"NSweep": a.nsweep, "NBin": "1"}); e.update(ex)
    jobs.sort(key=lambda j: -int(j[1]["L1"]))     # big first: balance the pool
    with ThreadPoolExecutor(a.jobs) as pool:
        futs = [pool.submit(do, a.exe, d, a.start, e, env, a.timeout) for d, e, env in jobs]
        for f in futs: f.result()
    rows, files = [], []
    for m, L, base in cells:
        def rc(n): return int(open(os.path.join(base, n, ".rc")).read().split()[0])
        if rc("cpu-std") != 0 or rc("cpu-nw5") != 0:
            rows.append([a.machine, m, "-", L, "-", "CPU RUN FAILED rc=%d/%d" % (rc("cpu-std"), rc("cpu-nw5"))] + [""] * 4); continue
        nw, _ = fdiffs(os.path.join(base, "cpu-std"), os.path.join(base, "cpu-nw5")); nwc = cmax(nw)
        pgc = e2e.green_precision(os.path.join(base, "cpu-std")); 
        for s in SCHEMES:
            g, c = os.path.join(base, "gpu-" + s), os.path.join(base, "cpu-" + s)
            if rc("gpu-" + s) != 0 or rc("cpu-" + s) != 0:
                rows.append([a.machine, m, s, L, "-", "RUN FAILED gpu rc=%d cpu rc=%d" % (rc("gpu-" + s), rc("cpu-" + s))] + [""] * 4); continue
            st = state(g); df, conf = fdiffs(g, c); gc = cmax(df)
            sc = cmax(fdiffs(c, os.path.join(base, "cpu-std"))[0]) if s != "std" else {}
            pg = e2e.green_precision(g); pc = e2e.green_precision(c)
            for k in sorted(set(gc) | set(nwc)):
                rows.append([a.machine, m, s, L, k, "%.3e" % gc.get(k, float("nan")), "%.3e" % nwc.get(k, float("nan")),
                             ("%.3e" % sc[k]) if k in sc else "-", "%s;confout_%s;PG gpu %s cpu %s" % (st, "id" if conf else "DIFF",
                             "%.2e" % pg[0] if pg else "-", "%.2e" % pc[0] if pc else "-")])
            for f, v in sorted(df.items()): files.append([a.machine, m, s, L, f, "%.3e" % v, "%.3e" % nw.get(f, float("nan"))])
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "tol_%s.tsv" % a.machine), "w") as o:
        o.write("machine\tmodel\tscheme\tL\tclass\tgpu_cpu_diff\tnull_nwrap\tnull_scheme\tstate;confout;PG means (mean of 'Precision Green')\n")
        for r in rows: o.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(a.out, "tol_%s_files.tsv" % a.machine), "w") as o:
        o.write("machine\tmodel\tscheme\tL\tfile\tgpu_cpu_diff\tnull_nwrap\n")
        for r in files: o.write("\t".join(map(str, r)) + "\n")
    print("tol done", len(rows), "rows")

def bench(a):
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, "bench_%s.tsv" % a.machine); new = not os.path.exists(path)
    o = open(path, "a")
    if new: o.write("machine\tpath\tmodel\tL\tscheme\tstate\tt_2sweeps_s\tt_6sweeps_s\ts_per_sweep\n"); o.flush()
    env = {"ALF_GPU_DISABLE": "1"} if a.cpu else {}
    for m in a.models.split(","):
        for L in [int(x) for x in a.sizes.split(",")] if a.sizes else BENCH_SIZES[m]:
            for s in SCHEMES:
                t = {}
                for n in (2, 6):
                    d = os.path.join(a.work, "%s-L%d-%s-n%d" % (m, L, s, n))
                    rc, tt = do(a.exe, d, a.start, edits(m, L, s, {"NSweep": str(n), "NBin": "1", "Ltau": "0"}), env, a.timeout)
                    t[n] = tt if rc == 0 else None
                    stt = state(d)
                    if rc != 0: stt = "rc=%d %s" % (rc, stt)
                sps = "%.3f" % ((t[6] - t[2]) / 4) if t[2] is not None and t[6] is not None else "-"
                o.write("\t".join([a.machine, "cpu" if a.cpu else "gpu", m, str(L), s, stt] + ["%.2f" % t[n] if t[n] else "-" for n in (2, 6)] + [sps]) + "\n"); o.flush()
    print("bench done")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["tol", "bench"]); ap.add_argument("--exe", required=True); ap.add_argument("--machine", required=True)
    ap.add_argument("--work", required=True); ap.add_argument("--out", required=True); ap.add_argument("--cpu", action="store_true")
    ap.add_argument("--start", default=os.path.join(REPO, "Scripts_and_Parameters_files", "Start"))
    ap.add_argument("--jobs", type=int, default=1); ap.add_argument("--models", default="hub-mz,hub-su2,kondo,sp,tV")
    ap.add_argument("--sizes", default=""); ap.add_argument("--maxL", type=int, default=0); ap.add_argument("--nsweep", default="2")
    ap.add_argument("--timeout", type=int, default=3600); ap.add_argument("--extra", default="", help="k=v,k=v extra parameter edits")
    a = ap.parse_args(); a.extra = dict(x.split("=", 1) for x in a.extra.split(",") if x)
    if a.mode == "bench": a.models = a.models.replace("hub-su2,", "").replace(",hub-su2", "")
    {"tol": tol, "bench": bench}[a.mode](a)
