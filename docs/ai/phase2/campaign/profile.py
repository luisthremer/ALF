#!/usr/bin/env python3
"""Job 5 (P2.9.5/6): profile of the GPU path, Hubbard Mz, L = 24 and 32, delayed N_delay = 32 / 64 (as E2 of CG_ALF_PLAN.md).
Two runs per size under the profiler, NSweep = 2 and NSweep = 4 (NBin 1, Ltau 0): the difference B - A is two sweeps after warm-up
(set-up, mirror build and first-call costs cancel).  usage: profile.py nsys|rocprof EXE WORK OUT [TAG]
Output OUT/profile_<tag>.tsv: per size: wall s per sweep, kernel instances (= launches) per sweep, kernel time s per sweep, GPU busy fraction
(kernel + memcopy time / wall), host sync calls per sweep (cudaStreamSynchronize / cudaDeviceSynchronize / hipStreamSynchronize ...), and
kernel-time shares of the top 10 kernels (B - A)."""
import csv, glob, os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "docs", "ai", "phase2", "p25_runs")); sys.path.insert(0, HERE)
import p25_run as p25, campaign as camp
def col(row, *keys):
    for k in row:
        if all(x in k.lower() for x in keys): return row[k]
def read(path):
    return list(csv.DictReader(open(path))) if path and os.path.exists(path) else []
def find(d, pat):
    r = glob.glob(os.path.join(d, "**", pat), recursive=True); return r[0] if r else None
def profile(tool, exe, d, ed):
    p25.prepare(d, os.path.join(REPO, "Scripts_and_Parameters_files", "Start"), ed)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    for k in [k for k in env if k.startswith("ALF_GPU_")]: del env[k]
    t0 = time.time()
    if tool == "nsys":
        cmd = ["nsys", "profile", "-o", "prof", "--force-overwrite", "true", "-t", "cuda,cublas,cusolver,openmp,osrt", "--sample=none", exe]
    else:
        cmd = ["rocprofv3", "--kernel-trace", "--memory-copy-trace", "--hip-trace", "--stats", "--output-format", "csv", "-d", d + "/out", "-o", "run", "--", exe]
    with open(d + "/stdout", "w") as o: rc = subprocess.call(cmd, cwd=d, env=env, stdout=o, stderr=subprocess.STDOUT)
    wall = time.time() - t0
    if tool == "nsys":
        subprocess.call(["nsys", "stats", "--report", "cuda_gpu_kern_sum,cuda_api_sum,cuda_gpu_mem_time_sum", "--format", "csv", "-o", "stats", "prof.nsys-rep"],
                        cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        K = read(find(d, "stats_cuda_gpu_kern_sum.csv")); A = read(find(d, "stats_cuda_api_sum.csv")); M = read(find(d, "stats_cuda_gpu_mem_time_sum.csv"))
        ker = {col(r, "name") if col(r, "name") else list(r.values())[-1]: (float(col(r, "total time")), int(col(r, "instances"))) for r in K}
        mem = sum(float(col(r, "total time")) for r in M)
        sync = sum(int(col(r, "num calls")) for r in A if "synchronize" in (col(r, "name") or "").lower())
    else:
        K = read(find(d + "/out", "*kernel_stats.csv")); M = read(find(d + "/out", "*memory_copy_stats.csv")); A = read(find(d + "/out", "*hip_api_stats.csv"))
        ker = {col(r, "name"): (float(col(r, "totalduration") or col(r, "total")), int(col(r, "calls"))) for r in K}
        mem = sum(float(col(r, "totalduration") or col(r, "total")) for r in M)
        sync = sum(int(col(r, "calls")) for r in A if "synchronize" in (col(r, "name") or "").lower())
    return rc, wall, ker, mem, sync
def main():
    tool, exe, W, out = sys.argv[1:5]; tag = sys.argv[5] if len(sys.argv) > 5 else tool
    rows = ["tag\tL\tN_delay\tsweeps\twall_s_per_sweep\tlaunches_per_sweep\tkernel_s_per_sweep\tgpu_busy_fraction\tsync_calls_per_sweep\ttop kernels (share of kernel time, 2-sweep difference)"]
    for L, nd in ((24, 32), (32, 64)):
        res = {}
        for n in (2, 4):
            ed = camp.edits("hub-mz", L, "dly16", {"N_delay": str(nd), "NSweep": str(n), "NBin": "1", "Ltau": "0"})
            res[n] = profile(tool, exe, os.path.join(W, "L%d-n%d" % (L, n)), ed)
        (ra, ta, ka, ma, sa), (rb, tb, kb, mb, sb) = res[2], res[4]
        if ra or rb: rows.append("%s\t%d\t%d\t2\tFAILED rc=%d/%d\t\t\t\t\t" % (tag, L, nd, ra, rb)); continue
        unit = 1e-9; dk = {k: (kb[k][0] - ka.get(k, (0, 0))[0], kb[k][1] - ka.get(k, (0, 0))[1]) for k in kb}
        kt = sum(v[0] for v in dk.values()) * unit; kn = sum(v[1] for v in dk.values()); wall = tb - ta; busy = (kt + (mb - ma) * unit) / wall if wall > 0 else float("nan")
        top = sorted(dk.items(), key=lambda kv: -kv[1][0])[:10]
        rows.append("%s\t%d\t%d\t2\t%.3f\t%.0f\t%.3f\t%.2f\t%.0f\t%s" % (tag, L, nd, wall / 2, kn / 2, kt / 2, busy, (sb - sa) / 2,
                    "; ".join("%s %.1f%% (%d)" % (k[:50], 100 * v[0] * unit / kt, v[1]) for k, v in top)))
        print(rows[-1], flush=True)
    os.makedirs(out, exist_ok=True); open(os.path.join(out, "profile_%s.tsv" % tag), "w").write("\n".join(rows) + "\n")
main()
