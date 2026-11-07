#!/usr/bin/env python3
"""GPU benchmark driver (runs on the GPU machine, one GPU process at a time).
usage: bench.py OUT.csv "L=16 scheme=submatrix ndelay=4 nwrap=10" ["..." ...]
Each argument is one configuration (run_cmp.py keys; beta defaults to L).  Reports wall per sweep (from the info-file CPU time minus the
mirror setup), peak GPU memory (nvidia-smi sampling, MiB), Precision Green mean/max, acceptance."""
import sys, os, re, subprocess, threading, time
W = os.environ.get("W", os.path.expanduser("~/alf-gpu/work/fg1"))
out = sys.argv[1]
def sample(stop, res):
    while not stop.is_set():
        r = subprocess.run(["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"], capture_output=True, text=True)
        try: res[0] = max(res[0], int(r.stdout.split()[0]))
        except Exception: pass
        time.sleep(1.0)
rows = []
for i, cfg in enumerate(sys.argv[2:]):
    kv = dict(x.split("=") for x in cfg.split())
    tag = f"bench_{os.getpid()}_{i}"
    nsweep = int(kv.get("nsweep", 2)); kv["nsweep"] = str(nsweep); kv.setdefault("nbin", "1"); kv.setdefault("only", "gpu")
    stop = threading.Event(); res = [0]; th = threading.Thread(target=sample, args=(stop, res)); th.start()
    t0 = time.time()
    r = subprocess.run([sys.executable, os.path.join(W, "run_cmp.py"), tag] + [f"{k}={v}" for k, v in kv.items()], capture_output=True, text=True)
    wall = time.time() - t0
    stop.set(); th.join()
    info = open(f"{W}/runs/{tag}/gpu/info").read() if os.path.exists(f"{W}/runs/{tag}/gpu/info") else ""
    def grab(p, g=1):
        m = re.search(p, info); return m.group(g) if m else "nan"
    setup = grab(r"GPU mirror setup time[^:]*:\s*([0-9.]+)")
    cput = grab(r"CPU Time\s*:\s*([0-9.Ee+-]+)")
    pg = re.search(r"Precision Green\s+Mean, Max :\s*([0-9.Ee+-]+)\s+([0-9.Ee+-]+)", info)
    acc = grab(r"Acceptance\s*:\s*([0-9.Ee+-]+)")
    active = "ACTIVE (mirror" in info
    try: per = (float(cput) - float(setup)) / nsweep
    except Exception: per = float("nan")
    rows.append((cfg, per, res[0], pg.group(1) if pg else "nan", pg.group(2) if pg else "nan", acc, active, round(wall, 1)))
    print(f"{cfg} | s/sweep {per:.2f} | peak GPU MiB {res[0]} | PrecG mean/max {rows[-1][3]}/{rows[-1][4]} | acc {acc} | gpu_active {active} | wall {wall:.0f}s", flush=True)
    if not active: print(r.stdout[-500:], r.stderr[-500:])
    subprocess.run(["rm", "-rf", f"{W}/runs/{tag}/gpu/data.h5"])
with open(out, "w") as f:
    f.write("config,s_per_sweep,peak_gpu_MiB,precG_mean,precG_max,acceptance,gpu_active,wall_s\n")
    for r in rows: f.write(",".join(str(x).replace(",", ";") for x in r) + "\n")
