#!/usr/bin/env python3
"""usage: bench_summarize.py RUNDIR NGPU RPG PCT ND L NSWEEPS  -> one CSV line
s/sweep per chain = (FIN TIME - START TIME)/NSWEEPS from Temp_i/info (rank 0 of each group; START is printed after the Hamiltonian setup,
the mirror setup (~0.3 s) is before it); sweeps/s per GPU = sum over the GPU's chains of 1/(s/sweep); GPU memory = max fb over dmon lines; sm = mean sm% over lines with sm > 0."""
import sys, re, os, glob, datetime as dt
d, ngpu, rpg, pct, nd, L, ns = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6], int(sys.argv[7])
def ts(s): return dt.datetime.fromisoformat(re.search(r"(\d{4}-\d\d-\d\dT[\d:.]+)", s).group(1))
t, pg, acc = [], [], []
for i in range(ngpu):
    try:
        info = open(os.path.join(d, "Temp_%d" % i, "info")).read()
        a = ts(re.search(r"START TIME:.*", info).group(0)); b = ts(re.search(r"FIN TIME:.*", info).group(0))
        t.append((b - a).total_seconds() / ns)
        pg.append(float(re.search(r"Precision Green\s+Mean, Max\s*:.*\n\s*(\S+)", info).group(1).replace("D", "E")))
    except Exception as e:
        t.append(float("nan"))
sm, fb = [], 0
try:
    for l in open(os.path.join(d, "dmon.log")):
        p = l.split()
        if len(p) >= 8 and p[1].isdigit():   # time gpu sm mem enc dec fb bar1 ... (-s um -o T)
            if int(p[2]) > 0: sm.append(int(p[2]))
            fb = max(fb, int(p[-2]) if p[-2].isdigit() else 0)
except Exception: pass
ok = [x for x in t if x == x]
mean = sum(ok) / len(ok) if ok else float("nan")
thr = sum(1.0 / x for x in ok) / max(1, len(ok)) * rpg if ok else float("nan")   # per GPU
wall = ""
try: wall = "%.1f" % (float(open(os.path.join(d, "t_end")).read()) - float(open(os.path.join(d, "t_start")).read()))
except Exception: pass
print("L=%s,ngpu=%d,ranks_per_gpu=%d,mps=%s,ndelay=%s,sweeps=%d,s_per_sweep_chain=%.2f,sweeps_per_s_per_gpu=%.3f,fb_max_MiB=%d,sm_mean_pct=%.0f,PrecGreenMax=%s,wall_s=%s" %
      (L, ngpu, rpg, pct, nd, ns, mean, thr, fb, (sum(sm) / len(sm) if sm else 0), "%.1e" % max(pg) if pg else "nan", wall))
