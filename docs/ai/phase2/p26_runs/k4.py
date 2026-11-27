import sys, os, filecmp
sys.path.insert(0, sys.argv[1] + "/testsuite/GPU.tests")
import gpu_e2e as e
w = sys.argv[2]
for m in ("kondo", "tV", "sp", "hub"):
    base = os.path.join(w, ("%s-std" % m) if m != "hub" else "hub-mz-4x4", "active")
    for k in ("dly4", "dly16", "sub4"):
        d = os.path.join(w, "%s-%s" % (m, k), "active")
        if not os.path.isdir(d): continue
        worst, cf = 0.0, []
        for f in e.files(base):
            a, b = os.path.join(base, f), os.path.join(d, f)
            if f.startswith("info") or f in ("parameters","stdout","stderr") or not os.path.exists(b): continue
            if f.startswith("confout"): cf.append(filecmp.cmp(a, b, shallow=False)); continue
            if filecmp.cmp(a, b, shallow=False): continue
            dd = e.numdiff(a, b)
            worst = max(worst, float("inf") if dd is None else dd[1])
        print("%-6s std vs %-6s max scaled diff %.2e  confout identical %s" % (m, k, worst, all(cf)))
