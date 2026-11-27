#!/usr/bin/env python3
"""Device-path gate (manager, 2026-11-26, after the AMD offload finding): one binary, same seeds, GPU path vs
ALF_GPU_DISABLE=1 (CPU path), plus the check-mode lane ALF_GPU_DEBUG=1 ALF_GPU_CGR_LU=2.
usage: gate.py OUT ALF.out START [env:NAME=VALUE ...] key=value ...   (env: only for the GPU lanes, e.g. a fault)
PASS when: all three runs end with rc 0, the GPU path is ACTIVE, no NaN/Inf in any output file, and every number of every
observable file (scalar, equal-time, local, time-displaced) and of confout_0 agrees with the CPU path within
tol = 100 x the CPU path's Precision Green mean (the rule of testsuite/GPU.tests/gpu_e2e.py), relative to the file's
largest entry; the check mode reports a phase difference below 1e-8."""
import os, re, glob, shutil, subprocess, sys
out, exe, start = sys.argv[1:4]
text = open(os.path.join(start, "parameters")).read()
genv = dict(a[4:].split("=", 1) for a in sys.argv[4:] if a.startswith("env:"))
for a in [a for a in sys.argv[4:] if not a.startswith("env:")]:
    k, v = a.split("=", 1)
    text, n = re.subn(r"(?mi)^(\s*%s\s*=\s*)[^!\n]*" % re.escape(k), lambda m: m.group(1) + v + " ", text, count=1)
    if n != 1:
        text, n = re.subn(r"(?mi)^(\s*&VAR_QMC\b.*\n)", lambda m: m.group(1) + "%s = %s\n" % (k, v), text, count=1)
        if n != 1: sys.exit("cannot place key " + k)
lanes = {"cpu": {"ALF_GPU_DISABLE": "1"}, "gpu": dict(genv), "lu2": dict(genv, ALF_GPU_DEBUG="1", ALF_GPU_CGR_LU="2")}
rc = {}
for lab, env in lanes.items():
    d = os.path.join(out, lab); shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    open(os.path.join(d, "parameters"), "w").write(text); shutil.copy(os.path.join(start, "seeds"), d)
    with open(os.path.join(d, "stdout"), "w") as f:
        rc[lab] = subprocess.call([os.path.abspath(exe)], cwd=d, env=dict(os.environ, **env), stdout=f, stderr=subprocess.STDOUT)
num = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?|[-+]?nan|[-+]?inf(?:inity)?", re.I)
def nums(p): return [float(x.replace("D", "E").replace("d", "e")) for x in num.findall(open(p, errors="replace").read())]
def info(lab): return open(os.path.join(out, lab, "info"), errors="replace").read() if os.path.exists(os.path.join(out, lab, "info")) else ""
skip = ("info", "parameters", "seeds", "stdout")
files = sorted(os.path.basename(p) for p in glob.glob(os.path.join(out, "cpu", "*")) if os.path.basename(p) not in skip and os.path.isfile(p))
m = re.search(r"Precision Green\s+Mean, Max :\s*([-+0-9.Ee]+)", info("cpu"))
pg = float(m.group(1)) if m else float("nan")
tol = 100 * pg
bad, worst, nan = [], {}, []
for lab in ("gpu", "lu2"):
    w = 0.0
    for f in files:
        a, b = os.path.join(out, "cpu", f), os.path.join(out, lab, f)
        if not os.path.exists(b): bad.append("%s missing %s" % (lab, f)); continue
        x, y = nums(a), nums(b)
        if any(v != v or abs(v) == float("inf") for v in y): nan.append("%s:%s" % (lab, f))
        if len(x) != len(y): bad.append("%s shape %s" % (lab, f)); continue
        s = max([abs(u) for u in x], default=0) or 1.0
        w = max(w, max([abs(u - v) / s for u, v in zip(x, y)], default=0))
    worst[lab] = w
active = "GPU path (gpu_state): ACTIVE" in info("gpu")
mp = re.search(r"max \|phase\(tracked\) - phase\(LU\)\|\s*([-+0-9.Ee]+)", info("lu2"))
phd = float(mp.group(1)) if mp else float("nan")
ok = all(v == 0 for v in rc.values()) and active and not nan and not bad and tol == tol and \
     worst["gpu"] <= tol and worst["lu2"] <= tol and phd < 1e-8
print("%s %s rc=%s active=%s tol=%.2e worst gpu=%.2e lu2=%.2e phase=%.1e %s%s" % (
    "PASS" if ok else "FAIL", out, ",".join("%s:%d" % kv for kv in rc.items()), active, tol, worst["gpu"], worst["lu2"],
    phd, ("NaN in " + " ".join(nan) + " ") if nan else "", " ".join(bad)))
