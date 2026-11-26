#!/usr/bin/env python3
"""T5 driver: run ALF (GPU=host build) in check mode ALF_GPU_DEBUG=1 ALF_GPU_CGR_LU=2 and report the CGR phase check line.
usage: t5run.py <build> <case> [KEY=VAL ...] [--np N] [--env K=V] [--dirs Temp_0,Temp_1] [--restart] [--omp N]"""
import os, re, shutil, subprocess, sys, time
T5 = os.path.expanduser("~/work/ALF_P2/t5")
args = sys.argv[1:]
build, case = args[0], args[1]
edits, env_extra, np_, dirs, restart, omp, blas = {}, {}, 0, [""], False, "1", None
i = 2
while i < len(args):
    a = args[i]
    if a == "--np": np_ = int(args[i+1]); i += 2; continue
    if a == "--env": k, v = args[i+1].split("=", 1); env_extra[k] = v; i += 2; continue
    if a == "--dirs": dirs = args[i+1].split(","); i += 2; continue
    if a == "--restart": restart = True; i += 1; continue
    if a == "--omp": omp = args[i+1]; i += 2; continue
    k, v = a.split("=", 1); v = v.strip("\"")
    if k in ("Sequential_scheme", "Lattice_type", "ham_name", "Model"): v = "\"%s\"" % v
    edits[k] = v; i += 1
src = os.path.join(T5, "src_" + build)
exe = os.path.join(src, "Prog", "ALF.out")
start = os.path.join(src, "Scripts_and_Parameters_files", "Start")
rd = os.path.join(T5, "runs", build, case)
def mkparams(d):
    os.makedirs(d, exist_ok=True)
    shutil.copy(os.path.join(start, "seeds"), d)
    t = open(os.path.join(start, "parameters")).read()
    for k, v in edits.items():
        t, n = re.subn(r"(?m)^(\s*%s\s*=\s*)[^\s!]+" % re.escape(k), lambda m: m.group(1) + v, t)
        if n == 0:
            t, n = re.subn(r"(?m)^(Ltau\s*=)", "%s = %s\n\\g<1>" % (k, v), t, count=1)
        if n == 0: sys.exit("cannot set " + k)
    open(os.path.join(d, "parameters"), "w").write(t)
def run(tag):
    env = dict(os.environ, OMP_NUM_THREADS=omp, ALF_GPU_DEBUG="1", ALF_GPU_CGR_LU="2")
    env.update(env_extra)
    cmd = ["mpirun", "--oversubscribe", "-np", str(np_), exe] if np_ else [exe]
    t0 = time.time()
    with open(os.path.join(rd, "stdout" + tag), "w") as o:
        rc = subprocess.call(cmd, cwd=rd, env=env, stdout=o, stderr=subprocess.STDOUT)
    return rc, time.time() - t0
def report(tag, rc, dt):
    for d in dirs:
        inf = os.path.join(rd, d, "info")
        txt = open(inf, errors="replace").read() if os.path.exists(inf) else ""
        act = re.search(r"GPU path \(gpu_state\): (.*)", txt)
        chk = re.findall(r"CGR calls compared (\d+), max \|phase\(tracked\) - phase\(LU\)\| +(\S+)", txt)
        ph = re.search(r"Precision Phase, Max\s*:\s*(\S+)", txt)
        print("%-6s %-26s %-8s rc=%d %6.0fs %s | check %s | PrecPhase %s | %s" % (build, case + tag, d or ".", rc, dt,
              (act.group(1)[:60] if act else "no GPU line"), chk if chk else "none", ph.groups() if ph else "-",
              "abort: " + re.sub(r"\s+", " ", open(os.path.join(rd, "stdout" + tag), errors="replace").read()[-300:]) if rc else ""))
if os.path.exists(rd): shutil.rmtree(rd)
os.makedirs(rd)
mkparams(rd)
if dirs != [""]:
    for d in dirs:
        mkparams(os.path.join(rd, d))
rc, dt = run(""); report("", rc, dt)
if restart and rc == 0:
    for d in dirs:
        dd = os.path.join(rd, d)
        for f in os.listdir(dd):
            if f.startswith("confout_"):
                shutil.copy(os.path.join(dd, f), os.path.join(dd, f.replace("confout", "confin")))
    rc, dt = run("_restart"); report("_restart", rc, dt)
