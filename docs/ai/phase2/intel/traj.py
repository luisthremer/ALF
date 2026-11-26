#!/usr/bin/env python3
"""Trajectory check GPU=intel vs GPU=none (same ifx, same parameters and seeds), plus the S1 check mode.
usage: traj.py OUT ALF-intel.out ALF-none.out START_DIR TOL key=value ...
Runs OUT/none, OUT/gpu (default) and OUT/gpu_lu2 (ALF_GPU_DEBUG=1 ALF_GPU_CGR_LU=2) and compares gpu and gpu_lu2
with none (compare.py, tolerance TOL); prints the GPU info lines, Precision Green and the check-mode phase line."""
import os, re, shutil, subprocess, sys
out, gpu, none, start, tol = sys.argv[1:6]
pars = dict(a.split("=", 1) for a in sys.argv[6:])
here = os.path.dirname(os.path.abspath(__file__))
text = open(os.path.join(start, "parameters")).read()
for k, v in pars.items():   # first occurrence of the key only (the generic namelists come first)
    text, n = re.subn(r"(?mi)^(\s*%s\s*=\s*)[^!\n]*" % re.escape(k), lambda m: m.group(1) + v + " ", text, count=1)
    if n != 1: sys.exit("unknown key " + k)
runs = {"none": (none, {}), "gpu": (gpu, {}), "gpu_lu2": (gpu, {"ALF_GPU_DEBUG": "1", "ALF_GPU_CGR_LU": "2"})}
for lab, (exe, env) in runs.items():
    d = os.path.join(out, lab)
    shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    open(os.path.join(d, "parameters"), "w").write(text)
    shutil.copy(os.path.join(start, "seeds"), d)
    e = dict(os.environ, **env)
    if lab == "none": e.pop("OMP_TARGET_OFFLOAD", None)
    with open(os.path.join(d, "stdout"), "w") as f:
        rc = subprocess.call([os.path.abspath(exe)], cwd=d, env=e, stdout=f, stderr=subprocess.STDOUT)
    print("%s rc=%d" % (lab, rc))
for lab in ("gpu", "gpu_lu2"):
    for l in open(os.path.join(out, lab, "info"), errors="replace"):
        if re.search(r"GPU path|Precision Green|phase diff|CGR calls|oneMKL", l): print(lab, "|", l.rstrip())
    r = subprocess.run([sys.executable, os.path.join(here, "compare.py"), os.path.join(out, "none"), os.path.join(out, lab),
                        "--tol", tol], capture_output=True, text=True)
    print(lab, "vs none:", r.stdout.strip().splitlines()[-1])
for l in open(os.path.join(out, "none", "info"), errors="replace"):
    if "Precision Green" in l: print("none |", l.rstrip())
