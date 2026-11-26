#!/usr/bin/env python3
"""One ALF run with edited parameters; prints CPU time, acceptance and the precision lines of `info`.
usage: run1.py OUT ALF.out START_DIR key=value ...   (the environment is passed on unchanged)
A key that is not in START_DIR/parameters is added to &VAR_QMC (e.g. Sequential_scheme, N_delay)."""
import os, re, shutil, subprocess, sys, time
out, exe, start = sys.argv[1:4]
text = open(os.path.join(start, "parameters")).read()
for a in sys.argv[4:]:
    k, v = a.split("=", 1)
    text, n = re.subn(r"(?mi)^(\s*%s\s*=\s*)[^!\n]*" % re.escape(k), lambda m: m.group(1) + v + " ", text, count=1)
    if n != 1:
        text, n = re.subn(r"(?mi)^(\s*&VAR_QMC\b.*\n)", lambda m: m.group(1) + "%s = %s\n" % (k, v), text, count=1)
        if n != 1: sys.exit("cannot place key " + k)
shutil.rmtree(out, ignore_errors=True); os.makedirs(out)
open(os.path.join(out, "parameters"), "w").write(text)
shutil.copy(os.path.join(start, "seeds"), out)
t0 = time.time()
with open(os.path.join(out, "stdout"), "w") as f:
    rc = subprocess.call([os.path.abspath(exe)], cwd=out, stdout=f, stderr=subprocess.STDOUT)
wall = time.time() - t0
lines = open(os.path.join(out, "info"), errors="replace").read().splitlines() if os.path.exists(os.path.join(out, "info")) else []
pick = []
for i, l in enumerate(lines):
    if re.search(r"Precision|Acceptance  |CPU Time|GPU path|GPU observables|GPU tau pass", l):
        pick.append(l.strip() + ((" " + lines[i+1].strip()) if "Precision" in l and i + 1 < len(lines) else ""))
print("%s rc=%d wall=%.1fs | %s" % (out, rc, wall, " | ".join(pick)))
