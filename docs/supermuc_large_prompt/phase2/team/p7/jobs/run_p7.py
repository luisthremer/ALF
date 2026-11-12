#!/usr/bin/env python3
"""P7 copy of run_ref_sng.py (team/jobs): runs the L=24 reference parameters for a SHORT time.  Additions: env SCHEME (Sequential_scheme) and NDELAY (N_delay)
override the parameters file; env LSIZE as before.  usage: run_p7.py TAG BIN NR NSWEEP NBIN [wrap-script]   env: W, SEEDS, PARAMS, OMPT"""
import os, re, sys, shutil, subprocess, time
tag, binary, nr, nsweep, nbin = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
wrap = [sys.argv[6]] if len(sys.argv) > 6 else []
d = os.path.join(os.environ["W"], "runs", tag)
if os.path.exists(d): shutil.rmtree(d)
os.makedirs(d)
P = open(os.environ["PARAMS"]).read()
P = re.sub(r"NSweep\s*=\s*\d+", f"NSweep               = {nsweep}", P)
P = re.sub(r"NBin\s*=\s*\d+", f"NBin                 = {nbin}", P)
P = re.sub(r"CPU_MAX\s*=\s*[\d.]+", "CPU_MAX              = 0.0", P)
if os.environ.get("SCHEME"): P = re.sub(r'Sequential_scheme\s*=\s*"[^"]*"', f'Sequential_scheme    = "{os.environ["SCHEME"]}"', P)
if os.environ.get("NDELAY"): P = re.sub(r"N_delay\s*=\s*\d+", f"N_delay              = {os.environ['NDELAY']}", P)
open(os.path.join(d, "parameters"), "w").write(P)
shutil.copy(os.environ["SEEDS"], os.path.join(d, "seeds"))
t0 = time.time()
with open(os.path.join(d, "run.log"), "w") as lg:
    r = subprocess.run(["mpiexec", "-n", str(nr)] + wrap + [binary], cwd=d, stdout=lg, stderr=subprocess.STDOUT,
                       env=dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMPT", "1")))
w = time.time() - t0
print(f"{tag}: rc={r.returncode} wall={w:.1f}s  ({nsweep} sweeps/rank, {nr} rank(s))", flush=True)
