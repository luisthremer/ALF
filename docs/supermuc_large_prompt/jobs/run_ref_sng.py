#!/usr/bin/env python3
"""Run the L=24 reference parameter set (docs/supermuc_large_prompt/ref_L24/parameters: L=24, Beta=24, U=5.3, h0=1, N_FL=2, Mz=F,
submatrix N_delay=32, Nwrap=10, Dtau=0.1) for a SHORT time on SuperMUC-NG Phase 2.
usage: run_ref_sng.py TAG BIN NR NSWEEP NBIN [wrap-script]   env: W, SEEDS, PARAMS; OMP_TARGET_OFFLOAD from the caller"""
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
open(os.path.join(d, "parameters"), "w").write(P)
shutil.copy(os.environ["SEEDS"], os.path.join(d, "seeds"))
t0 = time.time()
with open(os.path.join(d, "run.log"), "w") as lg:
    r = subprocess.run(["mpiexec", "-n", str(nr)] + wrap + [binary], cwd=d, stdout=lg, stderr=subprocess.STDOUT,
                       env=dict(os.environ, OMP_NUM_THREADS="1"))
w = time.time() - t0
sweeps = int(nsweep) * int(nbin)
print(f"{tag}: rc={r.returncode} wall={w:.1f}s  ({sweeps} sweeps/rank, {nr} rank(s); wall/sweep incl. setup {w/sweeps:.2f}s)", flush=True)
