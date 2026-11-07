#!/usr/bin/env python3
"""Fixed-seed GPU-vs-CPU comparison driver (runs on the GPU machine).
usage: run_cmp.py TAG [key=value ...]   keys: L mz scheme ndelay ltau nsweep nbin U h0 symm nwrap checker nphi beta
Creates $W/runs/TAG/{cpu,gpu}/ (fresh, empty), runs ALF-cpu-ref.out and ALF-gpu.out with the same parameters and seeds,
then compare.py.  Env: W (default ~/alf-gpu/work/fg1), CPUBIN, GPUBIN."""
import sys, os, shutil, subprocess, time
W = os.environ.get("W", os.path.expanduser("~/alf-gpu/work/fg1"))
cpubin = os.environ.get("CPUBIN", W + "/bin/ALF-cpu-ref.out")
gpubin = os.environ.get("GPUBIN", W + "/bin/ALF-gpu.out")
tag = sys.argv[1]
o = dict(L=4, mz="T", scheme="standard", ndelay=4, ltau=0, nsweep=10, nbin=2, U="5.0d0", h0="1.d0", symm="T", nwrap=10,
         checker="F", nsun=2, nfl=2, beta=None, nphi=None, dtau="0.1d0", only=None)
for kv in sys.argv[2:]:
    k, v = kv.split("="); o[k] = v
L = int(o["L"]); beta = o["beta"] or f"{L}.d0"; nphi = o["nphi"] or str(L*L//2)
sched = "" if o["scheme"] == "standard" else f'Sequential_scheme = "{o["scheme"]}"\nN_delay = {o["ndelay"]}\n'
P = f"""&VAR_ham_name
ham_name = "Hubbard"
/
&VAR_lattice
L1 = {L}
L2 = {L}
Lattice_type = "Square"
Model = "Hubbard"
/
&VAR_Model_Generic
Checkerboard = .{o['checker']}.
Symm = .{o['symm']}.
N_SUN = {o['nsun']}
N_FL = {o['nfl']}
Phi_X = 0.d0
Phi_Y = 0.d0
Bulk = .T.
N_Phi = {nphi}
Dtau = {o['dtau']}
Beta = {beta}
Projector = .F.
Theta = 10.d0
/
&VAR_QMC
Nwrap = {o['nwrap']}
NSweep = {o['nsweep']}
NBin = {o['nbin']}
Ltau = {o['ltau']}
LOBS_ST = 0
LOBS_EN = 0
CPU_MAX = 0.0
Propose_S0 = .F.
Global_moves = .F.
N_Global = 1
Global_tau_moves = .F.
N_Global_tau = 1
Sequential = .T.
Nt_sequential_start = 0
Nt_sequential_end = -1
Langevin = .F.
HMC = .F.
{sched}/
&VAR_errors
n_skip = 1
N_rebin = 1
N_Cov = 0
Extended_Zone = .F.
N_BZ_Zones = 1
/
&VAR_Hubbard
Mz = .{o['mz']}.
Continuous = .F.
ham_T = 1.d0
ham_chem = 0.d0
ham_U = {o['U']}
ham_h0 = {o['h0']}
/
! slash terminates namelist statement - DO NOT REMOVE
"""
base = os.path.join(W, "runs", tag)
if os.path.exists(base): shutil.rmtree(base)
times = {}
which = [o["only"]] if o["only"] else ["cpu", "gpu"]
for kind in which:
    d = os.path.join(base, kind); os.makedirs(d)
    open(os.path.join(d, "parameters"), "w").write(P)
    shutil.copy(W + "/src_gpu/Scripts_and_Parameters_files/Start/seeds", os.path.join(d, "seeds"))
    t0 = time.time()
    env = dict(os.environ, TMPDIR=W + "/tmp")
    with open(os.path.join(d, "run.log"), "w") as lg:
        r = subprocess.run([cpubin if kind == "cpu" else gpubin], cwd=d, stdout=lg, stderr=subprocess.STDOUT, env=env)
    times[kind] = time.time() - t0
    print(f"{tag} {kind}: rc={r.returncode} wall={times[kind]:.1f}s", flush=True)
    if r.returncode: print(open(os.path.join(d, "run.log")).read()[-1500:])
if len(which) == 2:
    r = subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), "compare.py"),
                        base + "/cpu", base + "/gpu"], capture_output=True, text=True)
    print(r.stdout[-3000:], r.stderr[-500:])
