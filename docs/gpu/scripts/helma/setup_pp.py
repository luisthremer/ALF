#!/usr/bin/env python3
"""Create one PARALLEL_PARAMS run directory (stdlib only; replaces pp_setup.py/prepare_sim.py/f90nml).

usage: setup_pp.py DIR L NSWEEP NDELAY CPU_MAX WALKERS U1,U2,...   [--nbin N] [--master TAG]
DIR/parameters   top level file (VAR_TEMP with mpi_per_parameter_set = WALKERS; the rest is read from Temp_i/parameters)
DIR/Temp_i/parameters  one per U (i = position in the list), derived from parameters_ref_L24_U5.3_h1_4w (Jakob's file, /jakob_params)
DIR/seeds        one integer per MPI rank, rank = i*WALKERS + walker; seed = f(master, L, U, walker) so it is unique and
                 does not depend on how the U list is split over jobs.
Only L1 = L2, Beta = L, N_Phi = L*L/2, NSweep (10), NBin, CPU_MAX, N_delay, ham_U differ from the reference.
Existing DIR/Temp_i are never overwritten (restart safety): the script refuses if DIR exists."""
import sys, os, re, hashlib
here = os.path.dirname(os.path.abspath(__file__))
a = [x for x in sys.argv[1:] if not x.startswith("--")]
opt = dict(zip([x for x in sys.argv[1:] if x.startswith("--")][::1], [None]*9))
nbin, master = 100000, "helma2026"
args = sys.argv[1:]
for k in ("--nbin", "--master"):
    if k in args:
        i = args.index(k); v = args[i+1]; del args[i:i+2]
        if k == "--nbin": nbin = int(v)
        else: master = v
d, L, nsweep, nd, cpumax, nw, us = args[0], int(args[1]), int(args[2]), int(args[3]), args[4], int(args[5]), args[6].split(",")
if os.path.exists(d): sys.exit("refuse: %s exists" % d)
ref = open(os.path.join(here, "parameters_ref_L24_U5.3_h1_4w")).read()
def sub(t, key, val):
    t2, n = re.subn(r"(?m)^(%s\s*=\s*)[^\n]*" % key, lambda m: m.group(1) + val, t)
    assert n == 1, key
    return t2
def fortran_d(u):
    s = u if "." in u else u + ".0"
    return s + "d0" if "d" not in s else s
t = ref
t = sub(t, "L1", str(L)); t = sub(t, "L2", str(L)); t = sub(t, "Beta", "%d.d0" % L); t = sub(t, "N_Phi", str(L*L//2))
t = sub(t, "NSweep", str(nsweep)); t = sub(t, "NBin", str(nbin)); t = sub(t, "CPU_MAX", cpumax); t = sub(t, "N_delay", str(nd))
os.makedirs(d)
vt = "&VAR_TEMP\nN_exchange_steps = 0\nN_Tempering_frequency = 100000\nmpi_per_parameter_set = %d\nTempering_calc_det = .F.\n/\n" % nw
open(os.path.join(d, "parameters"), "w").write(vt + sub(t, "ham_U", fortran_d(us[0])))
seeds, used = [], set()
for i, u in enumerate(us):
    os.makedirs(os.path.join(d, "Temp_%d" % i))
    open(os.path.join(d, "Temp_%d" % i, "parameters"), "w").write(sub(t, "ham_U", fortran_d(u)))
    for w in range(nw):
        s = int(hashlib.sha256(("%s|L%d|U%s|w%d" % (master, L, u, w)).encode()).hexdigest()[:8], 16) % (2**30) + 1
        assert s not in used; used.add(s); seeds.append(s)
open(os.path.join(d, "seeds"), "w").write("\n".join("%d" % s for s in seeds) + "\n")
print("%s: %d groups x %d walkers = %d ranks, L=%d NSweep=%d N_delay=%d CPU_MAX=%s" % (d, len(us), nw, len(us)*nw, L, nsweep, nd, cpumax))
