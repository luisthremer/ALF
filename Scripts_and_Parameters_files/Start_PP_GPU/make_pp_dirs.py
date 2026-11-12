#!/usr/bin/env python3
"""Set up a PARALLEL_PARAMS (PP) run directory for the point-defect Hubbard model: one parameter set per U, several MPI walkers per set.
Needs only the python standard library (no f90nml). Layout as written by prepare_sim.py:
    DIR/parameters, DIR/seeds, DIR/Temp_<i>/parameters  (i = 0..nU-1, ham_U = U[i]);  ranks = nU * mpi_per_parameter_set

usage: make_pp_dirs.py DIR L U1,U2,... MPI_PER_SET NSWEEP NBIN [N_DELAY=32] [CPU_MAX_HOURS=0.0] [SEED_OFFSET=0] [TEMPLATE]
Beta = L (beta = L, as in the reference run), N_Phi = L^2/2. Seeds: one random integer per rank (distinct), generated from SEED_OFFSET
so that different runs can be given different, reproducible seed streams.
Restart: put the confout_<rank> files of Temp_<i>/ into confin_<rank> (out_to_in_pp.sh) and run again; data.h5 bins are appended
in a new directory per segment (see JobfileSuperMUC-NG-Phase2-GPU-PP.sh)."""
import os, random, re, sys
d, L, Us, nper, nsweep, nbin = sys.argv[1], int(sys.argv[2]), sys.argv[3].split(','), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
nd = int(sys.argv[7]) if len(sys.argv) > 7 else 32
cpu_max = sys.argv[8] if len(sys.argv) > 8 else "0.0"
offset = int(sys.argv[9]) if len(sys.argv) > 9 else 0
tmpl = sys.argv[10] if len(sys.argv) > 10 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "parameters_template")
t = open(tmpl).read()
def sub(key, val):
    global t
    t, n = re.subn(r"(?m)^(\s*%s\s*=).*$" % key, r"\g<1> %s" % val, t); assert n == 1, key
for k in ("L1", "L2"): sub(k, L)
sub("N_Phi", L * L // 2); sub("Beta", "%d.d0" % L); sub("NSweep", nsweep); sub("NBin", nbin); sub("N_delay", nd); sub("CPU_MAX", cpu_max)
tp = t.replace("&VAR_errors", "&VAR_TEMP\nN_exchange_steps = 0\nN_Tempering_frequency = 100000\nmpi_per_parameter_set = %d\nTempering_calc_det = .F.\n/\n&VAR_errors" % nper)
os.makedirs(d, exist_ok=True)
t = tp; sub("ham_U", "%sd0" % Us[0] if "d" not in Us[0] else Us[0])
open(os.path.join(d, "parameters"), "w").write(t)
for i, u in enumerate(Us):
    os.makedirs(os.path.join(d, "Temp_%d" % i), exist_ok=True)
    sub("ham_U", "%sd0" % u)
    open(os.path.join(d, "Temp_%d" % i, "parameters"), "w").write(t)
nr = nper * len(Us)
rng = random.Random(1000003 * (offset + 1) + L)
seeds = rng.sample(range(1, 999999), nr)
open(os.path.join(d, "seeds"), "w").write("".join("%12d\n" % s for s in seeds))
print("wrote %s: L=%d, %d sets (U=%s) x %d ranks = %d ranks, NSweep=%d NBin=%d N_delay=%d" % (d, L, len(Us), ",".join(Us), nper, nr, nsweep, nbin, nd))
