#!/usr/bin/env python3
"""Set up a PARALLEL_PARAMS (PP) run: one parameter set per U value, several MPI ranks (independent walkers) per set.

Point-defect pi-flux Hubbard model as in the L=24 production run (Beta = L, Dtau = 0.1, N_Phi = L*L/2, N_FL = 2, ham_h0 = 1,
submatrix updates). The point defect exists only on branch gpu-port-unified; on gpu-port ham_h0 is
upstream's static pinning field (docs/gpu/BRANCHES.md). Needs only the Python standard library (no f90nml), so it can run on the compute cluster.

usage: make_PP_Uscan.py OUTDIR --L 28 --U 5.0 5.1 5.2 [--ranks-per-set 16] [--nsweep 50] [--n-delay 128] [--cpu-max 47.0]
                        [--seed0 1000] [--nbin 100000] [--mz F|T]
Creates OUTDIR/parameters (with &VAR_TEMP), OUTDIR/Temp_<i>/parameters (set i has U = U[i]) and OUTDIR/seeds (one integer per rank).
With --jobfile TEMPLATE --var NAME=x --var TIME=... the template (e.g. JobfileSuperMUC-NG-Phase2-GPU-PP.sh) is filled in and written to
OUTDIR/job.sh; Nnodes is set automatically (ranks / 8) unless given. NOTHING is submitted.
Run with the PP build (configure.sh ... PP ...): mpiexec -n <ranks-per-set * number of U values> ALF.out in OUTDIR.
Rank r belongs to set r // ranks_per_set; its output goes to Temp_<set>/ (data.h5, info, confout_*).
"""
import argparse, os

TEMPLATE = """&VAR_ham_name
ham_name = "Hubbard"
/
&VAR_lattice
L1           = {L}
L2           = {L}
Lattice_type = "Square"
Model        = "Hubbard"
/
&VAR_Model_Generic
Checkerboard = .T.
Symm         = .T.
N_SUN        = 2
N_FL         = 2
Phi_X        = 0.d0
Phi_Y        = 0.d0
Bulk         = .T.
N_Phi        = {NPHI}
Dtau         = 0.1d0
Beta         = {L}.d0
Projector    = .F.
Theta        = 0.1d0
/
&VAR_QMC
Nwrap                = 10
NSweep               = {NSWEEP}
NBin                 = {NBIN}
Ltau                 = 0
LOBS_ST              = 0
LOBS_EN              = 0
CPU_MAX              = {CPU_MAX}
Propose_S0           = .F.
Global_moves         = .F.
N_Global             = 1
Global_tau_moves     = .F.
N_Global_tau         = 1
Sequential           = .T.
Nt_sequential_start  = 0
Nt_sequential_end    = -1
Sequential_scheme    = "submatrix"
N_delay              = {NDELAY}
Langevin             = .F.
Delta_t_Langevin_HMC = 0.1
Max_Force            = 5.0
HMC                  = .F.
N_HMC_sweeps         = 1
Leapfrog_steps       = 0
Amplitude            = 1.d0
/
&VAR_errors
n_skip        = 1
N_rebin       = 1
N_Cov         = 0
Extended_Zone = .F.
N_BZ_Zones    = 1
/
&VAR_Hubbard
Mz         = .{MZ}.
Continuous = .F.
ham_T      = 1.d0
ham_chem   = 0.d0
ham_U      = {U}d0
ham_h0     = 1.d0
ham_T2     = 1.d0
ham_U2     = 4.d0
ham_Tperp  = 1.d0
/
"""
TEMP = """&VAR_TEMP
N_exchange_steps = 0
N_Tempering_frequency = 100000
mpi_per_parameter_set = {RPS}
Tempering_calc_det = .F.
/
"""

ap = argparse.ArgumentParser()
ap.add_argument("outdir")
ap.add_argument("--L", type=int, required=True)
ap.add_argument("--U", nargs="+", required=True, help="U values, e.g. 5.0 5.1 5.25 (written as 5.25d0)")
ap.add_argument("--ranks-per-set", type=int, default=16)
ap.add_argument("--nsweep", type=int, default=50)
ap.add_argument("--nbin", type=int, default=100000, help="ignored when --cpu-max > 0 (the run stops at CPU_MAX hours)")
ap.add_argument("--n-delay", type=int, default=128)
ap.add_argument("--cpu-max", type=float, default=47.0, help="hours; the run writes its last bin and stops (0 = use NBin)")
ap.add_argument("--seed0", type=int, default=1000)
ap.add_argument("--mz", default="F", choices=["F", "T"])
ap.add_argument("--ranks-per-node", type=int, default=8, help="8 = one rank per GPU tile; 16/32 = 2/4 ranks per tile")
ap.add_argument("--jobfile", help="jobscript template with ##VAR## placeholders to fill in")
ap.add_argument("--var", action="append", default=[], help="NAME=value for a ##NAME## placeholder (repeatable)")
a = ap.parse_args()

def text(u):
    return TEMPLATE.format(L=a.L, NPHI=a.L * a.L // 2, NSWEEP=a.nsweep, NBIN=a.nbin, CPU_MAX=a.cpu_max, NDELAY=a.n_delay,
                           MZ=a.mz, U=u) + "! slash terminates namelist statement - DO NOT REMOVE\n"

os.makedirs(a.outdir, exist_ok=True)
for i, u in enumerate(a.U):
    d = os.path.join(a.outdir, "Temp_%d" % i)
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "parameters"), "w").write(text(u))
open(os.path.join(a.outdir, "parameters"), "w").write(text(a.U[0]).replace("! slash", TEMP.format(RPS=a.ranks_per_set) + "! slash"))
nranks = a.ranks_per_set * len(a.U)
open(os.path.join(a.outdir, "seeds"), "w").write("".join("%d\n" % (a.seed0 + r) for r in range(nranks)))
print("L=%d, %d parameter sets (U = %s), %d ranks per set, %d ranks in total" % (a.L, len(a.U), " ".join(a.U), a.ranks_per_set, nranks))
if a.jobfile:
    kv = dict(v.split("=", 1) for v in a.var)
    kv.setdefault("NtaskPnode", str(a.ranks_per_node))
    kv.setdefault("Nthreads", str(112 // a.ranks_per_node))
    kv.setdefault("Nnodes", str(-(-nranks // a.ranks_per_node)))
    kv.setdefault("CPUMAX", str(a.cpu_max))
    out = []
    for l in open(a.jobfile).read().split("\n"):
        if ("mail-user" in l or "mail-type" in l) and "EMAIL" not in kv:
            continue
        for k, v in kv.items():
            l = l.replace("##%s##" % k, v)
        out.append(l)
    t = "\n".join(out)
    import re
    left = re.findall(r"^[^\n]*##\w+##[^\n]*$", t, re.M)
    open(os.path.join(a.outdir, "job.sh"), "w").write(t)
    print("jobscript written to %s/job.sh (%d nodes); unresolved placeholders: %s" % (a.outdir, int(kv["Nnodes"]), [x.strip()[:60] for x in left if not x.lstrip().startswith("#") or "SBATCH" in x]))
print("set i -> %s/Temp_<i>/ ; run: mpiexec -n %d ALF.out in %s" % (a.outdir, nranks, a.outdir))
