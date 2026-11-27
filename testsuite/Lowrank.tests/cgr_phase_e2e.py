#!/usr/bin/env python3
"""CGR determinant phase from the tracked det(U) against the LU (cgr1_mod, ALF_CGR_LU), whole ALF runs on fixed seeds.

Each case runs with ALF_CGR_LU = 0 (tracked det(U)), 1 (LU) and 2 (both, compared at every CGR), so the checks do not
depend on the default of the build:
  mode 2: stdout has "CGR phase check (ALF_CGR_LU=2): CGR calls compared N, max ... X" on every rank with N > 0 and
          X < 1e-8 (the run stops above), and every output file is byte-identical to mode 0 (the tracked phase is used);
  mode 0 against 1: confout byte-identical (same Markov chain), every other output number equal to TOL (the two phases
          differ at rounding level).
Negative controls: ALF_CGR_FAULT=det (sign of the tracked phase flipped) in mode 2 must stop the run; ALF_CGR_LU=7 must
be refused.  --no-tracked (STAB1/STAB2 builds, whose CGR keeps its own LU): mode 2 must report 0 comparisons, faults are
not checked.

Usage: cgr_phase_e2e.py --exe ALF.out --start Start/ --work DIR [--build serial|mpi|tempering|parallel_params]
                        [--mpiexec "mpirun ..."] [--no-tracked]
"""
import argparse, filecmp, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from compare_schemes import make_text, prepare, with_scheme, numbers, h5_text  # noqa: E402

TOL = 1.0e-10
TOL_DEV = 1.0e-13   # mode 2 against mode 0 on a device backend: two runs of the same binary differ at rounding level
CHECK = re.compile(r"CGR phase check \(ALF_CGR_LU=2\): CGR calls compared (\d+), max \|phase\(tracked\) - phase\(LU\)\| +(\S+)")
SKIP = {"info", "info_git", "parameters", "seeds", "stdout"}
BASE = {"L1": 4, "L2": 4, "Beta": 2.0, "Nsweep": 4, "NBin": 1, "Ltau": 1,
        "N_Tempering_frequency": 2, "N_exchange_steps": 1, "Tempering_calc_det": False}
# tempering builds: one exchange step without det moves the states between ranks (UDV_State MPI_Sendrecv, Exchange_Step)
CASES = [  # complex phases (flux, N_FL = 1), real sign problem with delayed updates, non-diagonal operators
    ("su2-flux", dict(BASE, Mz=False, N_Phi=1)),
    ("mz-chem-delayed", dict(BASE, Ham_chem=0.3, scheme=("delayed", 4))),
    ("kondo", dict(BASE, Model="Kondo", Lattice_type="Bilayer_square", Ham_Uf=1.0)),
]


def run(exe, d, env, a):
    ranks = 2 if a.build != "serial" else 1
    cmd = (a.mpiexec.split() + ["-np", str(ranks)] if a.mpiexec else []) + [exe]
    e = {k: v for k, v in os.environ.items() if not k.startswith("ALF_CGR_")}
    e.update(env, OMP_NUM_THREADS="1")
    with open(os.path.join(d, "stdout"), "w") as out:
        return subprocess.call(cmd, cwd=d, env=e, stdout=out, stderr=subprocess.STDOUT)


def outputs(d):
    res = []
    for root, _, names in os.walk(d):
        res += [os.path.relpath(os.path.join(root, n), d) for n in names if n not in SKIP]
    return sorted(res)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--start", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--build", default="serial")
    ap.add_argument("--mpiexec", default="")
    ap.add_argument("--no-tracked", action="store_true")
    ap.add_argument("--test-models", action="store_true")  # passed by CMake with the lowrank options; unused
    a = ap.parse_args()
    template = open(os.path.join(a.start, "parameters")).read()
    seeds = os.path.join(a.start, "seeds")
    nrep = 2 if a.build in ("tempering", "parallel_params") else 1
    bad = []

    def fresh(label, pars, mode_env):
        d = os.path.join(a.work, label)
        pars = dict(pars)
        scheme = pars.pop("scheme", None)
        text, missing = make_text(template, pars)
        if missing:
            sys.exit("cannot set %s" % missing)
        if scheme:
            text = with_scheme(text, *scheme)
        prepare(d, [text] * nrep, seeds, 1, a.build)
        return d, run(a.exe, d, mode_env, a)

    for name, pars in CASES:
        dirs = {}
        for m in ("0", "1", "2"):
            dirs[m], rc = fresh("%s/mode%s" % (name, m), pars, {"ALF_CGR_LU": m})
            if rc:
                bad.append("%s mode %s: exit %d" % (name, m, rc))
        if any(b.startswith(name + " ") for b in bad):
            continue
        found = CHECK.findall(open(os.path.join(dirs["2"], "stdout"), errors="replace").read())
        ranks = 2 if a.build != "serial" else 1
        if len(found) != ranks:
            bad.append("%s: %d check lines, expected %d" % (name, len(found), ranks))
        for n, x in found:
            ok = int(n) == 0 if a.no_tracked else (int(n) > 0 and float(x) < 1e-8)
            if not ok:
                bad.append("%s: check line compared %s, max %s" % (name, n, x))
        print("%-16s mode 2: %s" % (name, ", ".join("%s CGRs max %s" % f for f in found)))
        info0 = os.path.join(dirs["0"], "info")                 # missing in tempering builds (Temp_<n>/info): no GPU path there
        m_dev = re.search(r"GPU build: backend (\w+)", open(info0, errors="replace").read()) if os.path.exists(info0) else None
        dev = bool(m_dev) and m_dev.group(1) != "HOST"
        f0 = outputs(dirs["0"])
        if f0 != outputs(dirs["1"]) or f0 != outputs(dirs["2"]):
            bad.append("%s: output file lists differ" % name)
            continue
        worst = 0.0
        for f in f0:
            p0, p1, p2 = (os.path.join(dirs[m], f) for m in ("0", "1", "2"))
            # HDF5 files carry metadata (times): compare their contents (h5dump, as compare_schemes.py does)
            same = (lambda x, y: h5_text(x) == h5_text(y)) if f.endswith(".h5") else \
                   (lambda x, y: filecmp.cmp(x, y, shallow=False))
            if not same(p0, p2):
                # byte for byte on the CPU and GPU=host; on a device backend (CUDA, ONEMKL, HIP) re-runs are not bitwise
                # (device reductions), so there confout byte-identical and the other numbers to TOL_DEV
                n0, n2 = numbers(p0), numbers(p2)
                dev_ok = dev and not os.path.basename(f).startswith("confout") and not f.endswith(".h5") and len(n0) == len(n2) \
                    and max([abs(x - y)/max(1.0, abs(y)) for x, y in zip(n2, n0)] + [0.0]) <= TOL_DEV
                if not dev_ok:
                    bad.append("%s: %s mode 2 differs from mode 0" % (name, f))
            if same(p0, p1):
                continue
            if os.path.basename(f).startswith("confout"):
                bad.append("%s: %s differs between tracked and LU" % (name, f))
                continue
            v0, v1 = numbers(p0), numbers(p1)
            if len(v0) != len(v1):
                bad.append("%s: %s layout differs" % (name, f))
                continue
            d = max([abs(x - y) / max(1.0, abs(x), abs(y)) for x, y in zip(v0, v1)] + [0.0])
            worst = max(worst, d)
            if not d <= TOL:
                bad.append("%s: %s tracked vs LU differs by %.3g" % (name, f, d))
        print("%-16s tracked vs LU: max scaled difference %.3g" % (name, worst))

    if not a.no_tracked:
        for label, env, msg in (("fault-det", {"ALF_CGR_LU": "2", "ALF_CGR_FAULT": "det"}, "tracked det phase differs"),
                                ("invalid-mode", {"ALF_CGR_LU": "7"}, "ALF_CGR_LU=7 invalid")):
            d, rc = fresh(label, CASES[0][1], env)
            out = open(os.path.join(d, "stdout"), errors="replace").read()
            ok = rc != 0 and msg in out
            print("%-16s exit %d: %s" % (label, rc, "stopped as expected" if ok else "NOT stopped"))
            if not ok:
                bad.append(label)
    for b in bad:
        print("FAIL:", b)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
