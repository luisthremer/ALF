#!/usr/bin/env python3
"""Job 4 (P2.8.5): CPU differential runs.  usage: diff_runs.py SERIAL_EXE MPI_EXE IFX_EXE WORK OUT
  mpi      mpirun -np 1/2/3 of the MPI build against the serial build (same Start seeds)
  restart  one run of 2N bins against N bins + restart (confout -> confin) + N bins
  compiler gfortran build against ifx build (round-off expected)
Verdict per comparison: IDENTICAL (every file byte-identical, info/stdout excluded), ROUNDOFF max scaled diff x (confout identical), DIFFERENT."""
import filecmp, os, shutil, sys
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "testsuite", "GPU.tests")); sys.path.insert(0, os.path.join(REPO, "docs", "ai", "phase2", "p25_runs"))
import gpu_e2e as e2e, p25_run as p25
START = os.path.join(REPO, "Scripts_and_Parameters_files", "Start")
KONDO = {"ham_name": '"Kondo"', "Model": '"Kondo"', "Lattice_type": '"Bilayer_square"', "ham_Uf": "1.d0"}
TV = {"ham_name": '"tV"', "Model": '"tV"', "N_SUN": "1", "ham_V": "-1.d0", "ham_T2": "0.d0", "ham_Tperp": "0.d0", "Beta": "2.d0", "Ltau": "0"}
SP = {"ham_name": '"Spin_Peierls"', "Model": '"Spin_Peierls"', "L1": "6", "L2": "6", "Beta": "2.d0", "ham_U": "1.5d0", "ham_Lambda": "0.1d0", "ham_omega0": "0.5d0", "N_SUN": "2", "N_FL": "1"}
CASES = {"hub-mz": {}, "hub-su2": {"Mz": ".F."}, "hub-dly4": {"Sequential_scheme": '"delayed"', "N_delay": "4"}, "kondo": KONDO, "tV": TV, "sp": SP}
BASE = {"L1": "4", "L2": "4", "NSweep": "10", "NBin": "4"}
def ed(c, **kw):
    e = dict(BASE); e.update(CASES[c]); e.update(kw); return e
SEEDS = [l.strip() for l in open(os.path.join(START, "seeds")) if l.strip()]
def go(exe, d, e, mpi=0, line=1, env=None):
    p25.prepare(d, START, e)
    if line != 1: open(os.path.join(d, "seeds"), "w").write("\n".join(SEEDS[line - 1:]) + "\n")   # serial chain on the seed that rank r of np=n gets
    return p25.run(exe, d, env or {}, mpi, 3600)[0]
def verdict(a, b, only=None, fb_name=None):
    bad, mx, n = [], 0.0, 0
    for f in e2e.files(a):
        if f.startswith(".") or f in ("parameters", "seeds") or (only and f not in only): continue
        fb = os.path.join(b, f)
        if not os.path.exists(fb): bad.append(f + " missing"); continue
        n += 1
        if filecmp.cmp(os.path.join(a, f), fb, shallow=False): continue
        if f.startswith("confout"): bad.append(f + " differs"); continue
        dd = e2e.numdiff(os.path.join(a, f), fb)
        if dd is None: bad.append(f + " layout"); continue
        mx = max(mx, dd[1]); 
    if bad: return "DIFFERENT", "; ".join(bad[:4])
    return ("IDENTICAL", "%d files" % n) if mx == 0 else ("ROUNDOFF", "max scaled diff %.2e, %d files, confout identical" % (mx, n))
def main():
    se, me, ie, W, out = sys.argv[1:6]; rows = ["case\tcomparison\tverdict\tdetail"]
    def rec(c, cmp_, v): rows.append("%s\t%s\t%s\t%s" % (c, cmp_, v[0], v[1])); print(rows[-1], flush=True)
    if os.environ.get("DIFF_SHORT"):      # gfortran vs ifx on a 1-sweep run: the chain has not yet diverged, the observables show round-off
        for c in CASES:
            r1, r2 = os.path.join(W, c, "short_gnu"), os.path.join(W, c, "short_ifx")
            go(se, r1, ed(c, NSweep="1", NBin="2")); rc = go(ie, r2, ed(c, NSweep="1", NBin="2"), env={"LD_LIBRARY_PATH": "/opt/intel/oneapi/mkl/2026.1/lib:/opt/intel/oneapi/compiler/2026.1/lib"})
            rec(c, "ifx vs gfortran, NSweep 1 NBin 2", ("FAIL", "rc=%d" % rc) if rc else verdict(r1, r2))
        os.makedirs(out, exist_ok=True); open(os.path.join(out, "diff_runs_short.tsv"), "w").write("\n".join(rows) + "\n"); return
    for c in CASES:
        ref = os.path.join(W, c, "serial")
        if go(se, ref, ed(c)) != 0: rec(c, "serial run", ("FAIL", "rc != 0")); continue
        for n in (1, 2, 3):
            d = os.path.join(W, c, "mpi%d" % n)
            rc = go(me, d, ed(c), mpi=n)
            # ALF gives rank 0 the seed on line n of the seeds file and rank r > 0 the one on line n - r (Set_Random_number_Generator)
            for r in range(n):
                line = n - r; sref = os.path.join(W, c, "serial_seed%d" % line)
                if not os.path.exists(sref + "/.ok"):
                    go(se, sref, ed(c), line=line); open(sref + "/.ok", "w").write("1")
                if rc: rec(c, "MPI np=%d rank %d" % (n, r), ("FAIL", "rc=%d" % rc)); continue
                if r == 0: v = verdict(sref, d)       # rank 0 writes every file; observables of rank 0 only?  (compared as they are)
                else:
                    same = filecmp.cmp(os.path.join(sref, "confout_0"), os.path.join(d, "confout_%d" % r), shallow=False)
                    v = ("IDENTICAL", "confout_%d = serial chain of seed line %d" % (r, line)) if same else ("DIFFERENT", "confout_%d differs" % r)
                rec(c, "MPI np=%d rank %d (seed line %d) vs serial run of that seed" % (n, r, line), v)
        # restart
        full = os.path.join(W, c, "full"); go(se, full, ed(c, NBin="8"))
        half = os.path.join(W, c, "restart"); rc = go(se, half, ed(c, NBin="4"))
        for f in os.listdir(half):
            if f.startswith("confout_"): shutil.copy(os.path.join(half, f), os.path.join(half, "confin_" + f[8:]))
        rc2 = p25.run(se, half, {}, 0, 3600)[0]
        rec(c, "restart (4+4 bins) vs single run (8 bins)", ("FAIL", "rc=%d/%d" % (rc, rc2)) if rc or rc2 else verdict(full, half))
        # compiler
        d = os.path.join(W, c, "ifx"); rc = go(ie, d, ed(c), env={"LD_LIBRARY_PATH": "/opt/intel/oneapi/mkl/2026.1/lib:/opt/intel/oneapi/compiler/2026.1/lib"})
        rec(c, "ifx 2025.3 vs gfortran 14.2", ("FAIL", "rc=%d" % rc) if rc else verdict(ref, d))
    os.makedirs(out, exist_ok=True); open(os.path.join(out, "diff_runs.tsv"), "w").write("\n".join(rows) + "\n")
main()
