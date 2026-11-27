"""macOS BLAS benchmark: wall time per sweep and Precision Green for several ALF.out builds.
usage: mac_blas_bench.py Start/ workdir name=ALF.out [name=ALF.out ...]
Runs sequentially (one run at a time, so timings do not disturb each other), 1 and 4 threads."""
import os, re, shutil, subprocess, sys, time

start, work = sys.argv[1:3]
builds = [a.split("=", 1) for a in sys.argv[3:]]
CASES = {"Hubbard16": ("Hubbard", 16, 6), "Kondo12": ("Kondo", 12, 6)}
REPS = 2


def run(name, exe, case, nthr, rep):
    ham, L, nsw = CASES[case]
    d = os.path.join(work, "%s_%s_t%d_r%d" % (case, name, nthr, rep))
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    shutil.copy(os.path.join(start, "seeds"), d)
    t = open(os.path.join(start, "parameters")).read()
    t = re.sub(r'(?m)^(ham_name\s*=\s*)"\w+"', r'\g<1>"%s"' % ham, t, count=1)
    if ham == "Kondo":
        t = t.replace('Lattice_type = "Square" ', 'Lattice_type = "Bilayer_square"').replace(
            'Model        = "Hubbard"', 'Model        = "Kondo"')
    for k, v in (("L1", L), ("L2", L), ("NSweep", nsw), ("NBin", 1), ("Ltau", 0)):
        t = re.sub(r"(?m)^(%s\s*=\s*)\S+" % k, r"\g<1>%s" % v, t, count=1)
    open(os.path.join(d, "parameters"), "w").write(t)
    n = str(nthr)
    env = dict(os.environ, OMP_NUM_THREADS=n, OPENBLAS_NUM_THREADS=n, VECLIB_MAXIMUM_THREADS=n)
    t0 = time.time()
    rc = subprocess.call([exe], cwd=d, env=env, stdout=open(os.path.join(d, "stdout"), "w"),
                         stderr=subprocess.STDOUT)
    dt = (time.time() - t0) / nsw
    m = re.search(r"Precision Green\s+Mean, Max :\s*(\S+)\s+(\S+)", open(os.path.join(d, "info")).read()) \
        if os.path.exists(os.path.join(d, "info")) else None
    return rc, dt, (m.group(1), m.group(2)) if m else ("?", "?")


print("| case | build | threads | s/sweep (rep1, rep2) | Precision Green mean, max | rc |")
print("|---|---|---|---|---|---|")
for case in CASES:
    for nthr in (1, 4):
        for name, exe in builds:
            r = [run(name, exe, case, nthr, i) for i in range(REPS)]
            print("| %s | %s | %d | %.2f, %.2f | %s, %s | %s |" % (
                case, name, nthr, r[0][1], r[1][1], r[0][2][0], r[0][2][1], r[0][0]), flush=True)
