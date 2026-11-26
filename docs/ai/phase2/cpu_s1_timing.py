"""cpu-s1 timing: same ALF.out, ALF_CGR_LU=1 (LU) vs 0 (tracked), Hubbard Mz L = 8/12/16, standard and delayed.
Wall time per sweep (whole run / NSweep, start-up included), 3 repetitions, each run pinned to its own core.
usage: timing.py ALF.out Start/ workdir first_core [ncores]"""
import os, re, shutil, subprocess, sys, time, itertools, statistics
from concurrent.futures import ThreadPoolExecutor

exe, start, work = sys.argv[1:4]
c0 = int(sys.argv[4]); nc = int(sys.argv[5]) if len(sys.argv) > 5 else 6
cores = list(range(c0, c0 + nc))
SIZES = {8: 40, 12: 12, 16: 4}
jobs = list(itertools.product(range(3), SIZES, ("standard", "delayed"), ("1", "0")))


def go(i_job):
    i, (rep, L, scheme, mode) = i_job
    d = os.path.join(work, "L%d_%s_m%s_r%d" % (L, scheme, mode, rep))
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    shutil.copy(os.path.join(start, "seeds"), d)
    t = open(os.path.join(start, "parameters")).read()
    for k, v in (("L1", L), ("L2", L), ("NSweep", SIZES[L]), ("NBin", 1), ("Ltau", 0)):
        t = re.sub(r"(?m)^(%s\s*=\s*)\S+" % k, r"\g<1>%s" % v, t, count=1)
    t = re.sub(r"(?im)^(\s*&VAR_QMC\b.*)$",
               lambda m: '%s\nSequential_scheme = "%s"\nN_delay = 16' % (m.group(1), scheme), t)
    open(os.path.join(d, "parameters"), "w").write(t)
    env = dict(os.environ, OMP_NUM_THREADS="1", ALF_CGR_LU=mode)
    t0 = time.time()
    rc = subprocess.call(["taskset", "-c", str(cores[i % len(cores)]), exe], cwd=d, env=env,
                         stdout=open(os.path.join(d, "stdout"), "w"), stderr=subprocess.STDOUT)
    return (L, scheme, mode, rep, rc, (time.time() - t0) / SIZES[L])


with ThreadPoolExecutor(len(cores)) as ex:
    res = list(ex.map(go, enumerate(jobs)))
print("| L | N | scheme | LU s/sweep | tracked s/sweep | change |")
print("|---|---|---|---|---|---|")
for L in SIZES:
    for s in ("standard", "delayed"):
        lu = statistics.median([r[5] for r in res if r[0] == L and r[1] == s and r[2] == "1"])
        tr = statistics.median([r[5] for r in res if r[0] == L and r[1] == s and r[2] == "0"])
        print("| %d | %d | %s | %.4f | %.4f | %+.1f %% |" % (L, L * L, s, lu, tr, 100 * (tr - lu) / lu))
print("failed runs:", [r for r in res if r[4]])
