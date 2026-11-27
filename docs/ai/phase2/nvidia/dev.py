#!/usr/bin/env python3
"""Device checks of the CUDA backend on Helma (P2.4 step 2, not committed).

usage: dev.py <suite> <exe> <workdir> [--ref <exe>]
suites: traj, g2, phase, speed. Each run directory gets parameters/seeds from <start> with edits; results are printed as
one line per case. Compares: max |x-y|/max(1,|y|) over all numbers of each output file (as gpu_e2e.py), confout identity."""
import os, re, sys, shutil, subprocess, time, filecmp

START = os.environ.get("DEV_START", os.path.expanduser("~/p2_nvidia/start"))
NUM = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")


def prep(d, edits):
    if os.path.exists(d):
        shutil.rmtree(d)
    os.makedirs(d)
    for f in ("parameters", "seeds"):
        shutil.copy(os.path.join(START, f), os.path.join(d, f))
    p = os.path.join(d, "parameters")
    t = open(p).read()
    for key, val in edits.items():
        t, n = re.subn(r"(?m)^(%s\s*=\s*)[^\s!]+" % key, r"\g<1>%s" % val, t, count=1)
        if n != 1:
            t, n = re.subn(r"(?m)^(Ltau\s*=)", "%s = %s\n\\g<1>" % (key, val), t, count=1)
        if n != 1:
            sys.exit("cannot set %s" % key)
    open(p, "w").write(t)


def run(exe, d, env_extra, mpiexec=None):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    for k in [k for k in env if k.startswith("ALF_GPU_")]:
        del env[k]
    env.update(env_extra)
    cmd = (mpiexec + [exe]) if mpiexec else [exe]
    t0 = time.time()
    with open(os.path.join(d, "stdout"), "w") as out:
        rc = subprocess.call(cmd, cwd=d, env=env, stdout=out, stderr=subprocess.STDOUT)
    return rc, time.time() - t0


def nums(f):
    return [float(x.replace("D", "E").replace("d", "e")) for x in NUM.findall(open(f, errors="replace").read())]


def numdiff(fa, fb):
    a, b = nums(fa), nums(fb)
    if len(a) != len(b):
        return float("inf")
    return max([abs(x - y) / max(1.0, abs(y)) for x, y in zip(a, b)] or [0.0])


def outfiles(d):
    return [f for f in sorted(os.listdir(d)) if f not in ("info", "info_git", "stdout", "parameters", "seeds")
            and not os.path.isdir(os.path.join(d, f))]


def compare(a, b):
    """(confout identical, max numdiff, worst file, max over tau files)"""
    worst, wf, tau = 0.0, "", 0.0
    for f in outfiles(b):
        if f.startswith("confout") or not os.path.exists(os.path.join(a, f)):
            continue
        x = numdiff(os.path.join(a, f), os.path.join(b, f))
        if "_tau" in f:
            tau = max(tau, x)
        if x > worst:
            worst, wf = x, f
    conf = all(filecmp.cmp(os.path.join(a, f), os.path.join(b, f), shallow=False)
               for f in outfiles(b) if f.startswith("confout"))
    return conf, worst, wf, tau


def info(d, pat):
    t = open(os.path.join(d, "info"), errors="replace").read()
    m = re.search(pat, t)
    return m.group(1) if m else None


def prec(d):
    return info(d, r"Precision Green\s+Mean, Max\s*:\s*(\S+)")


def acc(d):
    return info(d, r"Acceptance\s*[:=]?\s*(\S+)")


def gpu_lines(d):
    return [l.strip() for l in open(os.path.join(d, "info"), errors="replace") if "GPU" in l and ("ACTIVE" in l or "INACTIVE" in l or "DEVICE" in l or "HOST" in l or "phase" in l)]


TRAJ = [  # label, edits
    ("mz12-std", {"L1": 12, "L2": 12, "Beta": "10.d0", "NSweep": 10, "NBin": 4}),
    ("mz12-delayed16", {"L1": 12, "L2": 12, "Beta": "10.d0", "NSweep": 10, "NBin": 4, "Sequential_scheme": '"delayed"', "N_delay": 16}),
    ("mz12-submatrix16", {"L1": 12, "L2": 12, "Beta": "10.d0", "NSweep": 10, "NBin": 4, "Sequential_scheme": '"submatrix"', "N_delay": 16}),
    ("mz12-ltau0", {"L1": 12, "L2": 12, "Beta": "10.d0", "NSweep": 10, "NBin": 4, "Ltau": 0}),
    ("su2-12-std", {"L1": 12, "L2": 12, "Beta": "10.d0", "NSweep": 10, "NBin": 4, "Mz": ".F."}),
    ("mz6-std", {"L1": 6, "L2": 6, "Beta": "10.d0", "NSweep": 20, "NBin": 4}),
    ("su2-6-nwrap25", {"L1": 6, "L2": 6, "Beta": "10.d0", "NSweep": 20, "NBin": 4, "Mz": ".F.", "Nwrap": 25}),
]


def suite_traj(exe, w, ref):
    print("case | confout(active=disabled) | max rel diff (file) | tau files | replay bitwise | CGR_LU=2 line | PrecG active/disabled")
    sel = os.environ.get("DEV_CASES")
    for lab, ed in [c for c in TRAJ if not sel or c[0] in sel.split(",")]:
        d = {k: os.path.join(w, lab, k) for k in ("active", "disabled", "replay", "check", "none")}
        envs = {"active": {}, "disabled": {"ALF_GPU_DISABLE": "1"}, "replay": {},
                "check": {"ALF_GPU_DEBUG": "1", "ALF_GPU_CGR_LU": "2"}}
        for k, e in envs.items():
            prep(d[k], ed)
            rc, _ = run(exe, d[k], e)
            if rc:
                print(lab, k, "rc", rc); break
        conf, worst, wf, tau = compare(d["active"], d["disabled"])
        rconf, rworst, _, _ = compare(d["replay"], d["active"])
        same = rconf and rworst == 0.0 and all(filecmp.cmp(os.path.join(d["replay"], f), os.path.join(d["active"], f), shallow=False) for f in outfiles(d["active"]))
        chk = [l for l in open(os.path.join(d["check"], "info"), errors="replace") if "compared" in l or "max |phase" in l]
        extra = ""
        if ref:
            prep(d["none"], ed); run(ref, d["none"], {})
            nconf, nworst, _, _ = compare(d["none"], d["disabled"])
            extra = " | GPU=none vs disabled: conf %s, max %.1e" % (nconf, nworst)
        print("%s | %s | %.2e (%s) | %.2e | %s | %s | %s/%s%s" % (lab, conf, worst, wf, tau, same, " ".join(x.strip() for x in chk)[:160],
                                                         prec(d["active"]), prec(d["disabled"]), extra), flush=True)
        if lab == TRAJ[0][0]:
            print("  info GPU lines:", gpu_lines(d["active"]))


def suite_g2(exe, w, ref):
    print("L Nwrap | PrecG disabled | CUDA HOPD=0 | CUDA HOPD=1")
    for L in (6, 12, 24):
        for nw in (10, 25):
            ed = {"L1": L, "L2": L, "Beta": "10.d0", "NSweep": 5, "NBin": 2, "Nwrap": nw, "Mz": ".F."}
            r = []
            for k, e in (("dis", {"ALF_GPU_DISABLE": "1"}), ("h0", {}), ("h1", {"ALF_GPU_DEBUG": "1", "ALF_GPU_HOPD": "1"})):
                dd = os.path.join(w, "L%d_nw%d" % (L, nw), k); prep(dd, ed); rc, _ = run(exe, dd, e)
                r.append(prec(dd) if rc == 0 else "rc%d" % rc)
            print("%d %d | %s | %s | %s" % (L, nw, *r), flush=True)


def suite_phase(exe, w, ref):
    print("L | acceptance | CGR check | phase line | time s")
    for L in (56, 60, 64):
        ed = {"L1": L, "L2": L, "Beta": "6.4d0", "NSweep": 2, "NBin": 1, "Ltau": 0, "ham_U": "5.d0", "Symm": ".F."}
        dd = os.path.join(w, "L%d" % L); prep(dd, ed)
        rc, t = run(exe, dd, {"ALF_GPU_DEBUG": "1", "ALF_GPU_CGR_LU": "2"})
        chk = [l.strip() for l in open(os.path.join(dd, "info"), errors="replace") if "compared" in l or "max |phase" in l]
        print("%d | %s | rc %d | %s | %.0f" % (L, acc(dd), rc, " ".join(chk)[:200], t), flush=True)


def suite_speed(exe, w, ref):
    """s per sweep from the info file (CPU Time / NSweep); Symm = .F. on both codes (start-up check P1 of Symm = .T. is slow)"""
    print("L case | CPU Time s | s/sweep | PrecG mean | wall s")
    ns = 3
    for L in (24, 32):
        ed = {"L1": L, "L2": L, "Beta": "32.d0", "NSweep": ns, "NBin": 1, "Ltau": 0, "Nwrap": 10, "Symm": ".F.",
              "Sequential_scheme": '"submatrix"', "N_delay": 32}
        cases = [("cuda", {}), ("cuda-lu", {"ALF_GPU_DEBUG": "1", "ALF_GPU_CGR_LU": "1"})]
        if L == 24:
            cases.append(("cuda-hopd", {"ALF_GPU_DEBUG": "1", "ALF_GPU_HOPD": "1"}))
        for k, e in cases:
            for rep in range(2):
                dd = os.path.join(w, "L%d_%s_%d" % (L, k, rep)); prep(dd, ed); rc, t = run(exe, dd, e)
                ct = info(dd, r"CPU Time\s*:\s*(\S+)") if rc == 0 else None
                print("%d %s | %s | %s | %s | %.0f rc %d" % (L, k, ct, "%.3f" % (float(ct)/ns) if ct else "-", prec(dd), t, rc), flush=True)


def suite_mpi(exe, w, ref):
    """2 MPI ranks on 2 GPUs (device binding: rank mod number of devices): active vs disabled, same seeds"""
    mpi = ["mpirun", "-np", "2", "--bind-to", "none"]
    ed = {"L1": 12, "L2": 12, "Beta": "10.d0", "NSweep": 10, "NBin": 4}
    d = {k: os.path.join(w, k) for k in ("active", "disabled")}
    for k, e in (("active", {}), ("disabled", {"ALF_GPU_DISABLE": "1"})):
        prep(d[k], ed); rc, t = run(exe, d[k], e, mpi); print(k, "rc", rc, "%.0f s" % t)
    conf, worst, wf, tau = compare(d["active"], d["disabled"])
    print("mpi2 | confout %s | max %.2e (%s) | tau %.2e | PrecG %s/%s" % (conf, worst, wf, tau, prec(d["active"]), prec(d["disabled"])))
    print("  info GPU lines:", gpu_lines(d["active"]))


if __name__ == "__main__":
    suite, exe, w = sys.argv[1], os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3])
    ref = os.path.abspath(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[4] == "--ref" else None
    os.makedirs(w, exist_ok=True)
    {"mpi": suite_mpi, "traj": suite_traj, "g2": suite_g2, "phase": suite_phase, "speed": suite_speed}[suite](exe, w, ref)
