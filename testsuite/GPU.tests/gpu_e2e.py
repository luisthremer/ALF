#!/usr/bin/env python3
"""End-to-end check of the GPU build (P2.4 steps 1b, 1c, 1d): Hubbard 4x4, same seeds, for each case

  active    ALF.out of the GPU build, GPU path active (device wraps, local updates and stabilisation; ALF_GPU_SELFTEST=1
            additionally runs the mirror self-test at start-up)
  disabled  the same binary with ALF_GPU_DISABLE=1 (CPU path, mirror not built)
  none      optional (--ref-exe): ALF.out of a GPU=none build (must be byte-identical to 'disabled': CPU path unchanged)

Requirements per case: the info file of the active run says "ACTIVE (mirror, wraps and updates on the device path)" (cases that
must stay inactive: an INACTIVE line with the expected reason, and then every output file byte-identical to 'disabled');
confout_0 byte-identical to 'disabled' (the Markov chain is the CPU one for these short runs); every other output file equal to
'disabled' to max(1e-11, 100 x Precision Green mean of the 'disabled' run) (|x - y| / max(1, |y|)).  The 1e-11 is the round-off of the device
kernels (e.g. FMA).  Since step 1d the device decomposition (ZGEQRF + ZUNMQR, column-norm pre-pivot) differs from the CPU one (ZGEQP3,
ZUNGQR) in the last bits of G at every stabilisation point, and the wrap instability amplifies that by the factor that also makes the
Precision Green line of the info file (mean deviation of the wrapped from the recomputed G) non-zero: 6e-11 for the Mz-Hubbard run
below (observable differences 1e-10 .. 1e-9), 1e-14 for the SU(2) runs (differences 1e-14).  The criterion ties the tolerance to that
noise floor of the run itself.  (The Precision Green means of 'active' and 'disabled' are printed: with Nwrap >= 25, where the
wrap is unstable, the device wraps and updates of step 1c give a mean 4 - 6x that of the CPU, already without the device stabilisation.)
Negative control (ALF_GPU_FAULT=upd, a deliberately wrong device
update): the observables must then differ by more than 1e-8, i.e. the comparison can see an error.
Cases: Mz-Hubbard (N_FL = 2) standard / delayed / submatrix with N_delay = 4, Ltau = 0, SU(2) Hubbard (Mz = .F., N_FL = 1) standard /
delayed, small, odd and large Nwrap (2, 7, 25, 50 of Ltrot = 50), timed N_delay = -1 (refused).  The run of the build is Ltau = 1
(time-displaced Green function: udvst exported to the host TAU_M) unless the case sets Ltau = 0.  .h5 files are compared with h5diff
if it is on the PATH, otherwise skipped with a note."""
import argparse, filecmp, os, re, shutil, subprocess, sys

TOL = 1.0e-11          # round-off of the device kernels: |x - y| / max(1, |y|) over all numbers of an output file
FAULT_SEEN = 1.0e-8    # a negative control must change the observables by more than this
NUM = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")
NUM_WS = re.compile(r"\s*" + NUM.pattern)      # a number with its leading blanks (a minus sign takes the place of a blank)
ACTIVE_TEXT = "ACTIVE (mirror; wraps, updates and stabilisation on the device path)"

# label, edits of the parameter file, env of the active run, expectation ("ok", "inactive:<text>", "fault" or "abort")
CASES = [
    ("mz-standard", {}, {"ALF_GPU_SELFTEST": "1"}, "ok"),
    ("mz-delayed-N4", {"Sequential_scheme": '"delayed"', "N_delay": "4"}, {}, "ok"),
    ("mz-submatrix-N4", {"Sequential_scheme": '"submatrix"', "N_delay": "4"}, {}, "ok"),
    ("su2-standard", {"Mz": ".F."}, {}, "ok"),
    ("su2-delayed-N4", {"Mz": ".F.", "Sequential_scheme": '"delayed"', "N_delay": "4"}, {}, "ok"),
    ("mz-ltau0", {"Ltau": "0"}, {}, "ok"),
    ("su2-nwrap2", {"Mz": ".F.", "Nwrap": "2"}, {}, "ok"),
    ("su2-nwrap7", {"Mz": ".F.", "Nwrap": "7"}, {}, "ok"),
    ("su2-nwrap25", {"Mz": ".F.", "Nwrap": "25"}, {}, "ok"),
    ("su2-nwrap50", {"Mz": ".F.", "Nwrap": "50"}, {}, "ok"),
    ("timed-N-1", {"Sequential_scheme": '"delayed"', "N_delay": "-1"}, {}, "inactive:N_delay = -1"),
    ("fault-upd", {}, {"ALF_GPU_DEBUG": "1", "ALF_GPU_FAULT": "upd"}, "fault"),
    # a wrong pivot permutation of the device decomposition (ALF_GPU_FAULT_STAB=perm) makes G wrong at the first stabilisation: the run
    # must stop with an error (expectation "abort": the active run exits with a non-zero return code)
    ("fault-stab-perm", {}, {"ALF_GPU_DEBUG": "1", "ALF_GPU_FAULT_STAB": "perm"}, "abort"),
    # ALF_GPU_FAULT=nsig (stale device copy of the field) cannot be seen on the host backend: host and "device" array are the
    # same memory there, so it is not a case here (to be added for the vendor backends, step 2).
]


def prepare(d, start, edits):
    if os.path.exists(d):
        shutil.rmtree(d)
    os.makedirs(d)
    for f in ("parameters", "seeds"):
        shutil.copy(os.path.join(start, f), os.path.join(d, f))
    p = os.path.join(d, "parameters")
    t = open(p).read()
    for key, val in (("L1", 4), ("L2", 4), ("NSweep", 6), ("NBin", 2)):
        t, n = re.subn(r"(?m)^(%s\s*=\s*)\d+" % key, r"\g<1>%d" % val, t, count=1)
        if n != 1:
            sys.exit("cannot set %s in parameters" % key)
    for key, val in edits.items():
        if key in ("Mz", "Nwrap", "Ltau"):       # existing keys: replace the value
            t, n = re.subn(r"(?m)^(%s\s*=\s*)[^\s!]+" % key, r"\g<1>%s" % val, t, count=1)
        else:   # new key of &VAR_QMC, inserted before Ltau
            t, n = re.subn(r"(?m)^(Ltau\s*=)", "%s = %s\n\\g<1>" % (key, val), t, count=1)
        if n != 1:
            sys.exit("cannot set %s in parameters" % key)
    open(p, "w").write(t)


def run(exe, d, env_extra, mpiexec):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    for k in [k for k in env if k.startswith("ALF_GPU_") or k.startswith("ALF_LA_TEST_")]:
        del env[k]
    env.update(env_extra)
    cmd = (mpiexec.split() + [exe]) if mpiexec else [exe]
    with open(os.path.join(d, "stdout"), "w") as out:
        rc = subprocess.call(cmd, cwd=d, env=env, stdout=out, stderr=subprocess.STDOUT)
    return rc


def files(d):
    return [f for f in sorted(os.listdir(d)) if f not in ("info", "info_git", "stdout") and not os.path.isdir(os.path.join(d, f))]


def numdiff(fa, fb):
    """(max abs difference, max of |x-y|/max(1,|y|)) over the numbers of two text files; None if the layout differs."""
    ta, tb = open(fa, errors="replace").read(), open(fb, errors="replace").read()
    if NUM_WS.sub("#", ta) != NUM_WS.sub("#", tb):
        return None
    va = [float(x.replace("D", "E").replace("d", "e")) for x in NUM.findall(ta)]
    vb = [float(x.replace("D", "E").replace("d", "e")) for x in NUM.findall(tb)]
    return (max([abs(x - y) for x, y in zip(va, vb)] + [0.0]),
            max([abs(x - y) / max(1.0, abs(y)) for x, y in zip(va, vb)] + [0.0]))


def green_precision(d):
    """(mean, max) of the 'Precision Green' line of the info file of run directory d, or None"""
    m = re.search(r"Precision Green\s+Mean, Max\s*:\s*(\S+)\s+(\S+)", open(os.path.join(d, "info"), errors="replace").read())
    return (float(m.group(1).replace("D", "E")), float(m.group(2).replace("D", "E"))) if m else None


def compare(a, b, strict, tol=TOL):
    """(files that are not within the rules, {file: (max abs diff, max scaled diff)} of the numeric files).  strict: every file
    byte-identical; else confout_* byte-identical and every other file equal to tol (scaled)."""
    bad, diffs = [], {}
    for f in files(a):
        fa, fb = os.path.join(a, f), os.path.join(b, f)
        if not os.path.exists(fb):
            bad.append(f + " (missing)")
        elif f.endswith(".h5"):
            if shutil.which("h5diff"):
                if subprocess.call(["h5diff", fa, fb], stdout=subprocess.DEVNULL) != 0:
                    bad.append(f)
            else:
                print("note: no h5diff, %s not compared" % f)
        elif filecmp.cmp(fa, fb, shallow=False):
            diffs[f] = (0.0, 0.0)
        elif strict or f.startswith("confout"):
            bad.append(f)
        else:
            dd = numdiff(fa, fb)
            if dd is None or dd[1] > tol:
                bad.append("%s (%s)" % (f, "layout differs" if dd is None else "max abs diff %.3g, scaled %.3g" % dd))
            diffs[f] = dd
    return bad, diffs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--start", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--ref-exe", default=None)
    ap.add_argument("--mpiexec", default="")
    a = ap.parse_args()
    fail = 0
    for label, edits, env, expect in CASES:
        dirs = {}
        runs = {"active": (a.exe, env), "disabled": (a.exe, {"ALF_GPU_DISABLE": "1"})}
        if a.ref_exe:
            runs["none"] = (a.ref_exe, {})
        for lab, (exe, e) in runs.items():
            dirs[lab] = os.path.join(a.work, label, lab)
            prepare(dirs[lab], a.start, edits)
            rc = run(exe, dirs[lab], e, a.mpiexec)
            if lab == "active":
                rc_active = rc
            if rc != 0 and not (expect == "abort" and lab == "active"):
                sys.exit("FAIL: %s: run '%s' exited with %d (see %s/stdout)" % (label, lab, rc, dirs[lab]))
        if expect == "abort":
            ok = rc_active != 0
            print("%-16s active run (fault) exit code %d: %s" % (label, rc_active, "OK (run stopped)" if ok else "FAILED (run NOT stopped)"))
            fail += not ok
            continue
        info = open(os.path.join(dirs["active"], "info"), errors="replace").read()
        out = open(os.path.join(dirs["active"], "stdout"), errors="replace").read()
        dinfo = open(os.path.join(dirs["disabled"], "info"), errors="replace").read()
        bad_case = []
        if "INACTIVE, reason: disabled by ALF_GPU_DISABLE=1" not in dinfo:
            bad_case.append("info of the disabled run has no INACTIVE line")
        if expect.startswith("inactive:"):
            m = re.search(r"INACTIVE, reason: (.*)", info)
            if not m or expect[9:] not in m.group(1):
                bad_case.append("expected INACTIVE with reason containing '%s'" % expect[9:])
            reason = m.group(1).strip() if m else "?"
        else:
            reason = ""
            if ACTIVE_TEXT not in info:
                bad_case.append("info of the active run has no ACTIVE line")
            if "ALF_GPU_SELFTEST" in env and "GPU self-test PASSED" not in out:
                bad_case.append("stdout of the active run has no 'GPU self-test PASSED'")
            if expect == "ok" and "GPU update path: ACTIVE" not in info:
                bad_case.append("info of the active run has no 'GPU update path: ACTIVE' line")
        if not os.path.exists(os.path.join(dirs["active"], "confout_0")) and not os.path.exists(os.path.join(dirs["active"], "data.h5")):
            bad_case.append("no confout_0 / data.h5 in the active run")
        gp_act, gp_dis = green_precision(dirs["active"]), green_precision(dirs["disabled"])
        tol = TOL if gp_dis is None else max(TOL, 100.0 * gp_dis[0])
        bad, diffs = compare(dirs["active"], dirs["disabled"], strict=expect.startswith("inactive:"), tol=tol)
        scal = ", ".join("%s %.3g / %.3g" % (f[:-5], diffs[f][0], diffs[f][1]) for f in ("Ener_scal", "Kin_scal", "Pot_scal")
                         if diffs.get(f) is not None)
        if expect == "fault":
            # negative control: the comparison must see the injected error
            seen = [f for f in files(dirs["active"]) if not f.startswith("confout") and (diffs.get(f) is None or diffs.get(f)[1] > FAULT_SEEN)]
            seen += [f for f in bad if not f.startswith("confout")]
            status = "OK (fault seen)" if seen else "FAILED (fault NOT seen)"
            if not seen:
                bad_case.append("fault not visible in the observables")
            print("%-16s active vs disabled: max abs / scaled diff (%s) %s" % (label, scal, status))
        else:
            status = "OK" if not bad and not bad_case else "FAILED"
            conf = "confout byte-identical" if "confout_0" in diffs and diffs["confout_0"][0] == 0.0 else "confout differs"
            gp = (", Precision Green mean active %.2g / disabled %.2g, tol %.2g" % (gp_act[0], gp_dis[0], tol)) if gp_act and gp_dis else ""
            print("%-16s active vs disabled: %s; max abs / scaled diff (%s)%s %s%s" % (label, conf, scal, gp, status,
                                                                                       (" [reason: %s]" % reason) if reason else ""))
            for b in bad:
                print("    not within the rules: " + b)
            fail += bool(bad)
        for b in bad_case:
            print("    " + b)
        fail += bool(bad_case)
        if a.ref_exe:
            bad_n = compare(dirs["none"], dirs["disabled"], strict=True)[0]
            print("%-16s none   vs disabled: %s" % (label, "byte-identical" if not bad_n else "DIFFER: " + " ".join(bad_n)))
            fail += bool(bad_n)
    print("PASSED" if not fail else "FAILED")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
