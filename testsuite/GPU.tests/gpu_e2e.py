#!/usr/bin/env python3
"""End-to-end check of the GPU build (P2.4 step 1b): Hubbard 4x4 with the same seeds, run

  active    ALF.out of the GPU build with ALF_GPU_SELFTEST=1 (mirror built and self-tested, not used by the sweep)
  disabled  the same binary with ALF_GPU_DISABLE=1 (CPU path, mirror not built)
  none      optional (--ref-exe): ALF.out of a GPU=none build

and require that every output file except the info files and stdout is byte-identical, that the info file of the
active run says "ACTIVE (mirror built, not yet used by the sweep)" and that stdout says "GPU self-test PASSED".
The sweep does not use the mirror in step 1b, so the agreement must be exact.  .h5 files are compared with h5diff
if it is on the PATH (HDF5 files carry no usable byte identity), otherwise skipped with a note."""
import argparse, filecmp, os, re, shutil, subprocess, sys


def prepare(d, start):
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


def same(a, b):
    """Names of the files that differ between two run directories (info, info_git and stdout are excluded)."""
    bad = []
    for f in sorted(os.listdir(a)):
        if f in ("info", "info_git", "stdout") or os.path.isdir(os.path.join(a, f)):
            continue
        fa, fb = os.path.join(a, f), os.path.join(b, f)
        if not os.path.exists(fb):
            bad.append(f + " (missing)")
        elif f.endswith(".h5"):
            if shutil.which("h5diff"):
                if subprocess.call(["h5diff", fa, fb], stdout=subprocess.DEVNULL) != 0:
                    bad.append(f)
            else:
                print("note: no h5diff, %s not compared" % f)
        elif not filecmp.cmp(fa, fb, shallow=False):
            bad.append(f)
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--start", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--ref-exe", default=None)
    ap.add_argument("--mpiexec", default="")
    a = ap.parse_args()
    runs = {"active": (a.exe, {"ALF_GPU_SELFTEST": "1"}), "disabled": (a.exe, {"ALF_GPU_DISABLE": "1"})}
    if a.ref_exe:
        runs["none"] = (a.ref_exe, {})
    dirs = {}
    for lab, (exe, env) in runs.items():
        dirs[lab] = os.path.join(a.work, lab)
        prepare(dirs[lab], a.start)
        rc = run(exe, dirs[lab], env, a.mpiexec)
        if rc != 0:
            sys.exit("FAIL: run '%s' exited with %d (see %s/stdout)" % (lab, rc, dirs[lab]))
    fail = 0
    info = open(os.path.join(dirs["active"], "info"), errors="replace").read()
    out = open(os.path.join(dirs["active"], "stdout"), errors="replace").read()
    if "ACTIVE (mirror built, not yet used by the sweep)" not in info:
        print("FAIL: info of the active run has no ACTIVE line"); fail += 1
    if "GPU self-test PASSED" not in out:
        print("FAIL: stdout of the active run has no 'GPU self-test PASSED'"); fail += 1
    dinfo = open(os.path.join(dirs["disabled"], "info"), errors="replace").read()
    if "INACTIVE, reason: disabled by ALF_GPU_DISABLE=1" not in dinfo:
        print("FAIL: info of the disabled run has no INACTIVE line"); fail += 1
    if not os.path.exists(os.path.join(dirs["active"], "confout_0")) and not os.path.exists(os.path.join(dirs["active"], "data.h5")):
        print("FAIL: no confout_0 / data.h5 in the active run"); fail += 1
    for lab in [l for l in runs if l != "active"]:
        bad = same(dirs["active"], dirs[lab])
        print("active vs %-8s %s" % (lab, "byte-identical" if not bad else "DIFFER: " + " ".join(bad)))
        fail += bool(bad)
    print("PASSED" if not fail else "FAILED")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
