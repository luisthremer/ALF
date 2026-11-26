#!/usr/bin/env python3
"""P2.5 run evidence, GPU=host column.  One case = one parameter file, run twice with the same seeds by the same binary:
  active    default environment (GPU path if gpu_supported accepts the run)
  disabled  ALF_GPU_DISABLE=1 (CPU path)
Verdict per case (printed as one TSV line, appended to <work>/results.tsv):
  works    info says "GPU path (gpu_state): ACTIVE", both runs exit 0, confout byte-identical, every other output file equal to
           max(1e-11, 100 x Precision Green mean of the 'disabled' run) (same rule as testsuite/GPU.tests/gpu_e2e.py)
  refused  info says "INACTIVE, reason: ..." (reason quoted), both runs exit 0, every output file byte-identical to 'disabled'
  FAIL     anything else (exit code /= 0, active but outside the tolerance, inactive but files differ, unexpected state)
The column 'expect' is what GPU_FEATURES.md predicted; MISMATCH marks a prediction that the run did not confirm.
Usage: p25_run.py --exe ALF.out --start Scripts_and_Parameters_files/Start --work DIR [--mpi N] [--only label,label] cases.py
cases.py defines CASES = [(label, {namelist key: value}, expect, {env}, note)]; expect = "works" | "refused" | "refused:<text>" | "stops".
"""
import argparse, importlib.util, os, re, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "..", "testsuite", "GPU.tests"))
import gpu_e2e as e2e   # compare(), green_precision(), numdiff()

def compare(a, b, strict, tol):
    """like gpu_e2e.compare, but .h5 files are compared with h5diff -d tol (absolute; the values are of order 1) unless strict or confout"""
    import filecmp
    bad, diffs = [], {}
    for f in e2e.files(a):
        fa, fb = os.path.join(a, f), os.path.join(b, f)
        if not os.path.exists(fb):
            bad.append(f + " (missing)")
        elif f.endswith(".h5"):   # HDF5 files carry creation times: never byte-compared
            exact = strict or f.startswith("confout")
            if not shutil.which("h5diff"):
                print("note: no h5diff, %s not compared" % f)
            elif subprocess.call(["h5diff"] + ([] if exact else ["-d", repr(tol)]) + [fa, fb], stdout=subprocess.DEVNULL) != 0:
                bad.append(f + (" (h5diff exact)" if exact else " (h5diff -d %.1e)" % tol))
            else:
                diffs[f] = (0.0, 0.0)   # equal (exact or within h5diff -d tol), no number to report
        elif filecmp.cmp(fa, fb, shallow=False):
            diffs[f] = (0.0, 0.0)
        elif strict or f.startswith("confout"):
            bad.append(f)
        else:
            dd = e2e.numdiff(fa, fb)
            if dd is None or dd[1] > tol:
                bad.append("%s (%s)" % (f, "layout differs" if dd is None else "max abs diff %.3g, scaled %.3g" % dd))
            diffs[f] = dd
    return bad, diffs

def edit_parameters(path, edits):
    t = open(path).read()
    for key, val in edits.items():
        t, n = re.subn(r"(?mi)^(%s\s*=\s*)[^\s!]+" % re.escape(key), r"\g<1>%s" % val, t)   # every namelist that has the key
        if n == 0:   # new key of &VAR_QMC
            t, n = re.subn(r"(?m)^(Ltau\s*=)", "%s = %s\n\\g<1>" % (key, val), t, count=1)
        if n == 0:
            sys.exit("cannot set %s" % key)
    open(path, "w").write(t)

def prepare(d, start, edits, groups=0):
    if os.path.exists(d):
        shutil.rmtree(d)
    os.makedirs(d)
    for f in ("parameters", "seeds"):
        shutil.copy(os.path.join(start, f), os.path.join(d, f))
    edit_parameters(os.path.join(d, "parameters"), edits)
    for k in range(groups):
        os.makedirs(os.path.join(d, "Temp_%d" % k))
        shutil.copy(os.path.join(d, "parameters"), os.path.join(d, "Temp_%d" % k, "parameters"))
        edit_parameters(os.path.join(d, "Temp_%d" % k, "parameters"), {"ham_U": "%.1fd0" % (4.0 - 0.5 * k)})
        if groups:
            edit_parameters(os.path.join(d, "parameters"), {"mpi_per_parameter_set": "1"})

def subdirs(d, groups):
    return [os.path.join(d, "Temp_%d" % k) for k in range(groups)] if groups else [d]

def run(exe, d, env_extra, mpi, timeout):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    for k in [k for k in env if k.startswith("ALF_GPU_") or k.startswith("ALF_LA_TEST_")]:
        del env[k]
    env.update(env_extra)
    cmd = (["mpirun", "--oversubscribe", "-np", str(mpi)] if mpi else []) + [exe]
    t0 = time.time()
    with open(os.path.join(d, "stdout"), "w") as out:
        try:
            rc = subprocess.call(cmd, cwd=d, env=env, stdout=out, stderr=subprocess.STDOUT, timeout=timeout)
        except subprocess.TimeoutExpired:
            rc = -999
    return rc, time.time() - t0

def gpu_lines(d):
    try:
        t = open(os.path.join(d, "info"), errors="replace").read()
    except OSError:
        return []
    return [l.strip() for l in t.splitlines() if l.strip().startswith("GPU path") or l.strip().startswith("GPU observables") or l.strip().startswith("GPU tau pass")]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True); ap.add_argument("--start", required=True); ap.add_argument("--work", required=True)
    ap.add_argument("--mpi", type=int, default=0); ap.add_argument("--build", default="serial", help="label of the build in the first column"); ap.add_argument("--groups", type=int, default=0, help="PARALLEL_PARAMS/TEMPERING layout: Temp_0..Temp_<n-1> with ham_U = 4 - 0.5 k"); ap.add_argument("--restart", action="store_true", help="after the first run copy confout_* to confin_* and run again"); ap.add_argument("--only", default=""); ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument("cases")
    a = ap.parse_args()
    spec = importlib.util.spec_from_file_location("cases", a.cases); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    only = [x for x in a.only.split(",") if x]
    os.makedirs(a.work, exist_ok=True)
    res = open(os.path.join(a.work, "results.tsv"), "a")
    for label, edits, expect, env, note in m.CASES:
        if only and label not in only:
            continue
        da, dd = os.path.join(a.work, label, "active"), os.path.join(a.work, label, "disabled")
        prepare(da, a.start, edits, a.groups); prepare(dd, a.start, edits, a.groups)
        rca, ta = run(a.exe, da, env, a.mpi, a.timeout)
        rcd, td = run(a.exe, dd, dict(env, ALF_GPU_DISABLE="1"), a.mpi, a.timeout)
        if a.restart and rca == 0 and rcd == 0:
            for d in (da, dd):
                for f in os.listdir(d):
                    if f.startswith("confout_"):
                        shutil.copy(os.path.join(d, f), os.path.join(d, "confin_" + f[8:]))
            rca, ta2 = run(a.exe, da, env, a.mpi, a.timeout)
            rcd, td2 = run(a.exe, dd, dict(env, ALF_GPU_DISABLE="1"), a.mpi, a.timeout)
        lines = gpu_lines(subdirs(da, a.groups)[0])
        state = "ACTIVE" if any("gpu_state): ACTIVE" in l or "ACTIVE (mirror" in l for l in lines) else ("INACTIVE" if any("INACTIVE" in l for l in lines) else "NOINFO")
        reason = next((l.split("reason:", 1)[1].strip() for l in lines if "INACTIVE, reason:" in l), "")
        detail = ""
        if rca != 0 or rcd != 0:
            verdict = "stops" if (rcd != 0 and rca != 0) else "FAIL"
            detail = "exit codes active %d, disabled %d" % (rca, rcd)
            tail = open(os.path.join(da, "stdout"), errors="replace").read().strip().splitlines()[-3:]
            detail += "; " + " | ".join(x.strip() for x in tail)[:300]
        else:
            bad, diffs, tol = [], {}, 1e-11
            for sa, sd in zip(subdirs(da, a.groups), subdirs(dd, a.groups)):
                gp = e2e.green_precision(sd)
                tol = max(tol, 100 * gp[0]) if gp else tol
                b, df = compare(sa, sd, state != "ACTIVE", tol)
                bad += b; diffs.update({os.path.basename(sa) + "/" + k: v for k, v in df.items()})
            mx = max([v[1] for v in diffs.values() if v] + [0.0])
            if state == "ACTIVE":
                verdict = "works" if not bad else "FAIL"
            elif state == "INACTIVE":
                verdict = "refused" if not bad else "FAIL"
            else:
                verdict = "FAIL"
            detail = "max scaled diff %.2e (tol %.1e), %d files%s" % (mx, tol, len(diffs), ("; BAD: " + ", ".join(bad)) if bad else "") + ("; .h5 via h5diff (exact for confout/refused, -d tol otherwise)" if any(k.endswith(".h5") for k in diffs) else "")
        ok = (verdict == expect) or (expect.startswith("refused:") and verdict == "refused" and expect[8:] in reason)
        line = "\t".join([a.build, label, expect, verdict, state, "OK" if ok else "MISMATCH", reason, detail, "%.1fs/%.1fs" % (ta, td), note])
        print(line, flush=True); res.write(line + "\n"); res.flush()

if __name__ == "__main__":
    main()
