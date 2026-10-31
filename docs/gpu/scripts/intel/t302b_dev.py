#!/usr/bin/env python3
"""T0 development driver (ticket 302b, run on intel with the venv python): prepare dev-seed configs, build, run chosen runs of
302b.json from a dev candidate root made by t302b_sync.sh, and judge them with the shipped verdict.py (partial run roots give
INPUT_ERROR for gates whose runs are missing: read the gates you ran).
usage: t302b_dev.py TAG prep | build NAME.. | run ID.. | verdict | ctest ID.. """
import json, os, subprocess, sys, time, hashlib
sys.dont_write_bytecode = True
TAG, cmd, args = sys.argv[1], sys.argv[2], sys.argv[3:]
C = os.path.expanduser(f"~/alf-gpu/work/t302b/dev/cand/dev/{TAG}")
SRC = os.path.join(C, "src"); ACC = os.path.join(SRC, "testsuite/gpu_reference/acceptance")
os.environ["TMPDIR"] = os.path.expanduser("~/alf-gpu/work/t302b/dev/tmp")
ACCN = os.environ.get("ACCFILE", "302b.json")   # 302b.nvidia.json on nividia
acc = json.load(open(os.path.join(ACC, ACCN)))
cand = dict(l.split(None, 1) for l in open(os.path.join(C, ".alf_candidate")).read().splitlines() if " " in l)
sys.path.insert(0, ACC)
LOCK = os.path.expanduser("~/alf-gpu/gpu.lock")
cdir = os.path.join(C, "configs"); rroot = os.path.join(C, "runs"); os.makedirs(rroot, exist_ok=True); os.makedirs(os.path.join(C, "logs"), exist_ok=True)
def busy():
    if os.path.exists("/usr/bin/nvidia-smi") or subprocess.call("command -v nvidia-smi >/dev/null", shell=True) == 0:
        p = subprocess.run(["nvidia-smi", "--query-compute-apps=pid", "--format=csv,noheader"], capture_output=True, text=True)
        return [l for l in p.stdout.splitlines() if l.strip()]
    p = subprocess.run(["ps", "-eo", "user=,pcpu="], capture_output=True, text=True)
    return [l for l in p.stdout.splitlines() if l.split()[0] not in (os.environ["USER"], "root") and float(l.split()[1]) > 10]
if cmd == "prep":
    import acc_driver as ad
    ad.prep_configs(acc, SRC, cdir, "agent", sys.stdout)
elif cmd == "build":
    for n in args:
        b = acc["builds"][n]
        env = dict(os.environ, ALF_FLAGS_EXT=b.get("flags_ext", ""))
        t0 = time.time()
        p = subprocess.run(["bash", os.path.join(C, "tools/acc_build.sh"), "302b", "dev", TAG, n, b["tool"]] + b["configure"].split(),
                           env=env, capture_output=True, text=True)
        print(f"build {n}: rc {p.returncode} {p.stdout.strip()[-200:]} {time.time()-t0:.0f}s", flush=True)
elif cmd == "run":
    renv = dict(os.environ, ACC_TIMEOUT=str(acc.get("timeout_s", 900)), ACC_GPU_LOCK=LOCK)
    pool = []
    def launch(rid, alias=None):   # ID@ALIAS: a second run of declared run ID in directory ALIAS (not judged by verdict.py)
        r = acc["runs"][rid]
        s = r.get("seeds", "dev"); s = "dev" + s[3:] if s.startswith("acc") else s
        c = ["bash", os.path.join(ACC, "acc_run.sh"), os.path.join(rroot, alias or rid), os.path.join(C, "bin", r["build"] + ".out"), os.path.join(cdir, s, r["config"])]
        if r.get("allow_debug"): c.append("--allow-debug")
        return subprocess.Popen(c + r.get("env", []), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=renv)
    for arg in args:
        rid, _, alias = arg.partition("@"); alias = alias or None
        r = acc["runs"][rid]
        if os.path.exists(os.path.join(rroot, alias or rid)): print(f"{arg}: exists, kept"); continue
        if acc["builds"][r["build"]].get("gpu", True):
            if busy(): print("STOP: other users busy", busy()); sys.exit(1)
            t0 = time.time(); p = launch(rid, alias); out = p.communicate()[0]
            print(f"gpu run {rid}: rc {p.returncode} wall {time.time()-t0:.0f}s {out.strip()[-150:]}", flush=True)
        else:
            while len([x for x in pool if x[1].poll() is None]) >= 2: time.sleep(1)
            pool.append((arg, launch(rid, alias), time.time())); print(f"cpu run {rid}: started", flush=True)
    for rid, p, t0 in pool:
        p.wait(); print(f"cpu run {rid}: rc {p.returncode} wall {time.time()-t0:.0f}s", flush=True)
elif cmd == "ctest":
    for tid in args:
        t = next(x for x in acc["ctests"] if x["id"] == tid)
        c = ["bash", os.path.join(C, "tools/acc_ctest.sh"), "302b", "dev", TAG, os.path.join(rroot, tid), t["tool"]]
        if t.get("regex"): c += ["--regex", t["regex"]]
        for kv in t.get("env", []): c += ["--env", kv]
        t0 = time.time(); q = subprocess.run(c + ["--"] + t["configure"].split(), capture_output=True, text=True)
        print(f"ctest {tid}: rc {q.returncode} {q.stdout.strip()[-200:]} {time.time()-t0:.0f}s", flush=True)
elif cmd == "verdict":
    vd = os.path.join(C, "verdict", time.strftime("%H%M%S")); os.makedirs(vd)
    subprocess.run([sys.executable, os.path.join(ACC, "verdict.py"), os.path.join(ACC, ACCN), rroot, "--out", vd,
                    "--candidate", cand["sha"], "--seeds-record", os.path.join(cdir, "seeds_record.json")])
    print("verdict dir", vd)
