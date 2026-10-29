#!/usr/bin/env python3
"""Acceptance driver (machine side): prepare configs, build, run and judge one frozen candidate, unattended.

usage: acc_driver.py --ticket NNN --sha12 SHA12 --who mgr|agent [--only RUN_ID,...] [--no-verdict]

Needs a candidate shipped by acc_ship.sh (~/alf-gpu/work/t<NNN>/cand/<sha12>/{src,base_src,.alf_candidate}) and the env
files in ~/alf-gpu/env/. The acceptance file and every tool come from the shipped, read-only candidate tree, which
acc_ship.sh has checked against the pre-registration commit; the driver re-checks the file's sha256 against
.alf_candidate. Layout: cand/<sha12>/{bin,tree_*,logs, configs/<who>, runs/<who>/<run_id>, verdict/<who>}.
Rules:
  * GPU runs one at a time, under flock ~/alf-gpu/gpu.lock (also against a second driver), after an idle check that
    must succeed and show no process of anybody; otherwise the driver stops;
  * CPU-only runs (builds with "gpu": false) run in a pool of at most cpu_parallel processes alongside;
  * seeds: --who mgr draws fresh acceptance seeds ("acc*" sets) with secrets.randbits(31); --who agent (tier T1) maps
    every "acc*" set to "dev" and never sees acceptance seeds; configs and seeds are separate per party;
  * an existing binary is reused only if its build record matches the acceptance file;
  * "mutations" (script-level positive controls) are applied after the runs with mutate.py, which must succeed;
  * "gpunone": true runs acc_gpunone.sh into runs/<who>/gpunone;
  * free disk >= disk_min_free_GB at start; run directories are never overwritten; TMPDIR stays in the ticket dir.
Exit code: verdict.py's exit code (0 = all gates PASS and all controls as expected); 1 on a driver failure.
"""
import argparse
import hashlib
import json
import os
import secrets
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
LOCK = os.path.expanduser("~/alf-gpu/gpu.lock")


def log(msg, fh):
    line = time.strftime("%Y-%m-%dT%H:%M:%S ") + msg
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


def gpu_idle():
    """(ok, detail): ok only if nvidia-smi succeeds and lists no compute process; machines without nvidia-smi
    (intel) report ok with a note, the run's provenance records it."""
    if shutil.which("nvidia-smi") is None:
        return True, "no nvidia-smi"
    p = subprocess.run(["nvidia-smi", "--query-compute-apps=pid,process_name,used_memory", "--format=csv,noheader"],
                       capture_output=True, text=True)
    apps = [l for l in p.stdout.splitlines() if l.strip()]
    return p.returncode == 0 and not apps, f"rc {p.returncode} {apps}"


def seed_set(r, who):
    s = r.get("seeds", "dev")
    return "dev" if (who == "agent" and s.startswith("acc")) else s


def prep_configs(acc, src, cdir, who, fh):
    needed = {}
    for r in acc["runs"].values():
        needed.setdefault(seed_set(r, who), set()).add(r["config"])
    code = r"""
import os, sys, json
sys.path.insert(0, sys.argv[2])
from py_alf import ALF_source, Simulation
from configs import CONFIGS
try:
    from configs_stat import CONFIGS as CS
    CONFIGS = {**CONFIGS, **CS}
except ImportError:
    pass
out, alf, spec = sys.argv[1], sys.argv[3], json.loads(sys.argv[4])
src = ALF_source(alf_dir=alf)
for name, c in spec.items():
    d = {**CONFIGS[c["base"]], **c.get("overrides", {})}
    Simulation(src, "Hubbard", d, sim_root=out, sim_dir=name, machine="GNU", mpi=True).run(only_prep=True)
"""
    record = {}
    for s in sorted(needed):
        out = os.path.join(cdir, s)
        if os.path.exists(out):
            raise SystemExit(f"refusing: {out} exists")
        if s != "dev" and not s.startswith("acc"):
            raise SystemExit(f"unknown seed set {s}")
        spec = {n: acc["configs"][n] for n in sorted(needed[s])}
        p = subprocess.run([sys.executable, "-I", "-c", code, out, os.path.join(src, "testsuite", "gpu_reference"), src,
                            json.dumps(spec)], capture_output=True, text=True)
        if p.returncode != 0:
            raise SystemExit(f"prep failed for seed set {s}: {p.stderr[-500:]}")
        for n in spec:
            sf = os.path.join(out, n, "seeds")
            if s.startswith("acc"):  # fresh, never-seen seeds (manager only)
                lines = open(sf).read().splitlines()
                open(sf, "w").write("\n".join(f"{secrets.randbits(31) % 999983 + 1:12d}" for _ in lines) + "\n")
            record[f"{s}/{n}"] = hashlib.sha256(open(sf, "rb").read()).hexdigest()
            ps = acc["configs"][n].get("parameters_sha256")
            if ps and hashlib.sha256(open(os.path.join(out, n, "parameters"), "rb").read()).hexdigest() != ps:
                raise SystemExit(f"config {n}: parameters differ from the pre-registered sha256")
        log(f"prepared seed set {s}: {sorted(spec)}", fh)
    if who == "agent":  # the T1 run uses dev seeds under the acc* names of the acceptance file
        for r in acc["runs"].values():
            s = r.get("seeds", "dev")
            if s.startswith("acc"):
                record[f"{s}/{r['config']}"] = record[f"dev/{r['config']}"]
    json.dump(record, open(os.path.join(cdir, "seeds_record.json"), "w"), indent=1)
    return record


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ticket", required=True)
    ap.add_argument("--sha12", required=True)
    ap.add_argument("--who", choices=["mgr", "agent"], required=True)
    ap.add_argument("--only", default="")
    ap.add_argument("--no-verdict", action="store_true")
    a = ap.parse_args()
    C = os.path.expanduser(f"~/alf-gpu/work/t{a.ticket}/cand/{a.sha12}")
    os.environ["TMPDIR"] = os.path.expanduser(f"~/alf-gpu/work/t{a.ticket}/tmp")
    os.makedirs(os.environ["TMPDIR"], exist_ok=True)
    src = os.path.join(C, "src")
    accf = os.path.join(src, "testsuite", "gpu_reference", "acceptance", f"{a.ticket}.json")
    raw = open(accf, "rb").read()
    cand = dict(l.split(None, 1) for l in open(os.path.join(C, ".alf_candidate")).read().splitlines() if " " in l)
    if hashlib.sha256(raw).hexdigest() != cand.get("acceptance_sha256"):
        raise SystemExit("the acceptance file differs from the pre-registered one recorded by acc_ship.sh")
    if os.path.abspath(HERE) != os.path.abspath(os.path.join(src, "testsuite", "gpu_reference", "acceptance")):
        raise SystemExit("run the driver from the shipped candidate tree (its tools were checked against the pre-registration)")
    acc = json.loads(raw)
    os.makedirs(os.path.join(C, "logs"), exist_ok=True)
    fh = open(os.path.join(C, "logs", f"driver_{a.who}_{time.strftime('%Y%m%d_%H%M%S')}.log"), "w")
    free = shutil.disk_usage(C).free / 2**30
    if free < acc.get("disk_min_free_GB", 5):
        raise SystemExit(f"only {free:.1f} GB free, need {acc.get('disk_min_free_GB', 5)}")
    log(f"candidate {cand['sha']} prereg {cand.get('prereg_sha')} acceptance sha256 {cand['acceptance_sha256']} "
        f"who {a.who} free {free:.1f} GB", fh)
    cdir = os.path.join(C, "configs", a.who)
    rec_path = os.path.join(cdir, "seeds_record.json")
    if not os.path.exists(rec_path):
        prep_configs(acc, src, cdir, a.who, fh)
    # builds
    for name, b in acc["builds"].items():
        bj = os.path.join(C, "bin", name + ".build.json")
        if os.path.exists(os.path.join(C, "bin", name + ".out")):
            rec = json.load(open(bj))
            if rec["configure"].split() != b["configure"].split() or rec["ALF_FLAGS_EXT"] != b.get("flags_ext", "") \
                    or rec["tool"] != b["tool"]:
                raise SystemExit(f"existing binary {name} was built differently from the acceptance file")
            log(f"build {name}: exists, build record matches, reused", fh)
            continue
        env = dict(os.environ, ALF_FLAGS_EXT=b.get("flags_ext", ""))
        p = subprocess.run(["bash", os.path.join(HERE, "acc_build.sh"), a.ticket, a.sha12, name, b["tool"]]
                           + b["configure"].split(), env=env, capture_output=True, text=True)
        log(f"build {name}: rc {p.returncode} {p.stdout.strip()[-200:]}", fh)
        if p.returncode != 0:
            raise SystemExit(f"build {name} failed")
    # runs
    only = [x for x in a.only.split(",") if x]
    rroot = os.path.join(C, "runs", a.who)
    os.makedirs(rroot, exist_ok=True)
    pool = []
    renv = dict(os.environ, ACC_TIMEOUT=str(acc.get("timeout_s", 900)), ACC_GPU_LOCK=LOCK)

    def launch(rid):
        r = acc["runs"][rid]
        cmd = ["bash", os.path.join(HERE, "acc_run.sh"), os.path.join(rroot, rid), os.path.join(C, "bin", r["build"] + ".out"),
               os.path.join(cdir, seed_set(r, a.who), r["config"])]
        if r.get("allow_debug"):
            cmd.append("--allow-debug")
        return subprocess.Popen(cmd + r.get("env", []), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=renv)

    def reap(block=False):
        for item in list(pool):
            rid, p = item
            if block:
                p.wait()
            if p.poll() is not None:
                pool.remove(item)
                log(f"cpu run {rid}: rc {p.returncode}", fh)

    mut_dst = {m["dst"] for m in acc.get("mutations", [])}  # declared for provenance, produced by mutate.py, never run
    for rid in [k for k in acc["runs"] if (not only or k in only) and k not in mut_dst]:
        r = acc["runs"][rid]
        if os.path.exists(os.path.join(rroot, rid)):
            log(f"run {rid}: directory exists, kept as is", fh)
            continue
        if acc["builds"][r["build"]].get("gpu", True):
            ok, detail = gpu_idle()
            if not ok:
                log(f"STOP: GPU not idle before {rid}: {detail}", fh)
                reap(block=True)
                return 1
            t0 = time.time()
            p = launch(rid)
            p.wait()
            log(f"gpu run {rid}: rc {p.returncode} wall {time.time() - t0:.0f} s", fh)
        else:
            if r.get("refcache"):
                if r["build"].startswith("base_"):
                    raise SystemExit(f"run {rid}: refcache is not allowed for base builds")
                q = subprocess.run([sys.executable, os.path.join(HERE, "refcache.py"), "get", os.path.join(rroot, rid), src,
                                    os.path.join(C, "bin", r["build"] + ".build.json"),
                                    os.path.join(cdir, seed_set(r, a.who), r["config"])], capture_output=True, text=True)
                log(f"cpu run {rid}: refcache {q.stdout.strip()} {q.stderr.strip()[-200:]}", fh)
                if q.returncode == 0:
                    continue
            while len([x for x in pool if x[1].poll() is None]) >= acc.get("cpu_parallel", 3):
                time.sleep(2)
                reap()
            pool.append((rid, launch(rid)))
            log(f"cpu run {rid}: started", fh)
        reap()
    reap(block=True)
    for m in acc.get("mutations", []):
        dst = os.path.join(rroot, m["dst"])
        if os.path.exists(dst):
            continue
        q = subprocess.run([sys.executable, os.path.join(HERE, "mutate.py"), os.path.join(rroot, m["src"]), dst, m["kind"]],
                           capture_output=True, text=True)
        log(f"mutation {m['dst']} ({m['kind']}): rc {q.returncode} {q.stdout.strip()[-200:]}", fh)
        if q.returncode != 0:
            raise SystemExit(f"mutation {m['dst']} failed: the control would test nothing")
    for sub, extra in (("gpunone", []), ("gpunone_plant", ["--plant"])):
        if acc.get("gpunone") and not os.path.exists(os.path.join(rroot, sub)):
            q = subprocess.run(["bash", os.path.join(HERE, "acc_gpunone.sh"), a.ticket, a.sha12, os.path.join(rroot, sub)] + extra,
                               capture_output=True, text=True)
            log(f"{sub}: rc {q.returncode} {q.stdout.strip()[-300:]}", fh)
    if a.no_verdict:
        return 0
    vdir = os.path.join(C, "verdict", a.who)
    os.makedirs(vdir, exist_ok=True)
    p = subprocess.run([sys.executable, os.path.join(HERE, "verdict.py"), accf, rroot, "--out", vdir,
                        "--candidate", cand["sha"], "--seeds-record", rec_path])
    log(f"verdict exit {p.returncode}", fh)
    return p.returncode


if __name__ == "__main__":
    sys.exit(main())
