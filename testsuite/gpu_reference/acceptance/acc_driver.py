#!/usr/bin/env python3
"""Acceptance driver (machine side): prepare configs, build, run and judge one frozen candidate, unattended.

usage: acc_driver.py ACCEPTANCE.json --ticket NNN --sha12 SHA12 --who mgr|agent [--only RUN_ID,...] [--no-verdict]

Needs a candidate shipped by acc_ship.sh (~/alf-gpu/work/t<NNN>/cand/<sha12>/{src,.alf_candidate}) and the env files in
~/alf-gpu/env/. Layout: cand/<sha12>/{bin,tree_*,configs,runs/<who>/<run_id>,logs}. Rules:
  * the acceptance file must be the one in the candidate's own tree (testsuite/gpu_reference/acceptance/NNN.json):
    its sha256 is compared, so it was frozen with the candidate;
  * GPU runs one at a time; before each, the GPU must be idle (also of other users' processes): otherwise stop;
  * CPU-only runs (builds with "gpu": false) run in a pool of at most cpu_parallel processes alongside;
  * fresh acceptance seeds ("acc*" seed sets) are drawn here with secrets.randbits(31), written to
    configs/seeds_acc.json, and never shown to the implementer before the candidate is frozen;
  * free disk >= disk_min_free_GB at start; run directories are never overwritten;
  * every run goes through acc_run.sh (env -i, whitelist, provenance.json, timeout).
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


def log(msg, fh):
    line = time.strftime("%Y-%m-%dT%H:%M:%S ") + msg
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


def gpu_apps():
    p = subprocess.run(["nvidia-smi", "--query-compute-apps=pid,process_name,used_memory", "--format=csv,noheader"],
                       capture_output=True, text=True)
    return [l for l in p.stdout.splitlines() if l.strip()] if p.returncode == 0 else ["nvidia-smi failed"]


def prep_configs(acc, src, cdir, fh):
    """pyALF prep of every config into configs/<seedset>/<name>; acc seed sets get fresh seeds."""
    sets = sorted({r.get("seeds", "dev") for r in acc["runs"].values()})
    needed = {}
    for r in acc["runs"].values():
        needed.setdefault(r.get("seeds", "dev"), set()).add(r["config"])
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
    for s in sets:
        out = os.path.join(cdir, s)
        if os.path.exists(out):
            raise SystemExit(f"refusing: {out} exists")
        spec = {n: acc["configs"][n] for n in sorted(needed[s])}
        p = subprocess.run([sys.executable, "-I", "-c", code, out, os.path.join(src, "testsuite", "gpu_reference"), src,
                            json.dumps(spec)], capture_output=True, text=True)
        if p.returncode != 0:
            raise SystemExit(f"prep failed for seed set {s}: {p.stderr[-500:]}")
        for n in spec:
            sf = os.path.join(out, n, "seeds")
            if s.startswith("acc"):  # fresh, never-seen seeds
                lines = open(sf).read().splitlines()
                new = [f"{secrets.randbits(31) % 999983 + 1:12d}" for _ in lines]
                open(sf, "w").write("\n".join(new) + "\n")
            elif s != "dev":
                raise SystemExit(f"unknown seed set {s}")
            record[f"{s}/{n}"] = hashlib.sha256(open(sf, "rb").read()).hexdigest()
        log(f"prepared seed set {s}: {sorted(spec)}", fh)
    json.dump(record, open(os.path.join(cdir, "seeds_acc.json"), "w"), indent=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("acceptance")
    ap.add_argument("--ticket", required=True)
    ap.add_argument("--sha12", required=True)
    ap.add_argument("--who", choices=["mgr", "agent"], required=True)
    ap.add_argument("--only", default="")
    ap.add_argument("--no-verdict", action="store_true")
    a = ap.parse_args()
    raw = open(a.acceptance, "rb").read()
    acc = json.loads(raw)
    C = os.path.expanduser(f"~/alf-gpu/work/t{a.ticket}/cand/{a.sha12}")
    src = os.path.join(C, "src")
    own = os.path.join(src, "testsuite", "gpu_reference", "acceptance", f"{a.ticket}.json")
    if not os.path.isfile(own) or open(own, "rb").read() != raw:
        raise SystemExit(f"the acceptance file is not the candidate's own frozen copy ({own})")
    cand = dict(l.split(None, 1) for l in open(os.path.join(C, ".alf_candidate")).read().splitlines() if " " in l)
    if not cand["sha"].startswith(acc["candidate_sha"]) and acc["candidate_sha"] != "self":
        raise SystemExit("candidate sha mismatch")
    os.makedirs(os.path.join(C, "logs"), exist_ok=True)
    fh = open(os.path.join(C, "logs", f"driver_{a.who}_{time.strftime('%Y%m%d_%H%M%S')}.log"), "w")
    free = shutil.disk_usage(C).free / 2**30
    if free < acc.get("disk_min_free_GB", 5):
        raise SystemExit(f"only {free:.1f} GB free, need {acc.get('disk_min_free_GB', 5)}")
    log(f"candidate {cand['sha']} acceptance sha256 {hashlib.sha256(raw).hexdigest()} who {a.who} free {free:.1f} GB", fh)
    cdir = os.path.join(C, "configs")
    if not os.path.exists(cdir):
        prep_configs(acc, src, cdir, fh)
    # builds
    for name, b in acc["builds"].items():
        if os.path.exists(os.path.join(C, "bin", name + ".out")):
            log(f"build {name}: exists, reused", fh)
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
    pool, failures = [], []
    order = [k for k in acc["runs"] if not only or k in only]

    def launch(rid):
        r = acc["runs"][rid]
        cmd = ["bash", os.path.join(HERE, "acc_run.sh"), os.path.join(rroot, rid), os.path.join(C, "bin", r["build"] + ".out"),
               os.path.join(cdir, r.get("seeds", "dev"), r["config"])]
        if r.get("allow_debug"):
            cmd.append("--allow-debug")
        return subprocess.Popen(cmd + r.get("env", []), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    def reap(block=False):
        for item in list(pool):
            rid, p = item
            if block:
                p.wait()
            if p.poll() is not None:
                pool.remove(item)
                log(f"cpu run {rid}: rc {p.returncode}", fh)

    for rid in order:
        r = acc["runs"][rid]
        if os.path.exists(os.path.join(rroot, rid)):
            log(f"run {rid}: directory exists, kept as is", fh)
            continue
        gpu = acc["builds"][r["build"]].get("gpu", True)
        if gpu:
            apps = gpu_apps()
            if apps:
                log(f"STOP: GPU not idle before {rid}: {apps}", fh)
                reap(block=True)
                return 1
            t0 = time.time()
            p = launch(rid)
            p.wait()
            log(f"gpu run {rid}: rc {p.returncode} wall {time.time() - t0:.0f} s", fh)
        else:
            if r.get("refcache"):
                q = subprocess.run([sys.executable, os.path.join(HERE, "refcache.py"), "get", os.path.join(rroot, rid), src,
                                    os.path.join(C, "bin", r["build"] + ".build.json"),
                                    os.path.join(cdir, r.get("seeds", "dev"), r["config"])], capture_output=True, text=True)
                log(f"cpu run {rid}: refcache {q.stdout.strip()}", fh)
                if q.returncode == 0:
                    continue
            while len([x for x in pool if x[1].poll() is None]) >= acc.get("cpu_parallel", 3):
                time.sleep(2)
                reap()
            pool.append((rid, launch(rid)))
            log(f"cpu run {rid}: started", fh)
        reap()
    reap(block=True)
    if a.no_verdict:
        return 0
    os.makedirs(os.path.join(C, "verdict", a.who), exist_ok=True)
    p = subprocess.run([sys.executable, os.path.join(HERE, "verdict.py"), a.acceptance, rroot,
                        "--out", os.path.join(C, "verdict", a.who), "--candidate", cand["sha"]])
    log(f"verdict exit {p.returncode}", fh)
    return p.returncode


if __name__ == "__main__":
    sys.exit(main())
