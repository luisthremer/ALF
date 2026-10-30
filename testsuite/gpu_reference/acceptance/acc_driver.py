#!/usr/bin/env python3
"""Acceptance driver (machine side): prepare configs, build, run and judge one frozen candidate, unattended.

usage: acc_driver.py --ticket NNN --sha12 SHA12 --who mgr|agent [--part NAME] [--reason TEXT] [--only RUN_ID,...]
                     [--no-verdict]

Needs a candidate shipped by acc_ship.sh (~/alf-gpu/work/t<NNN>/cand/<who>/<sha12>/{src,base_src,.alf_candidate}) and the env
files in ~/alf-gpu/env/. The acceptance file and every tool come from the shipped, read-only candidate tree, which
acc_ship.sh has checked against the pre-registration commit; the driver re-checks the file's sha256 against
.alf_candidate; the file must pass verdict.py's lint. Layout: cand/<who>/<sha12>/{bin,tree_*,logs,configs,runs/<run_id>,
verdict}; ~/alf-gpu/archive/attempts/t<NNN>.log records every attempt with its outcome (a second manager T2 needs --reason).
Rules:
  * GPU runs one at a time, under flock ~/alf-gpu/gpu.lock (also against a second driver), after an idle check that
    must succeed and show no process of anybody; otherwise the driver stops;
  * CPU-only runs (builds with "gpu": false) run in a pool of at most cpu_parallel processes alongside;
  * seeds: --who mgr draws fresh acceptance seeds ("acc*" sets) with secrets.randbits(31); --who agent (tier T1) maps
    every "accN" set to "devN" = the committed dev seeds rotated by 37*N lines (distinct per N, never acceptance seeds);
  * on a machine without nvidia-smi (intel) the GPU idle check is replaced by a user-load check (no other user's process
    above 10 % CPU);
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


def other_users_busy():
    p = subprocess.run(["ps", "-eo", "user=,pcpu="], capture_output=True, text=True)
    me = os.environ.get("USER", "")
    return [l for l in p.stdout.splitlines() if l.split()[0] not in (me, "root") and float(l.split()[1]) > 10]


def gpu_idle():
    """(ok, detail): ok only if nvidia-smi succeeds and lists no compute process; on machines without nvidia-smi
    (intel) only if no other user's process uses more than 10 % CPU."""
    if shutil.which("nvidia-smi") is None:
        busy = other_users_busy()
        return not busy, f"no nvidia-smi; other users busy: {busy}"
    p = subprocess.run(["nvidia-smi", "--query-compute-apps=pid,process_name,used_memory", "--format=csv,noheader"],
                       capture_output=True, text=True)
    apps = [l for l in p.stdout.splitlines() if l.strip()]
    return p.returncode == 0 and not apps, f"rc {p.returncode} {apps}"


def seed_set(r, who):
    s = r.get("seeds", "dev")
    return "dev" + s[3:] if (who == "agent" and s.startswith("acc")) else s


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
        if s != "dev" and not s.startswith(("acc", "dev")):
            raise SystemExit(f"unknown seed set {s}")
        spec = {n: acc["configs"][n] for n in sorted(needed[s])}
        p = subprocess.run([sys.executable, "-I", "-c", code, out, os.path.join(src, "testsuite", "gpu_reference"), src,
                            json.dumps(spec)], capture_output=True, text=True)
        if p.returncode != 0:
            raise SystemExit(f"prep failed for seed set {s}: {p.stderr[-500:]}")
        for n in spec:
            sf = os.path.join(out, n, "seeds")
            lines = open(sf).read().splitlines()
            if s.startswith("acc"):  # fresh, never-seen seeds (manager only)
                open(sf, "w").write("\n".join(f"{secrets.randbits(31) % 999983 + 1:12d}" for _ in lines) + "\n")
            elif s != "dev":  # devN: the committed dev seeds rotated by 37*N lines (deterministic, distinct per N)
                k = (37 * int(s[3:] or 1)) % len(lines)
                open(sf, "w").write("\n".join(lines[k:] + lines[:k]) + "\n")
            record[f"{s}/{n}"] = hashlib.sha256(open(sf, "rb").read()).hexdigest()
            ps = acc["configs"][n].get("parameters_sha256")
            if ps and hashlib.sha256(open(os.path.join(out, n, "parameters"), "rb").read()).hexdigest() != ps:
                raise SystemExit(f"config {n}: parameters differ from the pre-registered sha256")
        log(f"prepared seed set {s}: {sorted(spec)}", fh)
    if who == "agent":  # the T1 run uses devN seeds under the accN names of the acceptance file
        for r in acc["runs"].values():
            s = r.get("seeds", "dev")
            if s.startswith("acc"):
                record[f"{s}/{r['config']}"] = record[f"dev{s[3:]}/{r['config']}"]
    os.makedirs(cdir, exist_ok=True)  # an acceptance file with ctests only has no runs and no prepared configs
    json.dump(record, open(os.path.join(cdir, "seeds_record.json"), "w"), indent=1)
    return record


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ticket", required=True)
    ap.add_argument("--sha12", required=True)
    ap.add_argument("--who", choices=["mgr", "agent"], required=True)
    ap.add_argument("--only", default="")
    ap.add_argument("--part", default="", help="acceptance part: <ticket>.<part>.json")
    ap.add_argument("--reason", default="", help="required for a repeated manager T2 (recorded)")
    ap.add_argument("--no-verdict", action="store_true")
    a = ap.parse_args()
    C = os.path.expanduser(f"~/alf-gpu/work/t{a.ticket}/cand/{a.who}/{a.sha12}")
    os.environ["TMPDIR"] = os.path.expanduser(f"~/alf-gpu/work/t{a.ticket}/tmp")
    os.makedirs(os.environ["TMPDIR"], exist_ok=True)
    src = os.path.join(C, "src")
    name = f"{a.ticket}.{a.part}.json" if a.part else f"{a.ticket}.json"
    accf = os.path.join(src, "testsuite", "gpu_reference", "acceptance", name)
    raw = open(accf, "rb").read()
    lines = open(os.path.join(C, ".alf_candidate")).read().splitlines()
    cand = dict(l.split(None, 1) for l in lines if " " in l)
    parts = {l.split()[1]: l.split()[2] for l in lines if l.startswith("acceptance_part ")}
    if cand.get("who") != a.who:
        raise SystemExit(f"candidate root of {cand.get('who')}, not {a.who}")
    if hashlib.sha256(raw).hexdigest() != parts.get(name):
        raise SystemExit(f"{name} differs from the pre-registered file recorded by acc_ship.sh")
    if os.path.abspath(HERE) != os.path.abspath(os.path.join(src, "testsuite", "gpu_reference", "acceptance")):
        raise SystemExit("run the driver from the shipped candidate tree (its tools were checked against the pre-registration)")
    acc = json.loads(raw)
    dec = os.path.join(src, "docs", "gpu", "DECISIONS.md")
    if hashlib.sha256(open(dec, "rb").read()).hexdigest() != cand.get("decisions_sha256"):
        raise SystemExit("docs/gpu/DECISIONS.md differs from its state at the pre-registration (report_only approvals must predate it)")
    if not cand.get("base_sha", "").startswith(acc.get("base_sha", "-")) or len(acc.get("base_sha", "")) < 12:
        raise SystemExit(f"shipped base {cand.get('base_sha')} is not the acceptance file's base_sha {acc.get('base_sha')}")
    sys.path.insert(0, HERE)
    import verdict  # noqa: E402  (the shipped, pinned copy)
    errs = verdict.lint_errors(acc)
    if errs:
        raise SystemExit("the acceptance file fails lint: " + "; ".join(errs[:5]))
    if acc.get("machine") not in ("nividia", "intel"):
        raise SystemExit(f"unknown machine {acc.get('machine')}")
    if (acc["machine"] == "intel") != (shutil.which("nvidia-smi") is None):
        raise SystemExit(f"this acceptance file is for {acc['machine']}; run it there")
    attempts = os.path.expanduser(f"~/alf-gpu/archive/attempts/t{a.ticket}.log")  # manager-owned directory
    os.makedirs(os.path.dirname(attempts), exist_ok=True)
    earlier = [l for l in open(attempts)] if os.path.exists(attempts) else []
    if a.who == "mgr" and any(" mgr " in l and name in l for l in earlier) and not a.reason:
        raise SystemExit("a manager T2 of this ticket exists already: give --reason (re-running until green is not allowed)")
    with open(attempts, "a") as f:
        f.write(f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {a.who} {cand['sha']} {name} reason={a.reason or '-'}\n")
    os.makedirs(os.path.join(C, "logs"), exist_ok=True)
    fh = open(os.path.join(C, "logs", f"driver_{a.who}_{time.strftime('%Y%m%d_%H%M%S')}.log"), "w")
    free = shutil.disk_usage(C).free / 2**30
    if free < acc.get("disk_min_free_GB", 5):
        raise SystemExit(f"only {free:.1f} GB free, need {acc.get('disk_min_free_GB', 5)}")
    log(f"candidate {cand['sha']} prereg {cand.get('prereg_sha')} acceptance sha256 {cand['acceptance_sha256']} "
        f"who {a.who} free {free:.1f} GB", fh)
    cdir = os.path.join(C, "configs")
    rec_path = os.path.join(cdir, "seeds_record.json")
    if not os.path.exists(rec_path):
        prep_configs(acc, src, cdir, a.who, fh)
    # builds
    for bname, b in acc["builds"].items():
        bj = os.path.join(C, "bin", bname + ".build.json")
        if os.path.exists(os.path.join(C, "bin", bname + ".out")):
            rec = json.load(open(bj))
            if rec["configure"].split() != b["configure"].split() or rec["ALF_FLAGS_EXT"] != b.get("flags_ext", "") \
                    or rec["tool"] != b["tool"]:
                raise SystemExit(f"existing binary {bname} was built differently from the acceptance file")
            log(f"build {bname}: exists in this party's root, build record matches, reused", fh)
            continue
        env = dict(os.environ, ALF_FLAGS_EXT=b.get("flags_ext", ""))
        p = subprocess.run(["bash", os.path.join(HERE, "acc_build.sh"), a.ticket, a.who, a.sha12, bname, b["tool"]]
                           + b["configure"].split(), env=env, capture_output=True, text=True)
        log(f"build {bname}: rc {p.returncode} {p.stdout.strip()[-200:]}", fh)
        if p.returncode != 0:
            raise SystemExit(f"build {bname} failed")
    # runs
    only = [x for x in a.only.split(",") if x]
    rroot = os.path.join(C, "runs")
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
    for t in acc.get("ctests", []):  # unit-test gates and their controls (acc_ctest.sh), each in runs/<id>
        if os.path.exists(os.path.join(rroot, t["id"])):
            continue
        env_items = t.get("env", [])
        if any(kv.split("=")[0].startswith(("ALF_GPU_", "ALFSH_")) for kv in env_items) and not t.get("control"):
            raise SystemExit(f"ctest {t['id']}: debug/fault variables only in a declared control")
        if acc["machine"] != "intel" or t.get("gpu"):
            ok, detail = gpu_idle()
            if not ok:
                log(f"STOP: not idle before ctest {t['id']}: {detail}", fh)
                return 1
        cmd = ["bash", os.path.join(HERE, "acc_ctest.sh"), a.ticket, a.who, a.sha12, os.path.join(rroot, t["id"]), t["tool"]]
        if t.get("base"):
            cmd += ["--base"]
        if t.get("regex"):
            cmd += ["--regex", t["regex"]]
        if t.get("repeat"):
            cmd += ["--repeat", str(t["repeat"])]
        for kv in env_items:
            cmd += ["--env", kv]
        q = subprocess.run(cmd + ["--"] + t["configure"].split(), capture_output=True, text=True)
        log(f"ctest {t['id']}: rc {q.returncode} {q.stdout.strip()[-200:]}", fh)
    # GPU=none before the log mutations: a log mutation may copy gpunone (e.g. gpunone_m_ctestfail)
    for sub, extra in (("gpunone", []), ("gpunone_plant", ["--plant"])):
        if acc.get("gpunone") and not os.path.exists(os.path.join(rroot, sub)):
            q = subprocess.run(["bash", os.path.join(HERE, "acc_gpunone.sh"), a.ticket, a.who, a.sha12, os.path.join(rroot, sub)] + extra,
                               capture_output=True, text=True)
            log(f"{sub}: rc {q.returncode} {q.stdout.strip()[-300:]}", fh)
    for m in acc.get("log_mutations", []):  # script-level controls of log gates: a copy with one planted failure
        src_d, dst_d = os.path.join(rroot, m["src"]), os.path.join(rroot, m["dst"])
        if os.path.exists(dst_d) or not os.path.isdir(src_d):
            if not os.path.isdir(src_d):
                raise SystemExit(f"log mutation {m['dst']}: source {m['src']} missing")
            continue
        shutil.copytree(src_d, dst_d)
        f = os.path.join(dst_d, m["file"])
        t = open(f).read()
        if m["kind"] == "ctest_one_failure":
            import re as _re
            new = _re.sub(r"(\d+)% tests passed, 0 tests failed out of (\d+)", lambda x: f"99% tests passed, 1 tests failed out of {x.group(2)}", t)
            # the exit file stays as it was: the control must be killed by the gate's regex/equals criterion
        elif m["kind"] == "change_first_number":
            import re as _re
            new = _re.sub(r"(?<![\w.])(\d+\.\d{3,})", lambda x: x.group(1)[:-1] + str((int(x.group(1)[-1]) + 1) % 10), t, count=1)
        else:
            raise SystemExit(f"unknown log mutation {m['kind']}")
        if new == t:
            raise SystemExit(f"log mutation {m['dst']} changed nothing: the control would test nothing")
        open(f, "w").write(new)
        open(os.path.join(dst_d, "MUTATED"), "w").write(m["kind"] + "\n")
        log(f"log mutation {m['dst']} ({m['kind']} in {m['file']})", fh)
    if a.no_verdict:
        return 0
    vdir = os.path.join(C, "verdict", a.part or "main")
    os.makedirs(vdir, exist_ok=True)
    p = subprocess.run([sys.executable, os.path.join(HERE, "verdict.py"), accf, rroot, "--out", vdir,
                        "--candidate", cand["sha"], "--seeds-record", rec_path, "--attempts", attempts])
    log(f"verdict exit {p.returncode}", fh)
    with open(attempts, "a") as f:
        f.write(f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {a.who} {cand['sha']} {name} outcome_exit={p.returncode}\n")
    return p.returncode


if __name__ == "__main__":
    sys.exit(main())
