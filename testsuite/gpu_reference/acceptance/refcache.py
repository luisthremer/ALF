#!/usr/bin/env python3
"""Cache of CPU reference runs (GPU=none and test-only CPU variants), keyed by everything that determines their output.

usage: refcache.py key  RUN_DIR SRC_DIR            print the key of a finished CPU run
       refcache.py put  RUN_DIR SRC_DIR            store the run (manager only); refuses an existing key
       refcache.py get  RUN_DIR_NEW SRC_DIR BUILD_JSON CONFIG_DIR
                                                   copy a cached run whose key matches into RUN_DIR_NEW (exit 1 if none)
       refcache.py list | check                    list entries / verify every entry's sha256 manifest

Cache root: ~/alf-gpu/refcache/ (written only by the manager; size budget 2 GB, see BASELINE.md).
Key = sha256 of:
  * the CPU-path sources: every file under Prog/, Libraries/, Analysis/ of SRC_DIR except Prog/gpu_*.F90 and git.h
    (conservative: a change anywhere, even inside #ifdef ALF_GPU, is a miss);
  * the build: tool, configure line, ALF_FLAGS_EXT, compiler versions (from the run's build record);
  * the inputs: sha256 of parameters and seeds.
A hit is valid only if all of these are equal; nothing else is inferred. A run with a nonzero exit code, a timeout or
debug variables is never cached.
"""
import hashlib
import json
import os
import shutil
import sys

ROOT = os.path.expanduser("~/alf-gpu/refcache")
KEEP = ("parameters", "seeds", "info", "data.h5", "run.log", "env.txt", "provenance.json")


def src_hash(src):
    h = hashlib.sha256()
    for f in ("configure.sh", "Makefile"):
        h.update(f.encode() + b"\0" + open(os.path.join(src, f), "rb").read() + b"\0")
    for top in ("Prog", "Libraries", "Analysis"):
        for d, _, fs in sorted(os.walk(os.path.join(src, top))):
            for f in sorted(fs):
                rel = os.path.relpath(os.path.join(d, f), src)
                if (rel.startswith("Prog/gpu_") and rel.endswith(".F90")) or f == "git.h":
                    continue
                h.update(rel.encode() + b"\0" + open(os.path.join(d, f), "rb").read() + b"\0")
    return h.hexdigest()


def key_of(build, params_sha, seeds_sha, src):
    if build.get("source") == "base_src" or build.get("name", "").startswith("base_"):
        sys.exit("refcache is not used for base builds")
    b = {k: build.get(k) for k in ("tool", "configure", "ALF_FLAGS_EXT", "ALF_FC", "ALF_FLAGS_PROG", "ALF_LIB", "env_file_sha256")}
    b["toolchain"] = {k: build.get("toolchain", {}).get(k) for k in ("nvfortran", "gfortran", "ifx")}
    blob = json.dumps({"src": src_hash(src), "build": b, "params": params_sha, "seeds": seeds_sha}, sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest(), blob


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main(argv):
    cmd = argv[0] if argv else "list"
    if cmd in ("key", "put"):
        run, src = argv[1], argv[2]
        p = json.load(open(os.path.join(run, "provenance.json")))
        if p["exit_code"] != 0 or p.get("timed_out") or p.get("debug_env_present") or p.get("mutated"):
            sys.exit("not cacheable: failed, timed out, mutated or debug variables present")
        cand = os.path.join(os.path.dirname(os.path.abspath(src)), ".alf_candidate")
        sha = dict(l.split(None, 1) for l in open(cand).read().splitlines() if " " in l).get("sha") if os.path.isfile(cand) else None
        if sha != p["build"]["candidate"].get("sha"):
            sys.exit(f"SRC_DIR is not the source the run was built from ({sha} vs {p['build']['candidate'].get('sha')})")
        k, blob = key_of(p["build"], p["inputs"]["parameters_sha256"], p["inputs"]["seeds_sha256"], src)
        if cmd == "key":
            print(k)
            return 0
        dst = os.path.join(ROOT, k)
        if os.path.exists(dst):
            sys.exit(f"refusing: entry {k} exists")
        os.makedirs(dst)
        for f in KEEP:
            if os.path.isfile(os.path.join(run, f)):
                shutil.copy2(os.path.join(run, f), dst)
        open(os.path.join(dst, "KEY.json"), "w").write(blob)
        with open(os.path.join(dst, "MANIFEST.sha256"), "w") as m:
            for f in sorted(os.listdir(dst)):
                if f != "MANIFEST.sha256":
                    m.write(f"{sha(os.path.join(dst, f))}  {f}\n")
        print("stored", k)
        return 0
    if cmd == "get":
        new, src, bj, cfg = argv[1:5]
        build = json.load(open(bj))
        k, _ = key_of(build, sha(os.path.join(cfg, "parameters")), sha(os.path.join(cfg, "seeds")), src)
        e = os.path.join(ROOT, k)
        if not os.path.isdir(e):
            print("miss", k)
            return 1
        for line in open(os.path.join(e, "MANIFEST.sha256")):
            h, f = line.split()
            if sha(os.path.join(e, f)) != h:
                sys.exit(f"cache entry {k} corrupt: {f}")
        if os.path.exists(new):
            sys.exit(f"refusing: {new} exists")
        shutil.copytree(e, new)
        p = json.load(open(os.path.join(new, "provenance.json")))
        p["from_refcache"] = k
        json.dump(p, open(os.path.join(new, "provenance.json"), "w"), indent=1)
        print("hit", k)
        return 0
    if cmd in ("list", "check"):
        if not os.path.isdir(ROOT):
            return 0
        bad = 0
        for k in sorted(os.listdir(ROOT)):
            e = os.path.join(ROOT, k)
            size = sum(os.path.getsize(os.path.join(e, f)) for f in os.listdir(e))
            ok = True
            if cmd == "check":
                for line in open(os.path.join(e, "MANIFEST.sha256")):
                    h, f = line.split()
                    ok &= sha(os.path.join(e, f)) == h
            bad += not ok
            print(k[:16], f"{size / 2**20:8.1f} MB", "ok" if ok else "CORRUPT")
        return 1 if bad else 0
    sys.exit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
