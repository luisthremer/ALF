#!/usr/bin/env python3
"""Gate types `qref` and `h5_identical` of verdict.py (ticket 207-5, second round): the QREF accuracy gates of ticket 207 judged by the
approved 207-2 gate (docs/gpu/scripts/nividia/t207_d5_gate.py --spec --spec-sha --prog --report), called on the harness's own evidence.

This module decides nothing about accuracy itself.  It checks that what the gate is given is the stamped evidence of this candidate, fills
the instance fields of the pre-registered spec from the harness records, writes the spec, runs the gate and maps its exit status:
  PASS 0 -> PASS;  FAIL 1 -> FAIL (a detection);  INPUT_ERROR 2, DEV-ONLY 3, TEST-ONLY 4, anything else -> INPUT_ERROR (never a pass, never a
  detection; DEV-ONLY = a gate no control trips, TEST-ONLY = a spec marked test_only: not acceptance verdicts).

Gate `qref` (acceptance file):
  id, type "qref", spec {the 207-2 gate spec, see t207_d5_gate.py; the instance fields binary_sha256, prog_sha, parameters_sha256 and seeds
  may be the string "@harness"}, spec_template_sha256 (sha256 of canon(spec), canon = sorted keys, no blanks, + newline), runs [main run ids],
  ctl_runs {control role of the spec: [run ids]}, tools_sha256 {207-2 script: sha256} (TOOL_FILES), _reason.
  Instance fields (known only after the candidate is built and the seeds are drawn, VERIFICATION section 3): "@harness" is replaced by the
  value of the harness records: binary_sha256 = the ALF binary of the runs (build record from the verified tarball, = provenance of every
  run), prog_sha = the offline program of qbuilds.<path> (acc_qbuild.sh from the verified tarball), parameters_sha256 = the runs' parameters
  (= configs.<c>.parameters_sha256 when pre-registered), seeds = per role the seeds_sha256 of its runs (driver's seeds record, check_prov).  A
  concrete value in the template must equal the harness value.  The filled spec's sha256 is passed as --spec-sha and recorded in the verdict.
Before the gate runs, for every run of the gate (verdict-side checks, INPUT_ERROR on any failure):
  - check_prov (harness provenance, active lines, build from the candidate tarball, seeds record, MANDATORY offload);
  - env: the recorded environment holds exactly the declared ALF_GPU_* tokens, nothing with the prefix ALF (ALFSH_*, ALF_LA_TEST_*, ...), no
    LD_PRELOAD / LD_AUDIT, OMP_TARGET_OFFLOAD=MANDATORY, thread counts 1, and otherwise only the harness's own variables (ENV_BASE) and the
    names the acceptance file lists in qref_env_allowed (what its env file sets; the lint refuses debug-like names there): everything the
    gate cannot see is refused here;
  - no symbolic link anywhere in the run directory; stamp: exactly one qref_stamp line of the attempt log for the run directory, and every
    regular file of the run directory (records, provenance, qproc/: qref.txt, qref_proc.json, qref_manifest.txt, acc_qproc.json, the kept
    dumps, ...) has the stamped sha256, no file added or missing (checked again after the gate);
  - acc_qproc.json: processing exit 0, the spec's arguments (path, qs, nul, hostin, keep), the pre-registered script hashes, the build record
    used; machine (host, /etc/machine-id, libopenblas of the program) = qref_proc.json's = the machine running this verdict;
  - offline build record bin/qref_<path>.qbuild.json: from the candidate (sha, verified tarball), source hash = spec src_sha, program on
    disk = its prog_sha256, scripts = tools_sha256;
  - bridge records qproc/{parameters,info,run.log,provenance.json,stamp.txt} = the harness provenance (acc_qproc.sh wrote them from it).
Gate `h5_identical`: pairs [[on, off]] of runs with the same build, config and seeds, `on` with ALF_GPU_TAUDIAG_QREF, `off` without any
ALF_GPU_* variable; every dataset of data.h5 byte-identical.

CLI:  qref_gate.py machine PROGRAM    the machine record (JSON) used by acc_qproc.sh and checked by the verdict
      qref_gate.py stamp RUNDIR       the stamp object (JSON) of a run directory (the driver writes it into the attempt log)
"""
import hashlib
import importlib.util
import json
import math
import os
import re
import shutil
import socket
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.abspath(os.path.join(HERE, "..", "..", "..", "docs", "gpu", "scripts", "nividia"))
TOOL_FILES = ("t207_d5_proc.sh", "t207_d5_qparse.py", "t207_d5_agg.py", "t207_d5_gate.py", "t207_d5_qref.F90", "t207_d5_build.sh",
              "t207_d5_qbuild.sh")
STRIDE = "3"  # decision 35
HARNESS = "@harness"
INSTANCE = ("binary_sha256", "prog_sha", "parameters_sha256", "seeds")
GATE_EXIT = {0: "PASS", 1: "FAIL", 2: "INPUT_ERROR", 3: "DEV-ONLY", 4: "TEST-ONLY"}
# the variables acc_run.sh itself puts into a run's environment (env -i + its own exports) and those bash adds
ENV_BASE = ("HOME", "USER", "LANG", "LC_ALL", "PATH", "TMPDIR", "PWD", "OLDPWD", "SHLVL", "_", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS", "OMP_TARGET_OFFLOAD", "ONEAPI_DEVICE_SELECTOR")
# names that may never be allowed by qref_env_allowed (lint): debug / fault / steering variables the gate cannot see
ENV_NEVER = re.compile(r"^(ALF.*|LD_PRELOAD|LD_AUDIT|OMP_.*|LIBOMPTARGET_.*|KMP_.*|GOMP_.*|NV_ACC_.*|NVCOMPILER_.*|PGI_.*|ACC_.*|CUDA_.*|CUBLAS_.*|"
                       r"CUSOLVER_.*|ZE_.*|SYCL_.*|IGC_.*|NEO.*|MKL_.*|OPENBLAS_.*|.*DEBUG.*|.*FAULT.*|ONEAPI_DEVICE_SELECTOR)$")
STAMP_SKIP = ("tmp",)  # the run's own TMPDIR is not evidence


def canon(obj):
    return (json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False) + "\n").encode()


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def is_sha(x):
    return isinstance(x, str) and re.match(r"^[0-9a-f]{64}$", x) is not None


_MODS = {}


def _load(name):
    if name not in _MODS:
        sys.dont_write_bytecode = True
        sp = importlib.util.spec_from_file_location("qref_gate_" + name, os.path.join(TOOLS, name + ".py"))
        m = importlib.util.module_from_spec(sp)
        sp.loader.exec_module(m)
        _MODS[name] = m
    return _MODS[name]


# ------------------------------------------------------------------ machine record and stamps (shared with acc_qproc.sh / acc_driver.py)
def machine_record(prog):
    """host name, sha256 of /etc/machine-id, the libopenblas the program loads (ldd) with its sha256: what the exact re-processing depends on"""
    def h(p):
        try:
            return sha(p)
        except OSError:
            return None
    blas = None
    try:
        out = subprocess.run(["ldd", prog], capture_output=True, text=True, env={"PATH": "/usr/bin:/bin"}).stdout
        for ln in out.splitlines():
            t = ln.split()
            if t and t[0].startswith("libopenblas") and len(t) >= 3 and t[1] == "=>":
                blas = os.path.realpath(t[2])
                break
    except OSError:
        pass
    return {"host": socket.gethostname(), "machine_id_sha256": h("/etc/machine-id"), "openblas": blas, "openblas_sha256": h(blas) if blas else None}


def stamp_object(rd):
    """sha256 of every regular file of a run directory (relative path; the run's tmp/ excluded); a symbolic link is an error"""
    rd = os.path.abspath(rd)
    files = {}
    for top, dirs, fs in os.walk(rd):
        rel = os.path.relpath(top, rd)
        if rel == ".":
            dirs[:] = [d for d in dirs if d not in STAMP_SKIP]
        for d in dirs:
            if os.path.islink(os.path.join(top, d)):
                raise ValueError(f"symbolic link {os.path.join(top, d)}")
        for f in fs:
            p = os.path.join(top, f)
            if os.path.islink(p) or not os.path.isfile(p):
                raise ValueError(f"symbolic link or special file {p}")
            files[os.path.normpath(os.path.join(rel, f))] = sha(p)
    return {"run_dir": rd, "files": files}


def stamp_line(who, cand, rid, obj):
    return f"qref_stamp {who} {cand} {rid} " + json.dumps(obj, sort_keys=True, separators=(",", ":"))


def find_stamp(V, acc, rd):
    lines = acc.get("_attempt_lines")
    if lines is None:
        raise V.InputFail("qref gate: no attempt log given (--attempts): the QREF evidence is bound to the run by the stamps there")
    rd = os.path.abspath(rd)
    hits = []
    for ln in lines:
        t = ln.split(None, 5)          # TIME qref_stamp WHO CANDIDATE RUNID JSON (the JSON has no blanks)
        if len(t) == 6 and t[1] == "qref_stamp":
            try:
                o = json.loads(t[5])
            except ValueError:
                continue
            if o.get("run_dir") == rd:
                hits.append((t, o))
    if len(hits) != 1:
        raise V.InputFail(f"{os.path.basename(rd)}: {len(hits)} qref_stamp lines for {rd} in the attempt log (exactly one required)")
    t, o = hits[0]
    if not acc.get("candidate_sha") or not t[3].startswith(acc["candidate_sha"][:12]):
        raise V.InputFail(f"{os.path.basename(rd)}: stamp of candidate {t[3]}, not {acc.get('candidate_sha')}")
    return o


def check_stamp(V, acc, rd, stamp=None):
    rid = os.path.basename(os.path.normpath(rd))
    stamp = stamp or find_stamp(V, acc, rd)
    try:
        now = stamp_object(rd)
    except ValueError as e:
        raise V.InputFail(f"{rid}: {e} (the harness refuses symbolic links in its run directories)")
    want = stamp.get("files", {})
    if set(now["files"]) != set(want):
        raise V.InputFail(f"{rid}: files differ from the stamp: added {sorted(set(now['files']) - set(want))[:3]}, "
                          f"missing {sorted(set(want) - set(now['files']))[:3]}")
    bad = sorted(f for f in want if now["files"][f] != want[f])
    if bad:
        raise V.InputFail(f"{rid}: {len(bad)} files changed after the stamp (first {bad[0]})")
    return stamp


# ------------------------------------------------------------------ lint
def _num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def lint(V, acc):
    errs = []
    runs, builds = acc.get("runs", {}), acc.get("builds", {})
    draft = "DRAFT" in acc.get("_status", "")
    used_paths, used_runs = set(), set()
    for n in acc.get("qref_env_allowed", []):
        if not isinstance(n, str) or ENV_NEVER.match(n) or n in ENV_BASE:
            errs.append(f"qref_env_allowed: {n!r} may not be allowed (debug, fault or steering variable, or a harness variable)")
    controls = acc.get("controls", [])
    for g in acc.get("gates", []):
        if g.get("type") == "h5_identical":
            for p in g.get("pairs", []):
                if len(p) != 2 or any(x not in runs for x in p):
                    errs.append(f"gate {g['id']}: pair {p} must name two declared runs")
                    continue
                a, b = runs[p[0]], runs[p[1]]
                if (a["build"], a["config"], a.get("seeds", "dev")) != (b["build"], b["config"], b.get("seeds", "dev")):
                    errs.append(f"gate {g['id']}: pair {p}: build, config and seeds must be the same")
                if not any(kv.startswith("ALF_GPU_TAUDIAG_QREF=") for kv in a.get("env", [])):
                    errs.append(f"gate {g['id']}: run {p[0]} has no ALF_GPU_TAUDIAG_QREF (the 'on' run)")
                if any(kv.startswith(V.DEBUG_PREFIXES) for kv in b.get("env", [])):
                    errs.append(f"gate {g['id']}: run {p[1]} (the 'off' run) carries debug variables")
                used_runs.update(p)
            continue
        if g.get("type") != "qref":
            continue
        gid = g.get("id")
        spec = g.get("spec")
        if not isinstance(spec, dict):
            errs.append(f"gate {gid}: spec (the 207-2 gate spec) missing")
            continue
        try:
            if hashlib.sha256(canon(spec)).hexdigest() != g.get("spec_template_sha256"):
                errs.append(f"gate {gid}: spec_template_sha256 is not the sha256 of the canonical spec "
                            f"({hashlib.sha256(canon(spec)).hexdigest()})")
        except (TypeError, ValueError) as e:
            errs.append(f"gate {gid}: spec not serialisable: {e}")
            continue
        if "test_only" in spec and not draft:
            errs.append(f"gate {gid}: a spec marked test_only can never PASS: only in a DRAFT file")
        ts = g.get("tools_sha256")
        if not (isinstance(ts, dict) and set(ts) == set(TOOL_FILES) and all(is_sha(v) for v in ts.values())):
            errs.append(f"gate {gid}: tools_sha256 must give exactly the sha256 of {', '.join(TOOL_FILES)}")
        path = spec.get("path")
        if path not in ("def", "log"):
            errs.append(f"gate {gid}: spec.path must be def or log")
            continue
        used_paths.add(path)
        mains, ctl = g.get("runs", []), g.get("ctl_runs", {})
        if not (isinstance(mains, list) and mains and isinstance(ctl, dict)):
            errs.append(f"gate {gid}: runs (main run ids, a non-empty list) and ctl_runs ({{role: [run ids]}}) required")
            continue
        roles = {"main": list(mains), **{r: list(v) for r, v in ctl.items()}}
        if spec.get("n_runs") != {r: len(v) for r, v in roles.items()}:
            errs.append(f"gate {gid}: spec.n_runs {spec.get('n_runs')} differs from the runs per role {({r: len(v) for r, v in roles.items()})}")
        er = spec.get("env_roles", {})
        if not isinstance(er, dict) or er.get("main", {}).get("ALF_GPU_TAUDIAG_QREF") != STRIDE:
            errs.append(f"gate {gid}: env_roles.main ALF_GPU_TAUDIAG_QREF must be {STRIDE} (decision 35)")
        for k in INSTANCE:
            if k not in spec:
                errs.append(f"gate {gid}: spec.{k} missing (a value or \"{HARNESS}\")")
        # the gate's own spec check, with the instance fields filled by dummies (unique seeds per run)
        try:
            GATE = _load("t207_d5_gate")
            dummy = json.loads(json.dumps(spec))
            for k in ("binary_sha256", "prog_sha", "parameters_sha256"):
                if dummy.get(k) == HARNESS:
                    dummy[k] = "0" * 64
            if dummy.get("seeds") == HARNESS:
                dummy["seeds"] = {r: [hashlib.sha256(f"{r}/{x}".encode()).hexdigest() for x in v] for r, v in roles.items()}
            _s, table = GATE.check_spec(dummy)
        except Exception as e:  # InputError of the gate, or anything else: the spec cannot be evaluated
            errs.append(f"gate {gid}: the 207-2 gate refuses the spec: {e}")
            table = None
        try:
            QP = _load("t207_d5_qparse")
        except Exception as e:
            errs.append(f"gate {gid}: cannot load t207_d5_qparse.py: {e}")
            QP = None
        nul = ",".join(str(c) for c in spec.get("nul", [])) if isinstance(spec.get("nul"), list) else None
        for role, ids in roles.items():
            for rid in ids:
                used_runs.add(rid)
                r = runs.get(rid)
                if r is None:
                    errs.append(f"gate {gid}: run {rid} not declared")
                    continue
                if not builds.get(r["build"], {}).get("gpu", True):
                    errs.append(f"gate {gid}: run {rid} is a CPU build (the device under test must run)")
                if not r.get("allow_debug"):
                    errs.append(f"gate {gid}: run {rid} needs allow_debug")
                q = r.get("qref")
                if not isinstance(q, dict) or set(q) != {"path", "qs", "nul", "hostin", "keep"}:
                    errs.append(f"gate {gid}: run {rid}: qref entry {{path, qs, nul, hostin, keep}} required")
                    continue
                if (q["path"], q["qs"], q["nul"], q["hostin"], q["keep"]) != (path, spec.get("qs"), nul, spec.get("hostin"), spec.get("keep")):
                    errs.append(f"gate {gid}: run {rid}: qref {q} differs from the spec's path/qs/nul/hostin/keep")
                if not (isinstance(q["keep"], int) and q["keep"] >= 1):
                    errs.append(f"gate {gid}: run {rid}: keep must be >= 1 (retention, decision 33)")
                if any(not kv.startswith("ALF_GPU_") for kv in r.get("env", [])):
                    errs.append(f"gate {gid}: run {rid}: a QREF run may declare only ALF_GPU_* variables (the gate sees no other)")
                    continue
                if table is not None and QP is not None:
                    try:
                        got = QP.role_of_env(QP.parse_env_tokens(r.get("env", []), f"run {rid} env"), table)
                        if got != role:
                            errs.append(f"gate {gid}: run {rid} is listed for role {role}, its env is that of role {got}")
                    except Exception as e:
                        errs.append(f"gate {gid}: run {rid}: {e}")
        for c in controls:
            if c.get("like") == gid:
                rp = c.get("replace", {})
                if set(rp) - {"runs"}:
                    errs.append(f"control {c.get('id')}: a qref control may replace only runs (the main runs; e.g. mutate.py qref_swap copies)")
                for rid in rp.get("runs", []):
                    used_runs.add(rid)
                    if rid not in runs or not runs[rid].get("qref"):
                        errs.append(f"control {c.get('id')}: run {rid} is not a declared QREF run")
    for p in sorted(used_paths):
        qb = acc.get("qbuilds", {}).get(p)
        if not (isinstance(qb, dict) and isinstance(qb.get("configure"), str) and isinstance(qb.get("qflags"), str)):
            errs.append(f"qbuilds.{p}: configure and qflags required (offline program of path {p})")
        elif any(t.upper().startswith("GPU=") for t in qb["configure"].split()):
            errs.append(f"qbuilds.{p}: the oracle is a CPU build")
    for rid, r in runs.items():
        if r.get("qref") and rid not in used_runs:
            errs.append(f"run {rid}: has a qref entry, but no qref gate or control uses it")
    return errs


# ------------------------------------------------------------------ the gate
def check_tools(V, g):
    ts = g.get("tools_sha256", {})
    for f in TOOL_FILES:
        if not is_sha(ts.get(f)):
            raise V.InputFail(f"qref gate: tools_sha256 lacks {f}")
        if sha(V.need(os.path.join(TOOLS, f))) != ts[f]:
            raise V.InputFail(f"qref gate: {f} differs from the pre-registered script (sha256)")


def check_env(V, rid, p, decl, acc, bdecl):
    env = p.get("env", {})
    if not env:
        raise V.InputFail(f"{rid}: no recorded environment")
    tok = sorted(f"{k}={v}" for k, v in env.items() if k.startswith("ALF_GPU_"))
    if tok != sorted(decl.get("env", [])):
        raise V.InputFail(f"{rid}: recorded ALF_GPU_* {tok} differ from the declared env {sorted(decl.get('env', []))}")
    allowed = set(ENV_BASE) | set(acc.get("qref_env_allowed", []))
    bad = sorted(k for k in env if not k.startswith("ALF_GPU_") and (k.startswith("ALF") or k in ("LD_PRELOAD", "LD_AUDIT") or k not in allowed))
    if bad:
        raise V.InputFail(f"{rid}: variables the QREF gate cannot see are refused in a QREF run: {bad} (allowed besides the declared ALF_GPU_*: "
                          f"the harness variables and qref_env_allowed)")
    if env.get("OMP_TARGET_OFFLOAD") != "MANDATORY":
        raise V.InputFail(f"{rid}: OMP_TARGET_OFFLOAD={env.get('OMP_TARGET_OFFLOAD')!r} (MANDATORY required: another value can make 'dev' a host result)")
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        if env.get(k) != "1":
            raise V.InputFail(f"{rid}: {k}={env.get(k)!r}, the harness sets 1")
    if "ONEAPI_DEVICE_SELECTOR" in env and ("GPU=intel" not in bdecl["configure"].split() or env["ONEAPI_DEVICE_SELECTOR"] != "level_zero:0"):
        raise V.InputFail(f"{rid}: ONEAPI_DEVICE_SELECTOR={env['ONEAPI_DEVICE_SELECTOR']!r} in a run that is not a GPU=intel run pinned to level_zero:0")


def check_run(V, g, acc, root, rid, role, mach_cache):
    """all verdict-side checks of one QREF run; returns its instance values and the qproc directory"""
    spec = g["spec"]
    rd = os.path.join(root, rid)
    p = V.check_prov(rd, acc, 0)
    decl = acc["runs"][rid]
    bdecl = acc["builds"][decl["build"]]
    check_env(V, rid, p, decl, acc, bdecl)
    stamp = check_stamp(V, acc, rd)
    if stamp["files"].get("provenance.json") != sha(os.path.join(rd, "provenance.json")):
        raise V.InputFail(f"{rid}: provenance.json not stamped")
    q = os.path.join(rd, "qproc")
    if not os.path.isdir(q):
        raise V.InputFail(f"{rid}: no qproc/ directory (the QREF processing did not run)")
    ap = json.loads(V.read(os.path.join(q, "acc_qproc.json")))
    if ap.get("schema") != "alf-acc-qproc/1" or ap.get("proc_rc") != 0:
        raise V.InputFail(f"{rid}: acc_qproc.json: processing status {ap.get('proc_rc')}")
    want_args = {"path": spec["path"], "qs": spec["qs"], "nul": ",".join(str(c) for c in spec["nul"]), "hostin": spec["hostin"], "keep": spec["keep"]}
    if ap.get("args") != want_args or decl.get("qref") != want_args:
        raise V.InputFail(f"{rid}: processing arguments {ap.get('args')} / declared {decl.get('qref')} differ from the spec's {want_args}")
    if any(ap.get("tools_sha256", {}).get(f) != g["tools_sha256"][f] for f in TOOL_FILES):
        raise V.InputFail(f"{rid}: the scripts used at processing time differ from the pre-registered ones")
    qbf = os.path.join(root, "..", "bin", f"qref_{spec['path']}.qbuild.json")
    rec = json.loads(V.read(qbf))
    if sha(qbf) != ap.get("qbuild_sha256"):
        raise V.InputFail(f"{rid}: the offline build record changed after the processing, or is not the one used")
    cand = p["build"].get("candidate", {})
    if rec.get("schema") != "alf-acc-qbuild/1" or rec.get("path") != spec["path"] or not is_sha(rec.get("prog_sha256")):
        raise V.InputFail(f"{rid}: bad offline build record {qbf}")
    if not cand.get("sha") or rec.get("candidate", {}).get("sha") != cand.get("sha") or not cand["sha"].startswith(acc["candidate_sha"][:12]):
        raise V.InputFail(f"{rid}: the offline program was not built from the candidate {acc['candidate_sha']}")
    if not rec.get("tar_sha256_built") or rec.get("tar_sha256_built") != cand.get("tar_sha256") or rec["candidate"].get("tar_sha256") != cand.get("tar_sha256"):
        raise V.InputFail(f"{rid}: the offline program was not built from the verified tarball of the candidate")
    if any(rec.get("tools_sha256", {}).get(f) != g["tools_sha256"][f] for f in TOOL_FILES):
        raise V.InputFail(f"{rid}: the offline build used scripts other than the pre-registered ones")
    if rec.get("qref_src_sha256") != spec.get("src_sha") or rec.get("qref_src_sha256") != g["tools_sha256"]["t207_d5_qref.F90"]:
        raise V.InputFail(f"{rid}: offline build record: program source {rec.get('qref_src_sha256')} is not the spec's src_sha")
    prog = rec.get("program", "")
    if not os.path.isfile(prog) or sha(prog) != rec["prog_sha256"]:
        raise V.InputFail(f"{rid}: offline program {prog} missing or not the one of its build record")
    if prog not in mach_cache:
        mach_cache[prog] = machine_record(prog)
    mach = mach_cache[prog]
    pr = json.loads(V.read(os.path.join(q, "qref_proc.json")))
    if ap.get("machine") != mach or pr.get("machine") != mach:
        raise V.InputFail(f"{rid}: processed on {ap.get('machine')} / {pr.get('machine')}, this verdict runs on {mach}: processing and gate must "
                          f"run on the same machine (exact re-processing depends on its libopenblas)")
    if pr.get("python") != ap.get("python") or pr.get("prog") != rec["prog_sha256"]:
        raise V.InputFail(f"{rid}: qref_proc.json interpreter/program {pr.get('python')}/{pr.get('prog')} differ from acc_qproc.json / the build record")
    # bridge: the records the 207-2 scripts read are the harness's own
    for f, want in (("info", p["outputs"]["info"]), ("run.log", p["outputs"]["run.log"]), ("parameters", p["inputs"]["parameters_sha256"])):
        if sha(V.need(os.path.join(q, f))) != want:
            raise V.InputFail(f"{rid}: qproc/{f} differs from the run's harness provenance")
    bp = json.loads(V.read(os.path.join(q, "provenance.json")))
    tok = sorted(p["extra_env"])
    if bp.get("rc") != p["exit_code"] or sorted(bp.get("env", [])) != tok or bp.get("seeds_sha256") != p["inputs"]["seeds_sha256"] \
            or bp.get("parameters_sha256") != p["inputs"]["parameters_sha256"] or bp.get("run_log_sha256") != p["outputs"]["run.log"] \
            or bp.get("info_sha256") != p["outputs"]["info"] or bp.get("binary_sha256") != p["binary_sha256"]:
        raise V.InputFail(f"{rid}: qproc/provenance.json is not the harness provenance of the run")
    st = V.read(os.path.join(q, "stamp.txt")).splitlines()
    want_st = [f"binary {p['binary_sha256']}", "env " + " ".join(t for t in p["extra_env"]),
               f"seeds {p['inputs']['seeds_sha256'][:16]} parameters {p['inputs']['parameters_sha256'][:16]}"]
    if st != want_st:
        raise V.InputFail(f"{rid}: qproc/stamp.txt is not generated from the harness provenance")
    cfg = acc.get("configs", {}).get(decl["config"], {})
    if cfg.get("parameters_sha256") and cfg["parameters_sha256"] != p["inputs"]["parameters_sha256"]:
        raise V.InputFail(f"{rid}: parameters differ from the pre-registered parameters_sha256 of config {decl['config']}")
    return {"rid": rid, "role": role, "qdir": q, "stamp": stamp, "binary_sha256": p["binary_sha256"], "prog_sha": rec["prog_sha256"], "program": prog,
            "parameters_sha256": p["inputs"]["parameters_sha256"], "seeds": p["inputs"]["seeds_sha256"]}


def fill_spec(V, spec, infos):
    s = json.loads(json.dumps(spec))
    vals = {}
    for k in ("binary_sha256", "prog_sha", "parameters_sha256"):
        got = sorted({x[k] for x in infos})
        if len(got) != 1:
            raise V.InputFail(f"the runs of the gate have different {k}: {[x[:16] for x in got]}")
        vals[k] = got[0]
    roles = {}
    for x in infos:
        roles.setdefault(x["role"], []).append(x["seeds"])
    vals["seeds"] = {r: v for r, v in sorted(roles.items())}
    filled = {}
    for k in INSTANCE:
        if s.get(k) == HARNESS:
            s[k] = vals[k]
            filled[k] = vals[k]
        elif k == "seeds":
            if not isinstance(s[k], dict) or {r: sorted(v) for r, v in s[k].items()} != {r: sorted(v) for r, v in vals[k].items()}:
                raise V.InputFail(f"the pre-registered seeds {s[k]} differ from the seed sets of the runs {vals[k]}")
        elif s.get(k) != vals[k]:
            raise V.InputFail(f"the pre-registered {k} {s.get(k)} differs from the harness value {vals[k]}")
    return s, filled


def g_qref(V, g, root, acc):
    spec = g["spec"]
    if hashlib.sha256(canon(spec)).hexdigest() != g.get("spec_template_sha256"):
        raise V.InputFail("qref gate: the spec differs from spec_template_sha256 of the acceptance file")
    check_tools(V, g)
    roles = [("main", rid) for rid in g["runs"]] + [(r, rid) for r, ids in sorted(g["ctl_runs"].items()) for rid in ids]
    if len({rid for _, rid in roles}) != len(roles):
        raise V.InputFail("qref gate: a run listed twice")
    mach_cache = {}
    infos = [check_run(V, g, acc, root, rid, role, mach_cache) for role, rid in roles]
    filled, inst = fill_spec(V, spec, infos)
    raw = canon(filled)
    ssha = hashlib.sha256(raw).hexdigest()
    out = acc.get("_out") or V.tmpdir()
    tag = re.sub(r"[^A-Za-z0-9_.-]", "_", g["id"])
    sp, rp = os.path.join(out, f"qref_{tag}.spec.json"), os.path.join(out, f"qref_{tag}.report.json")
    for f in (sp, rp):
        if os.path.lexists(f):
            raise V.InputFail(f"refusing to overwrite {f}")
    with open(sp, "wb") as f:
        f.write(raw)
    env = {k: os.environ[k] for k in ("PATH", "HOME") if k in os.environ}
    env["OMP_NUM_THREADS"] = "4"
    cmd = [sys.executable, "-I", os.path.join(TOOLS, "t207_d5_gate.py"), "--spec", sp, "--spec-sha", ssha, "--prog", infos[0]["program"], "--report", rp] \
        + [x["qdir"] for x in infos]
    pr = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=out)
    lines = [l for l in pr.stdout.splitlines() if l.strip()]
    for x in infos:  # nothing changed while the gate ran
        check_stamp(V, acc, os.path.join(root, x["rid"]), x["stamp"])
    verdict = GATE_EXIT.get(pr.returncode, f"exit {pr.returncode}")
    try:
        rep = json.load(open(rp))
    except (OSError, ValueError) as e:
        raise V.InputFail(f"qref gate: no readable report ({type(e).__name__}: {e}); gate {verdict}: {' '.join(lines[-2:])[-300:]}")
    if rep.get("exit_code") != pr.returncode or rep.get("verdict") != verdict or not lines or not lines[-1].startswith("GATE " + verdict):
        raise V.InputFail(f"qref gate: report ({rep.get('verdict')}, {rep.get('exit_code')}) and exit status {pr.returncode} / last line disagree")
    if pr.returncode in (0, 1) and rep.get("spec_sha256") != ssha:
        raise V.InputFail("qref gate: the report's spec sha256 is not the one passed")
    res = {"gate_verdict": verdict, "spec_sha256": ssha, "spec_file": sp, "report": rp, "report_sha256": sha(rp), "filled": inst,
           "runs": {x["rid"]: x["role"] for x in infos}, "program": infos[0]["program"],
           "gates": [{k: x.get(k) for k in ("name", "measure", "statistic", "pool", "bound", "device_value", "passed")} for x in rep.get("gates", [])],
           "controls": [{"name": c["name"], "checks": [{k: ch.get(k) for k in ("gate", "must_trip", "value", "ok")} for ch in c.get("checks", [])]}
                        for c in rep.get("controls", [])],
           "authenticity": rep.get("authenticity"), "reasons": rep.get("reasons"), "last_line": lines[-1]}
    if pr.returncode == 0:
        return res
    summary = "; ".join(rep.get("reasons", []))[:600]
    if pr.returncode == 1:
        raise V.CriterionFail(f"gate FAIL: {summary} | " + " ".join(f"{x['name']}: device {x['device_value']:+.4f} bound {x['bound']:+.4f}"
                                                                     for x in res["gates"] if x.get("device_value") is not None)
                              + f" | report {rp}")
    raise V.InputFail(f"gate {verdict}: {summary} | report {rp}")


def g_h5_identical(V, g, root, acc):
    """data.h5 with QREF on and off (same build, config, seeds): every dataset byte-identical"""
    res = []
    for on, off in g["pairs"]:
        pa, pb = V.check_prov(os.path.join(root, on), acc, 0), V.check_prov(os.path.join(root, off), acc, 0)
        if not any(t.startswith("ALF_GPU_TAUDIAG_QREF=") for t in pa["extra_env"]) or any(t.startswith("ALF_GPU_") for t in pb["extra_env"]):
            raise V.InputFail(f"{on}/{off}: the 'on' run must have ALF_GPU_TAUDIAG_QREF and the 'off' run no ALF_GPU_* variable")
        if pa["inputs"] != pb["inputs"] or pa["build"].get("name") != pb["build"].get("name"):
            raise V.InputFail(f"{on}/{off}: build, parameters or seeds differ")
        da, db = V.h5_datasets(os.path.join(root, on, "data.h5")), V.h5_datasets(os.path.join(root, off, "data.h5"))
        if set(da) != set(db):
            raise V.InputFail(f"dataset sets differ: {sorted(set(da) ^ set(db))[:5]}")
        if len(da) < g.get("min_datasets", 1):
            raise V.InputFail(f"{len(da)} datasets, expected at least {g.get('min_datasets', 1)}")
        bad = [n for n in sorted(da) if da[n].shape != db[n].shape or da[n].dtype != db[n].dtype or da[n].tobytes() != db[n].tobytes()]
        res.append({"pair": [on, off], "datasets": len(da), "differing": len(bad)})
        if bad:
            raise V.CriterionFail(f"{on} vs {off}: {len(bad)} of {len(da)} datasets differ (first: {bad[0]})")
    return res


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "machine":
        print(json.dumps(machine_record(sys.argv[2]), sort_keys=True))
    elif len(sys.argv) == 3 and sys.argv[1] == "stamp":
        print(json.dumps(stamp_object(sys.argv[2]), sort_keys=True, separators=(",", ":")))
    else:
        sys.exit("usage: qref_gate.py machine PROGRAM | stamp RUNDIR")
