#!/usr/bin/env python3
"""Gate types `qref` and `h5_identical` of verdict.py (ticket 207-5): judging the QREF accuracy gates of ticket 207 from a pre-registered
acceptance file.  A thin layer over the 207-2 scripts (docs/gpu/scripts/nividia/t207_d5_*): the strict parser, the aggregation and the
decision of S_ens are theirs and are called, not re-implemented; this module checks that what they are given is the harness's own
evidence (provenance of the run, the offline build record, the sha256 of the scripts) and maps their exit statuses to verdicts.

Gate `qref` (acceptance file):
  id, type "qref", path "def"|"log", stride 3 (decision 35), qs [..], keep 39, nul "0,4,..." (the null codes of the offline run),
  runs [main run ids], ctl_runs [production control run ids: TAU_NS=0 / TAUSWAP=8 or a 207-3 fault], refs "def"|"log"|[names] (the
  frozen reference ensemble; the machine's member list is a parameter: on `intel` an MKL/ifx host member replaces nvh),
  bound (S_ens(pTT) <= bound), expect_points (points with B per run), required_controls ["var8", "ctl_dev"],
  pointwise {measure "emax"|"eF", bound number|null (null only in a DRAFT file), expect_points (all points per run),
             ctl_runs (default: the gate's), required_controls},
  tools_sha256 {script: sha256} = the 207-2 scripts as pre-registered (a different script is INPUT_ERROR),
  device "dev" (a control may replace it by "var8": the judged device becomes the planted stand-in defect, which must FAIL).
Every run of the gate: provenance and active lines as every harness run (check_prov), a qproc/ directory written by acc_qproc.sh whose
copies of the records equal the harness provenance, the offline build record of the candidate (tarball hash, candidate sha, program
hash identical in every run and equal to qref_proc.json's), `qparse check --dumps` after the thinning, the kept dumps by the keep rule.
Gate `h5_identical`: pairs [[on, off]] of runs with the same build, config and seeds, `on` with ALF_GPU_TAUDIAG_QREF, `off` without any
ALF_GPU_* variable; every dataset of data.h5 byte-identical.
"""
import hashlib
import importlib.util
import json
import math
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.abspath(os.path.join(HERE, "..", "..", "..", "docs", "gpu", "scripts", "nividia"))
TOOL_FILES = ("t207_d5_proc.sh", "t207_d5_qparse.py", "t207_d5_agg.py", "t207_d5_gate.py", "t207_d5_qref.F90", "t207_d5_build.sh",
              "t207_d5_qbuild.sh")
STRIDE = 3  # decision 35
_QP = None


def qparse():
    global _QP
    if _QP is None:
        sp = importlib.util.spec_from_file_location("qref_gate_qp", os.path.join(TOOLS, "t207_d5_qparse.py"))
        _QP = importlib.util.module_from_spec(sp)
        sp.loader.exec_module(_QP)
    return _QP


def sha(p):
    with open(p, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def is_sha(x):
    return isinstance(x, str) and len(x) == 64 and all(c in "0123456789abcdef" for c in x)


def ref_arg(g):
    r = g.get("refs", g["path"])
    return r if isinstance(r, str) else ",".join(r)


# ------------------------------------------------------------------ lint
def lint(V, acc):
    errs = []
    runs, builds = acc.get("runs", {}), acc.get("builds", {})
    draft = "DRAFT" in acc.get("_status", "")
    QP = qparse()
    used = set()
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
            continue
        if g.get("type") != "qref":
            continue
        gid = g["id"]
        if g.get("path") not in ("def", "log"):
            errs.append(f"gate {gid}: path must be def or log")
            continue
        used.add(g["path"])
        if g.get("stride", STRIDE) != STRIDE:
            errs.append(f"gate {gid}: stride must be {STRIDE} (decision 35)")
        for k in ("bound", "expect_points"):
            if k not in g:
                errs.append(f"gate {gid}: '{k}' missing")
        if "bound" in g and not (isinstance(g["bound"], (int, float)) and not isinstance(g["bound"], bool) and math.isfinite(g["bound"])):
            errs.append(f"gate {gid}: bound must be a finite number")
        if not (isinstance(g.get("qs"), list) and g["qs"] and all(isinstance(x, int) and x > 0 for x in g["qs"])):
            errs.append(f"gate {gid}: qs must be a non-empty list of positive integers")
        if not (isinstance(g.get("keep"), int) and g["keep"] >= 1):
            errs.append(f"gate {gid}: keep must be an integer >= 1 (retention rule of decision 33)")
        ts = g.get("tools_sha256")
        if not (isinstance(ts, dict) and all(is_sha(ts.get(f)) for f in TOOL_FILES)):
            errs.append(f"gate {gid}: tools_sha256 must give the sha256 of {', '.join(TOOL_FILES)}")
        pw = g.get("pointwise")
        if not isinstance(pw, dict) or pw.get("measure") not in ("emax", "eF"):
            errs.append(f"gate {gid}: pointwise.measure must be emax or eF")
        else:
            pb = pw.get("bound")
            if pb is None:
                if not draft:
                    errs.append(f"gate {gid}: pointwise.bound is null (allowed only in a DRAFT file)")
            elif not (isinstance(pb, (int, float)) and not isinstance(pb, bool) and math.isfinite(pb)):
                errs.append(f"gate {gid}: pointwise.bound must be a finite number or null")
            if not isinstance(pw.get("expect_points"), int):
                errs.append(f"gate {gid}: pointwise.expect_points missing")
        mains, ctls = g.get("runs", []), g.get("ctl_runs", [])
        if not mains:
            errs.append(f"gate {gid}: no main runs")
        if not ctls:
            errs.append(f"gate {gid}: no production control run (ctl_runs): the gate needs a TAU_NS=0 / TAUSWAP=8 or graded-fault run")
        for role, ids in (("main", mains), ("ctl", ctls)):
            for rid in ids:
                r = runs.get(rid)
                if r is None:
                    continue  # reported as undeclared by the generic reference check
                used_b = builds.get(r["build"], {})
                if not used_b.get("gpu", True):
                    errs.append(f"gate {gid}: run {rid} is a CPU build (the device under test must run)")
                if not r.get("allow_debug"):
                    errs.append(f"gate {gid}: run {rid} needs allow_debug")
                q = r.get("qref")
                if not isinstance(q, dict):
                    errs.append(f"gate {gid}: run {rid} has no 'qref' entry")
                    continue
                if q.get("path") != g["path"] or q.get("keep") != g.get("keep") or q.get("qs") not in g.get("qs", []):
                    errs.append(f"gate {gid}: run {rid}: qref path/keep/qs differ from the gate's")
                if q.get("nul") != g.get("nul"):
                    errs.append(f"gate {gid}: run {rid}: qref.nul differs from the gate's nul")
                try:
                    envd = QP.parse_env_tokens([kv for kv in r.get("env", [])], f"run {rid} env")
                    got = QP.env_role(envd, g["path"])
                    if got != role:
                        errs.append(f"gate {gid}: run {rid} is declared in the {role} list, but its env is that of a {got} run")
                    if envd.get("ALF_GPU_TAUDIAG_QREF") != str(g.get("stride", STRIDE)):
                        errs.append(f"gate {gid}: run {rid}: ALF_GPU_TAUDIAG_QREF must be {g.get('stride', STRIDE)}")
                except QP.QError as e:
                    errs.append(f"gate {gid}: run {rid}: {e}")
                if any(not kv.startswith("ALF_GPU_") for kv in r.get("env", [])):
                    errs.append(f"gate {gid}: run {rid}: a QREF run may declare only ALF_GPU_* variables")
    for p in sorted(used):
        qb = acc.get("qbuilds", {}).get(p)
        if not (isinstance(qb, dict) and isinstance(qb.get("configure"), str) and isinstance(qb.get("qflags"), str)):
            errs.append(f"qbuilds.{p}: configure and qflags required (offline program of path {p})")
        elif any(t.upper().startswith("GPU=") for t in qb["configure"].split()):
            errs.append(f"qbuilds.{p}: the oracle is a CPU build")
    for rid, r in runs.items():
        if r.get("qref") and rid not in {x for g in acc.get("gates", []) if g.get("type") == "qref" for x in g.get("runs", []) + g.get("ctl_runs", [])
                                         + (g.get("pointwise") or {}).get("ctl_runs", [])} and not any(
                rid in p for g in acc.get("gates", []) if g.get("type") == "h5_identical" for p in g.get("pairs", [])):
            errs.append(f"run {rid}: has a qref entry, but no qref gate uses it")
    return errs


# ------------------------------------------------------------------ helpers of the gates
def _py(script, *args, timeout=3600):
    p = subprocess.run([sys.executable, "-I", os.path.join(TOOLS, script)] + list(args), capture_output=True, text=True, timeout=timeout)
    return p.returncode, p.stdout + p.stderr


def check_tools(V, g):
    ts = g.get("tools_sha256", {})
    for f in TOOL_FILES:
        if not is_sha(ts.get(f)):
            raise V.InputFail(f"qref gate: tools_sha256 lacks {f}")
        if sha(V.need(os.path.join(TOOLS, f))) != ts[f]:
            raise V.InputFail(f"qref gate: {f} differs from the pre-registered script (sha256)")


def check_run(V, g, acc, root, rid, role):
    """provenance, bridge records, offline build record, retention of one QREF run; returns (qref.txt path, program sha256)"""
    QP = qparse()
    rd = os.path.join(root, rid)
    p = V.check_prov(rd, acc, 0)
    decl = acc["runs"][rid]
    q = os.path.join(rd, "qproc")
    if not os.path.isdir(q):
        raise V.InputFail(f"{rid}: no qproc/ directory (the QREF processing did not run)")
    ap = json.loads(V.read(os.path.join(q, "acc_qproc.json")))
    if ap.get("schema") != "alf-acc-qproc/1" or ap.get("proc_rc") != 0:
        raise V.InputFail(f"{rid}: acc_qproc.json: processing status {ap.get('proc_rc')}")
    qd = decl.get("qref") or {}
    a = ap.get("args", {})
    if (a.get("path"), a.get("qs"), a.get("nul"), a.get("hostin"), a.get("keep")) != (g["path"], qd.get("qs"), g.get("nul"), 0, g["keep"]) \
            or qd.get("path") != g["path"] or qd.get("qs") not in g["qs"]:
        raise V.InputFail(f"{rid}: processing arguments {a} differ from the gate/run declaration")
    if any(ap.get("tools_sha256", {}).get(f) != g["tools_sha256"][f] for f in TOOL_FILES):
        raise V.InputFail(f"{rid}: the scripts used at processing time differ from the pre-registered ones")
    # the offline build record of the candidate
    qb = os.path.join(root, "..", "bin", f"qref_{g['path']}.qbuild.json")
    rec = json.loads(V.read(qb))
    if sha(qb) != ap.get("qbuild_sha256"):
        raise V.InputFail(f"{rid}: the offline build record changed after the processing, or is not the one used")
    cand = p["build"].get("candidate", {})
    if rec.get("schema") != "alf-acc-qbuild/1" or rec.get("path") != g["path"] or not is_sha(rec.get("prog_sha256")):
        raise V.InputFail(f"{rid}: bad offline build record {qb}")
    if not acc["candidate_sha"] or not rec.get("candidate", {}).get("sha", "").startswith(acc["candidate_sha"][:12]) \
            or rec["candidate"].get("sha") != cand.get("sha"):
        raise V.InputFail(f"{rid}: the offline program was not built from the candidate {acc['candidate_sha']}")
    if rec.get("tar_sha256_built") != cand.get("tar_sha256") or rec["candidate"].get("tar_sha256") != cand.get("tar_sha256"):
        raise V.InputFail(f"{rid}: the offline program was not built from the verified tarball")
    for f in TOOL_FILES:
        if rec.get("tools_sha256", {}).get(f) != g["tools_sha256"][f]:
            raise V.InputFail(f"{rid}: the offline build used a different {f} than the pre-registered one")
    if rec.get("qref_src_sha256") != rec["tools_sha256"].get("t207_d5_qref.F90"):
        raise V.InputFail(f"{rid}: offline build record: source hash inconsistent")
    # bridge: the records the 207-2 scripts read are the harness's own
    for f, want in (("info", p["outputs"]["info"]), ("run.log", p["outputs"]["run.log"]), ("parameters", p["inputs"]["parameters_sha256"])):
        if sha(V.need(os.path.join(q, f))) != want:
            raise V.InputFail(f"{rid}: qproc/{f} differs from the run's harness provenance")
    bp = json.loads(V.read(os.path.join(q, "provenance.json")))
    tok = sorted(t for t in p["extra_env"])
    if bp.get("rc") != p["exit_code"] or sorted(bp.get("env", [])) != tok or bp.get("seeds_sha256") != p["inputs"]["seeds_sha256"] \
            or bp.get("parameters_sha256") != p["inputs"]["parameters_sha256"] or bp.get("run_log_sha256") != p["outputs"]["run.log"] \
            or bp.get("info_sha256") != p["outputs"]["info"]:
        raise V.InputFail(f"{rid}: qproc/provenance.json is not the harness provenance of the run")
    # the run's own configuration: allow-list role, stride, null codes
    try:
        envd = QP.parse_env_tokens(tok, f"{rid} env")
        got = QP.env_role(envd, g["path"])
    except QP.QError as e:
        raise V.InputFail(f"{rid}: {e}")
    if got != role:
        raise V.InputFail(f"{rid}: the run is a {got} run, listed as {role}")
    if envd.get("ALF_GPU_TAUDIAG_QREF") != str(g.get("stride", STRIDE)):
        raise V.InputFail(f"{rid}: ALF_GPU_TAUDIAG_QREF {envd.get('ALF_GPU_TAUDIAG_QREF')} is not the stride {g.get('stride', STRIDE)}")
    # retention: the manifest and the kept dumps after the thinning
    rc, out = _py("t207_d5_qparse.py", "check", q, "--dumps")
    if rc != 0:
        raise V.InputFail(f"{rid}: qparse check --dumps: {out.strip()[-300:]}")
    if not os.path.isfile(os.path.join(q, "qref_thinned.txt")):
        raise V.InputFail(f"{rid}: no qref_thinned.txt (retention not done)")
    prog = json.loads(V.read(os.path.join(q, "qref_proc.json"))).get("prog")
    if prog != rec["prog_sha256"]:
        raise V.InputFail(f"{rid}: program sha256 of the processing {prog} differs from the build record {rec['prog_sha256']}")
    return os.path.join(q, "qref.txt"), prog


def swap_device(agg, a="dev", b="var8"):
    """the same aggregate with the members a and b exchanged everywhere: the control 'device := stand-in defect'"""
    def sw(d):
        if not isinstance(d, dict):
            return d
        out = {}
        for k, v in d.items():
            out[b if k == a else a if k == b else k] = v
        return out
    g = json.loads(json.dumps(agg))
    g["pooled"] = sw(g["pooled"])
    g["runs"] = {l: sw(r) for l, r in g["runs"].items()}
    for i in g["inputs"].values():
        i["members"] = [b if m == a else a if m == b else m for m in i["members"]]
    return g


def statistic(V, g, acc, root, measure, bound, expect, mains, ctls, prog, device, refs, tmp, tag, required):
    """aggregate (t207_d5_agg.py) and decide (t207_d5_gate.py) one statistic; returns (verdict, lines, numbers)"""
    args = ["%s:%s=%s" % (role, rid, path) for role, ids in (("main", mains), ("ctl", ctls)) for rid, path in ids]
    aggf = os.path.join(tmp, f"agg_{tag}.json")
    rc, out = _py("t207_d5_agg.py", aggf, measure, "--ref", refs, "--prog-sha", prog, *args)
    if rc != 0:
        raise V.InputFail(f"{tag}: aggregation: {out.strip()[-400:]}")
    agg = json.load(open(aggf))
    pooled = agg.get("pooled", {})
    miss = [c for c in required if c not in pooled]
    if miss:
        raise V.InputFail(f"{tag}: required control member(s) {miss} not in the aggregate")
    if device != "dev":
        agg = swap_device(agg, "dev", device)
    # t207_d5_gate.py is written for measure pTT (it refuses another label); the decision (device <= bound, controls > bound, the input
    # checks) does not depend on the measure, so another measure is passed under the label pTT.  To be replaced by the gate's own
    # statistic option when the 207-2 fix provides it (manager's note of 2026-11-01).
    agg["measure_judged"], agg["measure"] = measure, "pTT"
    gf = os.path.join(tmp, f"gate_{tag}.json")
    json.dump(agg, open(gf, "w"))
    cmd = ["--bound", repr(float(bound)), "--expect-points", str(expect), "--prog-sha", prog, "--require-controls"]
    for q in g["qs"]:
        cmd += ["--qs", str(q)]
    rc, out = _py("t207_d5_gate.py", *cmd, gf)
    lines = [l.strip() for l in out.splitlines() if l.strip()]
    nums = {}
    for k in ("dev", "var8", "ctl_dev"):
        if k in agg["pooled"]:
            e = agg["pooled"][k]["S_ens"]
            nums[k] = {"S_ens": e["mean"], "n": e["n"]}
    nums["bound"] = bound
    if rc == 2:
        raise V.InputFail(f"{tag}: {' '.join(lines)[-500:]}")
    if rc not in (0, 1):
        raise V.InputFail(f"{tag}: gate exit status {rc}: {' '.join(lines)[-300:]}")
    return ("PASS" if rc == 0 else "FAIL"), lines, nums


def g_qref(V, g, root, acc):
    check_tools(V, g)
    device = g.get("device", "dev")
    if device not in ("dev", "var8"):
        raise V.InputFail(f"device {device!r}")
    pw = g["pointwise"]
    ctl_pw = pw.get("ctl_runs", g.get("ctl_runs", []))
    prog, paths = None, {}
    for role, ids in (("main", g["runs"]), ("ctl", sorted(set(g.get("ctl_runs", [])) | set(ctl_pw)))):
        for rid in ids:
            qtxt, pr = check_run(V, g, acc, root, rid, role)
            if prog is not None and pr != prog:
                raise V.InputFail(f"{rid}: offline program {pr} differs from the other runs' {prog}")
            prog = pr
            paths[rid] = qtxt
    tmp = V.tmpdir()
    refs = ref_arg(g)
    req = g.get("required_controls", ["var8", "ctl_dev"])
    res, fails = {"program_sha256": prog, "runs": list(paths)}, []
    mains = [(r, paths[r]) for r in g["runs"]]
    v, lines, nums = statistic(V, g, acc, root, "pTT", g["bound"], g["expect_points"], mains, [(r, paths[r]) for r in g.get("ctl_runs", [])],
                               prog, device, refs, tmp, "pTT", req)
    res["pTT"] = {"verdict": v, "numbers": nums, "gate_output": lines}
    if v != "PASS":
        fails.append("S_ens(pTT): " + "; ".join(l for l in lines if "FAIL" in l or "NOT KILLED" in l)[:300])
    if pw.get("bound") is None:
        res["pointwise"] = {"verdict": "NOT JUDGED", "reason": "pointwise.bound is null (DRAFT file): no decision"}
    else:
        v, lines, nums = statistic(V, g, acc, root, pw["measure"], pw["bound"], pw["expect_points"], mains, [(r, paths[r]) for r in ctl_pw],
                                   prog, device, refs, tmp, "pointwise", pw.get("required_controls", req))
        res["pointwise"] = {"verdict": v, "measure": pw["measure"], "numbers": nums, "gate_output": lines}
        if v != "PASS":
            fails.append(f"pointwise S_ens({pw['measure']}): " + "; ".join(l for l in lines if "FAIL" in l or "NOT KILLED" in l)[:300])
    if fails:
        raise V.CriterionFail(" | ".join(fails))
    return res


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
