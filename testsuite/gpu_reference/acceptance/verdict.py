#!/usr/bin/env python3
"""Acceptance verdicts from raw run outputs (ALF GPU port, standard of 2026-10-29).

usage: verdict.py ACCEPTANCE.json RUNROOT [--out DIR] [--candidate SHA] [--seeds-record JSON]
       verdict.py --selfcheck            positive controls of this script (synthetic data, CPU only, seconds)

RUNROOT holds one directory per run id of the acceptance file (written by acc_run.sh: provenance.json, info, run.log,
data.h5, ...). Writes verdict.json and verdict.md (refuses to overwrite them). Needs TMPDIR (confinement).

Outcomes per gate: PASS; FAIL (the gate's criterion is not met); INPUT_ERROR (the gate could not be evaluated:
missing run/file/line/dataset, provenance mismatch, crash, timeout, NaN/Inf, count mismatch, GPU path inactive,
undeclared or leaked debug variable, GPU busy at start, mutated copy outside a control). A gate passes only with PASS.
A control is "killed" only by its expected outcome: normally FAIL, i.e. the criterion itself caught the planted defect;
an INPUT_ERROR never counts as a detection unless the control expects INPUT_ERROR (e.g. a NaN or a dropped dataset).
Exit code: 0 only if every gate PASSes and every control gives its expected outcome and nothing is BLOCKED; 3 if
all of that holds but the acceptance file lists BLOCKED items (never a pass); 1 otherwise; 2 on a usage error.
"""
import argparse
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile

FLOAT = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?"


class InputFail(Exception):
    """The gate could not be evaluated (never a detection)."""


class CriterionFail(Exception):
    """The gate's criterion is not met (a detection)."""


def num(s):
    v = float(s.replace("D", "E").replace("d", "e"))
    if not math.isfinite(v):
        raise InputFail(f"non-finite number {s}")
    return v


def need(path):
    if not os.path.exists(path):
        raise InputFail(f"missing {path}")
    return path


def read(path):
    with open(need(path), errors="replace") as f:
        return f.read()


def tmpdir():
    base = os.environ.get("TMPDIR")
    if not base or not os.path.isdir(base):
        raise InputFail("TMPDIR is not set to an existing directory (confinement rule)")
    return tempfile.mkdtemp(dir=base)


# ------------------------------------------------------------------ provenance
DEBUG_PREFIXES = ("ALF_GPU_", "ALFSH_", "ALF_LA_TEST_")


def check_prov(rd, acc, expect_rc=0, gpu_exclusive=None):
    """Common checks; the run id is the directory name and must be declared in acc['runs']."""
    p = json.loads(read(os.path.join(rd, "provenance.json")))
    rid = os.path.basename(os.path.normpath(rd))
    if p.get("schema") != "alf-acc-prov/1":
        raise InputFail(f"{rid}: unknown provenance schema")
    rc = p.get("exit_code")
    if p.get("timed_out") or rc in (124, 137):
        raise InputFail(f"{rid}: timed out (rc {rc})")
    if expect_rc is not None and rc != expect_rc:
        raise InputFail(f"{rid}: exit code {rc}, expected {expect_rc}")
    decl = acc.get("runs", {}).get(rid)
    if decl is None:
        raise InputFail(f"{rid}: run not declared in the acceptance file")
    b = p.get("build", {})
    bdecl = acc.get("builds", {}).get(decl["build"])
    if bdecl is None or b.get("name") != decl["build"]:
        raise InputFail(f"{rid}: built as {b.get('name')}, declared {decl['build']}")
    if b.get("configure", "").split() != bdecl["configure"].split() or b.get("ALF_FLAGS_EXT", "") != bdecl.get("flags_ext", ""):
        raise InputFail(f"{rid}: build configuration differs from the acceptance file")
    if p.get("from_refcache"):
        if bdecl.get("gpu", True) or decl["build"].startswith("base_"):
            raise InputFail(f"{rid}: refcache hit for a GPU or base build")
    else:
        want = acc.get("base_sha") if b.get("source") == "base_src" else acc.get("candidate_sha")
        cand = b.get("candidate", {})
        have = cand.get("base_sha") if b.get("source") == "base_src" else cand.get("sha")
        if not want or len(want) < 12 or not have or not have.startswith(want):
            raise InputFail(f"{rid}: built from {have}, acceptance names {want} (>= 12 hex digits needed)")
    if p.get("binary_sha256") != b.get("binary_sha256"):
        raise InputFail(f"{rid}: binary sha256 differs from its build record")
    tkey = "base_tar_sha256" if b.get("source") == "base_src" else "tar_sha256"
    if not b.get("tar_sha256_built") or b.get("tar_sha256_built") != b.get("candidate", {}).get(tkey):
        raise InputFail(f"{rid}: binary not built from the verified tarball ({b.get('tar_sha256_built')} vs {b.get('candidate', {}).get(tkey)})")
    for f, h in (p.get("outputs") or {}).items():
        fp = os.path.join(rd, f)
        if h is None:
            if os.path.exists(fp):
                raise InputFail(f"{rid}: {f} appeared after the run (absent at run time)")
        elif hashlib.sha256(open(need(fp), "rb").read()).hexdigest() != h:
            raise InputFail(f"{rid}: {f} changed after the run (sha256 differs from provenance)")
    if not p.get("outputs") or p["outputs"].get("run.log") is None:
        raise InputFail(f"{rid}: no output hashes in provenance")
    if sorted(p.get("extra_env", [])) != sorted(decl.get("env", [])):
        raise InputFail(f"{rid}: environment {p.get('extra_env')} differs from the declared {decl.get('env', [])}")
    declared = {kv.split("=", 1)[0] for kv in decl.get("env", [])}
    # from the recorded environment as well, not only acc_run.sh's list (older provenance may lack a prefix)
    debug = sorted(k for k in set(p.get("debug_env_present", [])) | set(p.get("env", {}))
                   | {kv.split("=", 1)[0] for kv in p.get("extra_env", [])} if k.startswith(DEBUG_PREFIXES))
    if [k for k in debug if k not in declared]:
        raise InputFail(f"{rid}: undeclared debug variables {sorted(set(debug) - declared)}")
    if debug and not decl.get("allow_debug"):
        raise InputFail(f"{rid}: debug variables in a run not declared allow_debug")
    rec = acc.get("_seeds_record")
    if rec is not None:
        key = f"{decl.get('seeds', 'dev')}/{decl['config']}"
        if rec.get(key) != p.get("inputs", {}).get("seeds_sha256"):
            raise InputFail(f"{rid}: seeds differ from the driver's record for {key}")
    excl = acc.get("gpu_exclusive", True) if gpu_exclusive is None else gpu_exclusive
    if bdecl.get("gpu", True) and excl:
        for k in ("gpu_check_rc", "gpu_before", "load_ok"):
            if k not in p:
                raise InputFail(f"{rid}: provenance lacks {k}")
        if acc.get("machine") == "intel":  # no nvidia-smi: the user-load check stands in for the idle check
            if p["gpu_check_rc"] != -1 or not p["load_ok"]:
                raise InputFail(f"{rid}: intel load check failed or unexpected GPU tool (rc {p['gpu_check_rc']})")
        else:
            if p["gpu_check_rc"] != 0:
                raise InputFail(f"{rid}: GPU idle check failed (rc {p['gpu_check_rc']})")
            if [l for l in p["gpu_before"] if l.strip()]:
                raise InputFail(f"{rid}: GPU not idle at start: {p['gpu_before']}")
    if bdecl.get("gpu", True) and not p.get("from_refcache"):
        if p.get("env", {}).get("OMP_TARGET_OFFLOAD") != "MANDATORY":
            raise InputFail(f"{rid}: GPU run without OMP_TARGET_OFFLOAD=MANDATORY (host fallback possible)")
        if "GPU=intel" in bdecl["configure"].split():  # acc_run.sh pins the discrete GPU (the iGPU has no FP64)
            env = p.get("env", {})
            if env.get("ONEAPI_DEVICE_SELECTOR") != "level_zero:0" or "ZE_AFFINITY_MASK" in env:
                raise InputFail(f"{rid}: GPU=intel run without ONEAPI_DEVICE_SELECTOR=level_zero:0 or with ZE_AFFINITY_MASK")
        text = read(os.path.join(rd, "run.log"))
        if os.path.isfile(os.path.join(rd, "info")):
            text += read(os.path.join(rd, "info"))
        need_l = decl.get("active_lines", bdecl.get("active_lines", acc.get("active_lines", [])))
        if decl.get("active_lines") is not None and not decl.get("allow_debug"):
            raise InputFail(f"{rid}: only declared debug runs may override the active lines")
        missing = [l for l in need_l if l not in text]
        if missing or re.search(r"\bINACTIVE\b", text):
            raise InputFail(f"{rid}: GPU path not active (missing {missing[:1]}, or an INACTIVE line)")
    if p.get("mutated") and not acc.get("_allow_mutated"):
        raise InputFail(f"{rid}: mutated copy ({p['mutated']}) used outside a control")
    return p


# ------------------------------------------------------------------ parsers
def info_precision(rd):
    t = read(os.path.join(rd, "info")).replace("\n", " ")
    out = {}
    for key, label, n in (("green", "Precision Green", 2), ("tau", "Precision tau", 2), ("phase", "Precision Phase", 1)):
        m = re.search(label + r"[^:]*:\s+" + r"\s+".join([f"({FLOAT})"] * n), t)
        if not m:
            if key == "tau" and "Precision tau" not in t:
                continue
            raise InputFail(f"{rd}: cannot parse '{label}' in info")
        out[key] = [num(x) for x in m.groups()]
    return out


def h5_datasets(path):
    import h5py
    import numpy as np
    out = {}
    with h5py.File(need(path), "r") as f:
        f.visititems(lambda n, o: out.__setitem__(n, np.asarray(o[()])) if isinstance(o, h5py.Dataset) else None)
    if not out:
        raise InputFail(f"{path}: no datasets")
    return out


def taudiag(rd):
    rows = []
    for line in read(os.path.join(rd, "run.log")).splitlines():
        if not line.lstrip().startswith("TAUDIAG pt "):  # data lines only, not "GPU debug env: ...ALF_GPU_TAUDIAG=..." (212)
            continue
        t = line.split()
        try:
            i, j, k = t.index("dev"), t.index("hostfull"), t.index("P")
            r = {"pass": int(t[t.index("pass") + 1]), "nst": int(t[t.index("nst") + 1]),
                 "dev": [num(x) for x in t[i + 1:i + 5]], "host": [num(x) for x in t[j + 1:j + 5]],
                 "P": [num(x) for x in t[k + 1:k + 6]]}
        except (ValueError, IndexError) as e:
            raise InputFail(f"{rd}: malformed TAUDIAG line: {line[:80]} ({e})")
        if len(r["dev"]) != 4 or len(r["host"]) != 4 or len(r["P"]) != 5:
            raise InputFail(f"{rd}: short TAUDIAG line: {line[:80]}")
        rows.append(r)
    if not rows:
        raise InputFail(f"{rd}: no TAUDIAG lines")
    return rows


def numbers_after(line, after):
    """All numbers on `line` after the regex `after`, up to an optional '(tol' marker; parenthesised slice numbers
    such as '(worst slice 11)' or '(40)' are dropped."""
    m = re.search(after, line)
    if not m:
        return None
    rest = line[m.end():].split("(tol")[0]
    rest = re.sub(r"\([^)]*\)", " ", rest)
    rest = re.sub(r"\b[A-Z][A-Z0-9]*\d\b", " ", rest)  # labels such as G00, GT0, G0T carry digits
    rest = re.sub(r"\b[A-Za-z_]+=\d+\b", " ", rest)  # labels such as NVAR=2
    return [num(x) for x in re.findall(FLOAT, rest)]


# ------------------------------------------------------------------ gates
def g_selftest(g, root, acc):
    """Self-test runs: rc 0, the required PASSED line, no FAILED, every checked number finite and below its
    tolerance, integer conditions (e.g. both CGR2_2 branches hit) met."""
    res = []
    for r in g["runs"]:
        rd = os.path.join(root, r)
        p = check_prov(rd, acc, expect_rc=None)
        log = read(os.path.join(rd, "run.log"))
        fpat = g.get("fail_regex", r"GPU self-test FAILED: ")
        fl = [l.strip() for l in log.splitlines() if re.search(fpat, l)]
        if p["exit_code"] == 1 and fl:
            raise CriterionFail(f"{r}: self-test FAILED: {fl[:1]}")
        if p["exit_code"] != 0:
            raise InputFail(f"{r}: exit code {p['exit_code']} without a self-test FAILED line (crash or setup error)")
        if fl:
            raise CriterionFail(f"{r}: self-test FAILED line with rc 0")
        for pat in g.get("require_lines", []):
            if not re.search(pat, log):
                raise InputFail(f"{r}: required line missing: {pat}")
        if re.search(r"\b(?:NaN|Infinity)\b", log, re.I):
            raise InputFail(f"{r}: NaN/Inf in the log")
        out = {"run": r}
        res.append(out)
        for chk in g.get("number_checks", []):
            vals = []
            for line in log.splitlines():
                if re.search(chk["line"], line):
                    v = numbers_after(line, chk["after"])
                    if v is None:
                        raise InputFail(f"{r}: '{chk['after']}' not found in a line matching '{chk['line']}'")
                    vals += v
            if len(vals) < chk.get("min_count", 1):
                raise InputFail(f"{r}: {len(vals)} numbers for '{chk['line']}', expected >= {chk.get('min_count', 1)}")
            out[chk["name"]] = max(vals)
            if max(vals) > chk["max"]:
                raise CriterionFail(f"{r}: {chk['name']} {max(vals):.3e} > {chk['max']:.1e}")
        for chk in g.get("int_checks", []):
            m = re.search(chk["regex"], log)
            if not m:
                raise InputFail(f"{r}: no match for {chk['regex']}")
            vals = [int(x) for x in m.groups()]
            out[chk["name"]] = vals
            if min(vals) < chk["min"]:
                raise CriterionFail(f"{r}: {chk['name']} {vals} below {chk['min']}")
    return res


def g_faults(g, root, acc):
    """Fault controls: rc must be exactly 1 with a FAILED line naming the expected sub-test. rc 0 = fault not
    detected (criterion FAIL); any other rc (crash, timeout) = INPUT_ERROR, never a detection."""
    res = []
    for r, expect in g["cases"].items():
        rd = os.path.join(root, r)
        p = check_prov(rd, acc, expect_rc=None)
        log = read(os.path.join(rd, "run.log"))
        for pat in g.get("armed_lines", []):
            if not re.search(pat, log):
                raise InputFail(f"{r}: fault not shown as armed ({pat})")
        fl = [l for l in log.splitlines() if re.search(g.get("fail_regex", r"GPU self-test FAILED: "), l)]
        if p["exit_code"] == 0:
            raise CriterionFail(f"{r}: fault not detected (rc 0)")
        if p["exit_code"] != 1 or not fl:
            raise InputFail(f"{r}: exit code {p['exit_code']} without a self-test FAILED line (a crash or timeout is not a detection)")
        if not any(expect in l for l in fl):
            raise CriterionFail(f"{r}: FAILED line does not name '{expect}': {fl[:1]}")
        res.append({"run": r, "failed_line": fl[0].strip()[:120]})
    return res


def _dev(a, b, nbin):
    import numpy as np
    if a.shape != b.shape:
        raise InputFail(f"shape {a.shape} vs {b.shape}")
    if a.ndim == 0 or a.shape[0] != nbin:
        raise InputFail(f"{a.shape[0] if a.ndim else 0} bins, expected {nbin}")
    if not (np.isfinite(a).all() and np.isfinite(b).all()):
        raise InputFail("non-finite values")
    scale = max(float(np.abs(a).max()), 1e-300)
    d = np.abs(a - b).reshape(nbin, -1).max(axis=1) / scale
    per = (np.abs(a - b) / np.maximum(np.abs(a), 1e-6 * scale)).reshape(nbin, -1).max(axis=1)
    return d, per


def g_same_seed(g, root, acc):
    """GPU vs CPU oracle, same seeds: per-bin deviation of every dataset (obser, sign, back) of data.h5.
    pairs: [run_a, run_b, nbin]."""
    import numpy as np
    res = []
    for ra, rb, nbin in g["pairs"]:
        A, B = os.path.join(root, ra), os.path.join(root, rb)
        pa_ = check_prov(A, acc, 0)
        pb_ = check_prov(B, acc, 0, gpu_exclusive=False)
        if pa_["inputs"]["parameters_sha256"] != pb_["inputs"]["parameters_sha256"]:
            raise InputFail(f"{ra}/{rb}: parameters differ")
        if read(os.path.join(A, "seeds")) != read(os.path.join(B, "seeds")):
            raise InputFail(f"{ra}/{rb}: seeds differ (same-seed gate)")
        da, db = h5_datasets(os.path.join(A, "data.h5")), h5_datasets(os.path.join(B, "data.h5"))
        if set(da) != set(db):
            raise InputFail(f"dataset sets differ: {sorted(set(da) ^ set(db))}")
        if len(da) != g["n_datasets"]:
            raise InputFail(f"{len(da)} datasets, expected {g['n_datasets']}")
        tau = [n for n in da if "_tau/" in n and n.endswith("obser")]
        if len(tau) != g["n_tau_datasets"]:
            raise InputFail(f"{len(tau)} tau obser datasets, expected {g['n_tau_datasets']}")
        dev_all = dev_tau = per_all = None
        for n in sorted(da):
            d, per = _dev(da[n], db[n], nbin)
            dev_all = d if dev_all is None else np.maximum(dev_all, d)
            per_all = per if per_all is None else np.maximum(per_all, per)
            if n in tau:
                dev_tau = d if dev_tau is None else np.maximum(dev_tau, d)
        out = {"pair": [ra, rb], "overall_max": float(dev_all.max()),
               "tau_max": float(dev_tau.max()) if dev_tau is not None else None, "per_entry_max": float(per_all.max()),
               "first_bin_above_1e-9": int(np.argmax(dev_all > 1e-9)) + 1 if (dev_all > 1e-9).any() else None}
        res.append(out)
        if (dev_all > g["no_bin_above"]).any():
            raise CriterionFail(f"{ra}: a bin above {g['no_bin_above']:.0e} (divergent trajectory)")
        if out["overall_max"] > g["overall_max"]:
            raise CriterionFail(f"{ra}: overall max {out['overall_max']:.3e} > {g['overall_max']:.1e}")
        if out["tau_max"] is not None and out["tau_max"] > g["tau_max"]:
            raise CriterionFail(f"{ra}: tau max {out['tau_max']:.3e} > {g['tau_max']:.1e}")
        if "per_entry_max" in g and out["per_entry_max"] > g["per_entry_max"]:
            raise CriterionFail(f"{ra}: per-entry max {out['per_entry_max']:.3e} > {g['per_entry_max']:.1e}")
    return res


def precision_lines(text):
    """The Precision lines of info, each with its numeric continuation line (the Max of "Mean, Max" entries)."""
    L, out = text.splitlines(), []
    for i, l in enumerate(L):
        if "Precision" in l:
            out.append(l.strip())
            if i + 1 < len(L) and re.fullmatch(rf"\s*(?:{FLOAT}\s*)+", L[i + 1]):  # NaN/Inf: the precision gate's job
                out.append(L[i + 1].strip())
    return out


def g_replay(g, root, acc):
    """Bit identity of two runs (re-run of one binary, two builds of one commit, base vs candidate on GPU=none).
    "datasets": "equal_time" compares only non-tau datasets (non-interference gate)."""
    import numpy as np
    res = []
    for ra, rb in g["pairs"]:
        A, B = os.path.join(root, ra), os.path.join(root, rb)
        check_prov(A, acc, 0, gpu_exclusive=g.get("gpu_exclusive"))
        check_prov(B, acc, 0, gpu_exclusive=g.get("gpu_exclusive"))
        da, db = h5_datasets(os.path.join(A, "data.h5")), h5_datasets(os.path.join(B, "data.h5"))
        if g.get("datasets") == "equal_time":
            da = {k: v for k, v in da.items() if "_tau/" not in k}
            db = {k: v for k, v in db.items() if "_tau/" not in k}
        if set(da) != set(db):
            raise InputFail(f"dataset sets differ: {sorted(set(da) ^ set(db))}")
        if len(da) != g["n_datasets"]:
            raise InputFail(f"{len(da)} datasets, expected {g['n_datasets']}")
        for n in da:
            if np.issubdtype(da[n].dtype, np.number) and not (np.isfinite(da[n]).all() and np.isfinite(db[n]).all()):
                raise InputFail(f"{ra}: non-finite values in {n}")
        diff = [n for n in sorted(da) if da[n].shape != db[n].shape or da[n].tobytes() != db[n].tobytes()]
        if g.get("datasets") != "equal_time":
            pa, pb = precision_lines(read(os.path.join(A, "info"))), precision_lines(read(os.path.join(B, "info")))
            if not pa:
                raise InputFail(f"{ra}: no Precision lines in info")
            if pa != pb:
                diff.append("info:Precision")
        res.append({"pair": [ra, rb], "n": len(da), "differing": len(diff)})
        if diff:
            raise CriterionFail(f"{ra} vs {rb}: {len(diff)} differ, first {diff[0]}")
    return res


def g_paired(g, root, acc):
    """Paired same-trajectory test (ticket 207): GTT pair (index 1), passes >= from_pass; reference row
    host+1ulp(b)/host+1ulp(a) = P[2]/P[1]; counts over all points."""
    res = []
    for r in g["runs"]:
        rd = os.path.join(root, r)
        check_prov(rd, acc, 0)
        rows = taudiag(rd)
        if len(rows) != g["n_points"]:
            raise InputFail(f"{r}: {len(rows)} TAUDIAG lines, expected {g['n_points']} (incomplete run?)")
        sel = [x for x in rows if x["pass"] >= g["from_pass"]]
        if len(sel) != g["n_used"]:
            raise InputFail(f"{r}: {len(sel)} rows used, expected {g['n_used']}")
        for x in sel:
            if not all(v > 0 for v in (x["dev"][1], x["host"][1], x["P"][1], x["P"][2])):
                raise InputFail(f"{r}: non-positive value in a used row (pass {x['pass']})")
        lr = [math.log10(x["dev"][1] / x["host"][1]) for x in sel]
        ref = [math.log10(x["P"][2] / x["P"][1]) for x in sel]
        sd, sref, mean = statistics.pstdev(lr), statistics.pstdev(ref), statistics.mean(lr)
        if sref <= 0:
            raise InputFail(f"{r}: zero reference spread")
        cnt = {}
        for t in g["count_thresholds"]:
            nd, nh = sum(x["dev"][1] > t for x in rows), sum(x["host"][1] > t for x in rows)
            if nd + nh < g.get("count_min_events", 1):
                raise InputFail(f"{r}: only {nd + nh} events above {t}: the count condition has no power")
            cnt[str(t)] = [nd, nh, abs(nd - nh) <= g["count_k"] * math.sqrt(nd + nh)]
        out = {"run": r, "n": len(sel), "mean": mean, "sd": sd, "ref_sd": sref, "sd_ratio": sd / sref, "counts": cnt}
        res.append(out)
        bad = []
        if abs(mean) > g["abs_mean_max"]:
            bad.append(f"|mean| {abs(mean):.4f} > {g['abs_mean_max']}")
        if sd / sref > g["sd_ratio_max"]:
            bad.append(f"sd ratio {sd / sref:.3f} > {g['sd_ratio_max']}")
        bad += [f"counts above {t}: {v[0]} vs {v[1]}" for t, v in cnt.items() if not v[2]]
        if bad:
            raise CriterionFail(f"{r}: " + "; ".join(bad))
    return res


def g_precision(g, root, acc):
    """Precision mean ratio GPU / larger CPU mean <= limit (Green and tau); the max is reported only."""
    res = []
    for case in g["cases"]:
        gp_dir = os.path.join(root, case["gpu"])
        pg = check_prov(gp_dir, acc, 0)["inputs"]
        gp = info_precision(gp_dir)
        cps = []
        for c in case["cpu"]:
            pc = check_prov(os.path.join(root, c), acc, 0, gpu_exclusive=False)["inputs"]
            if pc != pg:
                raise InputFail(f"{case['gpu']} vs {c}: parameters or seeds differ")
            cps.append(info_precision(os.path.join(root, c)))
        out = {"gpu": case["gpu"]}
        res.append(out)
        for key in g.get("quantities", ["green", "tau"]):
            if key not in gp or any(key not in c for c in cps):
                raise InputFail(f"{case['gpu']}: Precision {key} missing")
            cm = max(c[key][0] for c in cps)
            if gp[key][0] <= 0:
                raise InputFail(f"{case['gpu']}: GPU mean {key} is not positive")
            if cm <= 0:
                raise InputFail(f"{case['gpu']}: CPU mean {key} is zero")
            out[key] = {"gpu_mean": gp[key][0], "gpu_max": gp[key][1], "cpu_mean_max": cm, "mean_ratio": gp[key][0] / cm,
                        "max_ratio_info": gp[key][1] / max(max(c[key][1] for c in cps), 1e-300)}
        for key in g.get("quantities", ["green", "tau"]):
            if out[key]["mean_ratio"] > g["mean_ratio_max"]:
                raise CriterionFail(f"{case['gpu']}: {key} mean ratio {out[key]['mean_ratio']:.3f} > {g['mean_ratio_max']}")
    return res


def g_compare2(g, root, acc):
    """compare2.py (fail-closed) on independent-seed runs; alpha, min-bins and the exact observable count from the file."""
    A, B = os.path.join(root, g["a"]), os.path.join(root, g["b"])
    check_prov(A, acc, 0)
    check_prov(B, acc, 0, gpu_exclusive=False)
    here = os.path.dirname(os.path.abspath(__file__))
    js = os.path.join(tmpdir(), "c2.json")
    cmd = [sys.executable, os.path.join(here, "..", "compare2.py"), A, B, "--json", js,
           "--expect", ";".join(g.get("expect", [])), "--expect-m-exact", str(g["m"]),
           "--alpha", str(g["alpha"]), "--min-bins", str(g["min_bins"])]
    p = subprocess.run(cmd, capture_output=True, text=True)
    tail = p.stdout.strip().splitlines()[-1:] if p.stdout else []
    if p.returncode == 1:  # a real difference only if the JSON says so; a traceback also exits 1
        try:
            c = json.loads(read(js))
            real = c.get("verdict") == "FAIL" and (any(r.get("differs") and not r.get("under") for r in c.get("rows", []))
                                                   or any(not e.get("ok") for e in c.get("exact", [])))
        except (InputFail, ValueError):
            real = False
        if real:
            raise CriterionFail(f"compare2: an observable differs {tail}")
        raise InputFail(f"compare2 exit 1 without a FAIL verdict in its JSON (crash?): {p.stderr.strip()[-300:]}")
    if p.returncode != 0:
        raise InputFail(f"compare2 exit {p.returncode} (2 = inconclusive, 3 = input error): {tail} {p.stderr.strip()[-200:]}")
    c = json.loads(read(js))
    return [{"m": c["m"], "sum_t2": c["sum_t2"], "max_abs_t": max(abs(r["t"]) for r in c["rows"])}]


def g_log_regex(g, root, acc):
    """Checks on logs written by harness scripts (ctest, ppcmp, t_order): the stamp file must name the candidate;
    every regex must match (else INPUT_ERROR); captured integers must equal the expected ones (else FAIL);
    a nonzero exit file is a FAIL."""
    res = []
    for chk in g["checks"]:
        stamp = read(os.path.join(root, chk["stamp_file"]))
        if not acc.get("candidate_sha") or acc["candidate_sha"] not in stamp:
            raise InputFail(f"{chk['stamp_file']}: not stamped with candidate {acc.get('candidate_sha')}")
        if os.path.exists(os.path.join(root, os.path.dirname(chk["file"]), "MUTATED")) and not acc.get("_allow_mutated"):
            raise InputFail(f"{chk['file']}: mutated copy outside a control")
        bl = [l.split()[1] for l in stamp.splitlines() if l.startswith("base_sha ")]
        if bl and (not acc.get("base_sha") or len(acc["base_sha"]) < 12 or not bl[0].startswith(acc["base_sha"])):
            raise InputFail(f"{chk['stamp_file']}: base {bl[0]} is not the acceptance file's base {acc.get('base_sha')}")
        for need_s in chk.get("stamp_contains", []):
            if need_s not in stamp:
                raise InputFail(f"{chk['stamp_file']}: stamp lacks '{need_s}'")
        t = read(os.path.join(root, chk["file"]))
        m = re.search(chk["regex"], t, re.M)
        if not m:
            raise InputFail(f"{chk['file']}: no match for {chk['regex']}")
        groups = [int(x) for x in m.groups()]
        res.append({"file": chk["file"], "groups": groups})
        if "exit_file" in chk:
            ev = read(os.path.join(root, chk["exit_file"])).strip()
            if ev not in ("0", "1"):
                raise InputFail(f"{chk['exit_file']}: exit value {ev} (crash, invalid run or failed plant)")
            if ev == "1":
                raise CriterionFail(f"{chk['exit_file']}: the check failed")
        if groups != chk["equals"]:
            raise CriterionFail(f"{chk['file']}: {groups} != {chk['equals']}")
    return res


def g_files_identical(g, root, acc):
    """Two stamped harness outputs (e.g. acc_ctest.sh of base and candidate) have identical content in `file` after
    dropping lines that match `ignore_regex` (timings). Both stamps must name the candidate; INVALID exits are input
    errors; a MUTATED marker is allowed only in controls."""
    res = []
    for da, db in g["pairs"]:
        texts = []
        for d in (da, db):
            stamp = read(os.path.join(root, d, "stamp"))
            if acc["candidate_sha"] not in stamp:
                raise InputFail(f"{d}: not stamped with the candidate")
            bl = [l.split()[1] for l in stamp.splitlines() if l.startswith("base_sha ")]
            if not bl or not acc.get("base_sha") or len(acc["base_sha"]) < 12 or not bl[0].startswith(acc["base_sha"]):
                raise InputFail(f"{d}: base not the acceptance file's base {acc.get('base_sha')}")
            if os.path.exists(os.path.join(root, d, "MUTATED")) and not acc.get("_allow_mutated"):
                raise InputFail(f"{d}: mutated copy outside a control")
            ev = read(os.path.join(root, d, g.get("exit_file", "ctest.exit"))).strip()
            if ev != "0":
                raise InputFail(f"{d}: exit value {ev} (the run itself did not pass)")
            lines = [l for l in read(os.path.join(root, d, g["file"])).splitlines()
                     if not (g.get("ignore_regex") and re.search(g["ignore_regex"], l))]
            if len(lines) < g.get("min_lines", 1):
                raise InputFail(f"{d}: only {len(lines)} lines to compare")
            texts.append(lines)
        st = [{l.split(None, 1)[0]: l.split(None, 1)[1] if " " in l else "" for l in read(os.path.join(root, d, "stamp")).splitlines()}
              for d in (da, db)]
        for k in ("tool", "configure", "regex", "repeat", "env", "offload"):
            if st[0].get(k) != st[1].get(k):
                raise InputFail(f"{da}/{db}: stamps differ in {k} ({st[0].get(k)} vs {st[1].get(k)})")
        if {st[0].get("source"), st[1].get("source")} != {"base_src", "src"}:
            raise InputFail(f"{da}/{db}: must compare a base build with a candidate build")
        diff = [i for i, (x, y) in enumerate(zip(*texts)) if x != y] + ([-1] if len(texts[0]) != len(texts[1]) else [])
        res.append({"pair": [da, db], "lines": len(texts[0]), "differing": len(diff)})
        if diff:
            raise CriterionFail(f"{da} vs {db}: {len(diff)} lines differ (first index {diff[0]})")
    return res


GATES = {"files_identical": g_files_identical, "selftest": g_selftest, "faults": g_faults, "same_seed": g_same_seed, "replay": g_replay,
         "paired": g_paired, "precision": g_precision, "compare2": g_compare2, "log_regex": g_log_regex}


def run_gate(g, root, acc):
    try:
        return "PASS", GATES[g["type"]](g, root, acc), None
    except CriterionFail as e:
        return "FAIL", None, str(e)
    except InputFail as e:
        return "INPUT_ERROR", None, str(e)
    except Exception as e:  # anything unexpected: never a pass, never a detection
        return "INPUT_ERROR", None, f"{type(e).__name__}: {e}"


def resolve_control(c, acc):
    """A control names its gate ('like') and replaces only runs/pairs/cases/a/b/checks: same criterion as the gate."""
    if "like" not in c:
        raise SystemExit(f"control {c['id']}: must name its gate with 'like' (same criterion as the gate)")
    base = next((g for g in acc["gates"] if g["id"] == c["like"]), None)
    if base is None:
        raise SystemExit(f"control {c['id']}: unknown gate {c['like']}")
    allowed = {"runs", "pairs", "cases", "a", "b", "checks"}
    rp = c.get("replace", {})
    extra = set(rp) - allowed
    if extra:
        raise SystemExit(f"control {c['id']}: may replace only {sorted(allowed)}, not {sorted(extra)}")
    if "checks" in rp:  # only the files may change; regex, equals and the number of checks stay the gate's
        if len(rp["checks"]) != len(base["checks"]):
            raise SystemExit(f"control {c['id']}: must replace every check of {c['like']}")
        merged = []
        for new, old in zip(rp["checks"], base["checks"]):
            if set(new) - {"file", "stamp_file", "exit_file", "stamp_contains"}:
                raise SystemExit(f"control {c['id']}: a check may replace only file/stamp_file/exit_file/stamp_contains")
            m = dict(old, **new)  # stamp requirements can only be added (e.g. "plant --plant"), "offload MANDATORY" stays
            m["stamp_contains"] = [x for x in old.get("stamp_contains", []) if x not in ("plant no",)] + \
                [x for x in new.get("stamp_contains", []) if x not in old.get("stamp_contains", []) or x == "plant no"]
            if "plant no" in old.get("stamp_contains", []) and not ({"plant no", "plant --plant"} & set(m["stamp_contains"])):
                raise SystemExit(f"control {c['id']}: a 'plant' stamp requirement cannot be dropped (keep 'plant no' or state 'plant --plant')")
            merged.append(m)
        rp = dict(rp, checks=merged)
    if "cases" in rp and isinstance(rp["cases"], dict):  # fault controls: only the runs change, not the expected sub-tests
        if not set(rp["cases"].values()) <= set(base["cases"].values()):
            raise SystemExit(f"control {c['id']}: expected sub-tests must be the gate's")
    return dict(base, **rp, id=c["id"])


def evaluate(acc, root):
    out = {"schema": "alf-acc-verdict/2", "ticket": acc.get("ticket"), "candidate_sha": acc.get("candidate_sha"),
           "acceptance_sha256": acc.get("_sha256"), "gates": [], "controls": [], "lint": lint_errors(acc),
           "blocked": acc.get("blocked", []), "attempts": acc.get("_attempts", [])}
    ok = not out["lint"]
    for g in acc["gates"]:
        v, r, why = run_gate(g, root, acc)
        out["gates"].append({"id": g["id"], "type": g["type"], "verdict": v, "result": r, "reason": why,
                             "report_only": bool(g.get("report_only"))})
        if not g.get("report_only"):  # report-only gates (owner decision, stated in _reason) never decide the verdict
            ok &= v == "PASS"
    gate_verdict = {g["id"]: g["verdict"] for g in out["gates"]}
    for c in acc.get("controls", []):
        cg = resolve_control(c, acc)
        v, r, why = run_gate(cg, root, dict(acc, _allow_mutated=True))
        exp = c.get("expect", "FAIL")
        good = v == exp
        if gate_verdict.get(c["like"]) != "PASS":  # a control proves something only against a gate that passes
            good = False
            why = f"uninformative: gate {c['like']} is {gate_verdict.get(c['like'])} on the clean runs" + (f"; {why}" if why else "")
        out["controls"].append({"id": c["id"], "gate": cg["type"], "like": c["like"], "expect": exp, "outcome": v,
                                "as_expected": good, "reason": why})
        ok &= good
    out["overall"] = ("BLOCKED" if out["blocked"] else "PASS") if ok else "FAIL"
    return out


def write_md(v, path):
    L = [f"# Verdict ticket {v['ticket']} candidate {v['candidate_sha']}", "",
         f"acceptance file sha256 `{v['acceptance_sha256']}`", "", f"**Overall: {v['overall']}**", ""]
    L += [f"- LINT: {e}" for e in v.get("lint", [])]
    L += [f"- BLOCKED: {b.get('item')} — {b.get('reason')}" for b in v.get("blocked", [])]
    L += [f"- earlier attempt: {a}" for a in v.get("attempts", [])] + ([""] if v.get("lint") or v.get("blocked") or v.get("attempts") else [])
    L += [
         "| gate | type | verdict | detail |", "|---|---|---|---|"]
    for g in v["gates"]:
        L.append(f"| {g['id']} | {g['type']} | {g['verdict']}{' (report only)' if g.get('report_only') else ''} | {(g['reason'] or json.dumps(g['result'], default=str))[:300]} |")
    L += ["", "| control | gate | expected | outcome | as expected | reason |", "|---|---|---|---|---|---|"]
    for c in v["controls"]:
        L.append(f"| {c['id']} | {c['like']} | {c['expect']} | {c['outcome']} | {c['as_expected']} | {(c['reason'] or '')[:200]} |")
    open(path, "w").write("\n".join(L) + "\n")


# ------------------------------------------------------------------ self-check
def selfcheck():
    """Synthetic run directories through every gate type: each case must give exactly its expected outcome."""
    import h5py
    import numpy as np
    tmp = tmpdir()
    rng = np.random.default_rng(12345)
    CAND = "c0ffee000000" + "0" * 28
    builds = {"gpu": {"tool": "nvhpc", "configure": "PGI noMPI HDF5 GPU=nvidia", "gpu": True},
              "cpu": {"tool": "nvhpc", "configure": "PGI noMPI HDF5", "gpu": False},
              "igpu": {"tool": "intel", "configure": "IntelLLVM noMPI HDF5 GPU=intel", "gpu": True}}
    acc = {"candidate_sha": CAND[:12], "gpu_exclusive": True, "builds": builds, "runs": {},
           "active_lines": ["GPU tau pass (207): ACTIVE"]}
    ACT = "GPU tau pass (207): ACTIVE\n"

    def mkrun(name, data, info, log=ACT, build="gpu", env=(), allow=False, rc=0, seeds="1 2 3\n", debug=None,
              cand=CAND, binsha="b" * 64, gpu_before=(), timed_out=False, mutated=None, decl=True,
              tar_built="t" * 64, offload="MANDATORY", params="p", sel="level_zero:0", envextra=None):
        d = os.path.join(tmp, name)
        os.makedirs(d)
        with h5py.File(os.path.join(d, "data.h5"), "w") as f:
            for k, v in data.items():
                f[k] = v
        open(os.path.join(d, "info"), "w").write(info)
        open(os.path.join(d, "run.log"), "w").write(log)
        open(os.path.join(d, "seeds"), "w").write(seeds)
        b = dict(name=build, configure=builds[build]["configure"], ALF_FLAGS_EXT="", binary_sha256="b" * 64,
                 candidate={"sha": cand, "tar_sha256": "t" * 64}, source="src", tar_sha256_built=tar_built)
        hh = lambda f: hashlib.sha256(open(os.path.join(d, f), "rb").read()).hexdigest()
        envd = {"OMP_TARGET_OFFLOAD": offload} if (builds[build]["gpu"] and offload) else {}
        if build == "igpu" and sel:
            envd["ONEAPI_DEVICE_SELECTOR"] = sel
        envd.update(envextra or {})
        prov = {"schema": "alf-acc-prov/1", "binary_sha256": binsha, "build": b, "exit_code": rc, "timed_out": timed_out,
                "debug_env_present": list(debug if debug is not None else [kv.split("=")[0] for kv in env]),
                "extra_env": list(env), "gpu_before": list(gpu_before), "gpu_check_rc": 0, "load_ok": True, "env": envd,
                "inputs": {"parameters_sha256": params, "seeds_sha256": "s"},
                "outputs": {f: hh(f) for f in ("info", "data.h5", "run.log")}}
        if mutated:
            prov["mutated"] = mutated
        json.dump(prov, open(os.path.join(d, "provenance.json"), "w"))
        if decl:
            acc["runs"][name] = {"build": build, "config": "x", "env": list(env), "allow_debug": allow}
        return name

    def info(gm, gx, tm, tx):
        return (f"  Precision Green  Mean, Max :    {gm:.16E} \n   {gx:.16E}\n  Precision Phase, Max       :     0.0     \n"
                f"  Precision tau    Mean, Max :    {tm:.16E} \n   {tx:.16E}\n")

    NB = 20
    base = {"Ener_scal/obser": rng.normal(size=(NB, 1, 2)), "Ener_scal/sign": np.ones(NB),
            "Green_tau/obser": rng.normal(size=(NB, 4, 3)) * 1e-3, "Green_eq/obser": rng.normal(size=(NB, 3))}
    base["Green_tau/obser"][:, 0, 0] = 2.0  # one large entry: the others are small relative to the dataset maximum
    ok_info = info(5e-11, 1e-6, 6e-12, 5e-7)
    mkrun("A", base, ok_info)
    mkrun("B", {k: v * (1 + 1e-12) for k, v in base.items()}, ok_info, build="cpu", log="")
    small = {k: v.copy() for k, v in base.items()}
    small["Green_tau/obser"][7, 2, 2] = 1e-4  # |e|/scale = 5e-5 (above the per-entry floor 1e-6)
    small_b = {k: v.copy() for k, v in small.items()}
    small_b["Green_tau/obser"][7, 2, 2] = 1e-4 * (1 + 1e-6)  # overall metric 5e-11 (< 1e-9), per-entry 1e-6
    mkrun("A_small", small, ok_info)
    mkrun("B_small", small_b, ok_info, build="cpu", log="")
    nanb = {k: v.copy() for k, v in base.items()}
    nanb["Ener_scal/obser"][3, 0, 0] = float("nan")
    mkrun("B_nan", nanb, ok_info, build="cpu", log="")
    mkrun("B_missing", {k: v for k, v in base.items() if k != "Green_tau/obser"}, ok_info, build="cpu", log="")
    mkrun("B_seed", base, ok_info, build="cpu", log="", seeds="9 9 9\n")
    mkrun("A_inactive", base, ok_info, log="GPU tau pass (207): INACTIVE\n")
    mkrun("A_leak", base, ok_info, debug=["ALF_GPU_NS"])
    mkrun("A_crash", base, ok_info, rc=137)
    mkrun("A_timeout", base, ok_info, rc=124, timed_out=True)
    mkrun("A_busy", base, ok_info, gpu_before=["12345, other, 300 MiB"])
    mkrun("A_othercand", base, ok_info, cand="deadbeef0000" + "0" * 28)
    mkrun("A_binsha", base, ok_info, binsha="c" * 64)
    mkrun("A_mut", base, ok_info, mutated="scale")
    mkrun("A_undecl", base, ok_info, decl=False)
    mkrun("A_twin", base, ok_info)
    mkrun("A_tar", base, ok_info, tar_built="x" * 64)
    mkrun("A_nooffload", base, ok_info, offload=None)
    mkrun("B_params", base, ok_info, build="cpu", log="", params="q")
    mkrun("A_tamper", base, ok_info)
    mkrun("A_in", base, ok_info, build="igpu")
    mkrun("A_in_nosel", base, ok_info, build="igpu", sel=None)
    mkrun("A_in_sel1", base, ok_info, build="igpu", sel="level_zero:1")
    mkrun("A_in_ze", base, ok_info, build="igpu", envextra={"ZE_AFFINITY_MASK": "0"})
    mkrun("A_laleak", base, ok_info, debug=["ALF_LA_TEST_X"])
    mkrun("A_ladecl", base, ok_info, env=["ALF_LA_TEST_X=1"])
    mkrun("A_laenv", base, ok_info, debug=[], envextra={"ALF_LA_TEST_X": "1"})
    open(os.path.join(tmp, "A_tamper", "run.log"), "a").write("edited later\n")
    nan_info = info(float("nan"), 1e-6, 6e-12, 5e-7)
    assert "NAN" in nan_info.upper()
    mkrun("A_nanprec", base, nan_info)
    mkrun("A_2x", base, info(1e-10, 1e-6, 6e-12, 5e-7))
    mkrun("A_2x_inactive", base, info(5e-11, 1e-6, 6e-12, 5e-7), log="GPU tau pass (207): INACTIVE\n")

    ss = {"type": "same_seed", "id": "ss", "n_datasets": 4, "n_tau_datasets": 1, "no_bin_above": 1e-4,
          "overall_max": 1e-9, "tau_max": 1e-9}
    sse = dict(ss, per_entry_max=1e-8)
    rp = {"type": "replay", "id": "rp", "n_datasets": 4}
    pr = {"type": "precision", "id": "pr", "mean_ratio_max": 1.5}
    cases = [("same_seed clean", dict(ss, pairs=[["A", "B", NB]]), "PASS"),
             ("same_seed small entry, no per-entry", dict(ss, pairs=[["A_small", "B_small", NB]]), "PASS"),
             ("same_seed small entry, per-entry", dict(sse, pairs=[["A_small", "B_small", NB]]), "FAIL"),
             ("same_seed NaN", dict(ss, pairs=[["A", "B_nan", NB]]), "INPUT_ERROR"),
             ("same_seed missing dataset", dict(ss, pairs=[["A", "B_missing", NB]]), "INPUT_ERROR"),
             ("same_seed other seeds", dict(ss, pairs=[["A", "B_seed", NB]]), "INPUT_ERROR"),
             ("same_seed wrong nbin", dict(ss, pairs=[["A", "B", NB + 1]]), "INPUT_ERROR"),
             ("same_seed GPU inactive", dict(ss, pairs=[["A_inactive", "B", NB]]), "INPUT_ERROR"),
             ("same_seed leaked debug var", dict(ss, pairs=[["A_leak", "B", NB]]), "INPUT_ERROR"),
             ("same_seed crash 137", dict(ss, pairs=[["A_crash", "B", NB]]), "INPUT_ERROR"),
             ("same_seed timeout", dict(ss, pairs=[["A_timeout", "B", NB]]), "INPUT_ERROR"),
             ("same_seed GPU busy", dict(ss, pairs=[["A_busy", "B", NB]]), "INPUT_ERROR"),
             ("same_seed other candidate", dict(ss, pairs=[["A_othercand", "B", NB]]), "INPUT_ERROR"),
             ("same_seed binary sha", dict(ss, pairs=[["A_binsha", "B", NB]]), "INPUT_ERROR"),
             ("same_seed mutated in a gate", dict(ss, pairs=[["A_mut", "B", NB]]), "INPUT_ERROR"),
             ("same_seed undeclared run", dict(ss, pairs=[["A_undecl", "B", NB]]), "INPUT_ERROR"),
             ("same_seed missing run dir", dict(ss, pairs=[["A_nowhere", "B", NB]]), "INPUT_ERROR"),
             ("same_seed binary not from tarball", dict(ss, pairs=[["A_tar", "B", NB]]), "INPUT_ERROR"),
             ("same_seed no MANDATORY offload", dict(ss, pairs=[["A_nooffload", "B", NB]]), "INPUT_ERROR"),
             ("same_seed parameters differ", dict(ss, pairs=[["A", "B_params", NB]]), "INPUT_ERROR"),
             ("same_seed output edited later", dict(ss, pairs=[["A_tamper", "B", NB]]), "INPUT_ERROR"),
             ("same_seed GPU=intel pinned", dict(ss, pairs=[["A_in", "B", NB]]), "PASS"),
             ("same_seed GPU=intel no selector", dict(ss, pairs=[["A_in_nosel", "B", NB]]), "INPUT_ERROR"),
             ("same_seed GPU=intel other selector", dict(ss, pairs=[["A_in_sel1", "B", NB]]), "INPUT_ERROR"),
             ("same_seed GPU=intel ZE_AFFINITY_MASK", dict(ss, pairs=[["A_in_ze", "B", NB]]), "INPUT_ERROR"),
             ("same_seed leaked ALF_LA_TEST_", dict(ss, pairs=[["A_laleak", "B", NB]]), "INPUT_ERROR"),
             ("same_seed ALF_LA_TEST_ no allow_debug", dict(ss, pairs=[["A_ladecl", "B", NB]]), "INPUT_ERROR"),
             ("same_seed ALF_LA_TEST_ only in env", dict(ss, pairs=[["A_laenv", "B", NB]]), "INPUT_ERROR"),
             ("replay identical", dict(rp, pairs=[["A", "A_twin"]]), "PASS"),
             ("replay 1e-12", dict(rp, pairs=[["A", "B"]], gpu_exclusive=False), "FAIL"),
             ("replay Precision Max only", dict(rp, pairs=[["A", mkrun("A_pmax", base, info(5e-11, 1.0000001e-6, 6e-12, 5e-7))]]), "FAIL"),
             ("replay Precision Max unwrapped", dict(rp, pairs=[[mkrun("A_unw", base, info(5e-11, 1e-6, 6e-12, 5e-7).replace(" \n   ", "  ")),
                                                                 mkrun("A_unw2", base, info(5e-11, 1.0000001e-6, 6e-12, 5e-7).replace(" \n   ", "  "))]]), "FAIL"),
             ("precision clean", dict(pr, cases=[{"gpu": "A", "cpu": ["B"]}]), "PASS"),
             ("precision NaN", dict(pr, cases=[{"gpu": "A_nanprec", "cpu": ["B"]}]), "INPUT_ERROR"),
             ("precision 2x mean", dict(pr, cases=[{"gpu": "A_2x", "cpu": ["B"]}]), "FAIL"),
             ("precision GPU inactive", dict(pr, cases=[{"gpu": "A_2x_inactive", "cpu": ["B"]}]), "INPUT_ERROR")]

    def taulog(n_pass, sd_dev, bias, nan_at=None):
        L = [ACT.strip()]
        for ps in range(1, n_pass + 1):
            for nst in range(1, 5):
                for nfe in (1, 2):
                    h = 10 ** rng.normal(-7.5, 0.5)
                    pa, pb = h * 10 ** rng.normal(0, 0.1), h * 10 ** rng.normal(0, 0.1)
                    dv = h * 10 ** (rng.normal(0, sd_dev) + bias)
                    if nan_at == (ps, nst, nfe):
                        dv = float("nan")
                    L.append(f" TAUDIAG pt 0 pass {ps} nst {nst} nfe {nfe} br 1 dev 1E-12 {dv:.3E} 1E-12 1E-12 | "
                             f"hostfull 1E-12 {h:.3E} 1E-12 1E-12 | P {h:.3E} {pa:.3E} {pb:.3E} {pa:.3E} {pb:.3E} |")
        return "\n".join(L) + "\n"
    E = ["ALF_GPU_TAUDIAG=2"]
    pg = {"type": "paired", "id": "pd", "n_points": 400, "n_used": 200, "from_pass": 26, "abs_mean_max": 0.05,
          "sd_ratio_max": 1.25, "count_thresholds": [1e-7], "count_k": 2.0, "count_min_events": 10}
    cnt_lines = []
    for line in taulog(50, 0.1, 0.0).splitlines():  # count-only failure: dev x30 only for passes < from_pass
        t = line.split()
        if line.lstrip().startswith("TAUDIAG pt ") and int(t[t.index("pass") + 1]) < 26:
            i = t.index("dev") + 2
            t[i] = f"{float(t[i]) * 30:.3E}"
            line = " " + " ".join(t)
        cnt_lines.append(line)
    cases += [("paired clean", dict(pg, runs=[mkrun("P_ok", base, ok_info, taulog(50, 0.1, 0.0), env=E, allow=True)]), "PASS"),
              ("paired variance 1.6x ref", dict(pg, runs=[mkrun("P_var", base, ok_info, taulog(50, 0.225, 0.0), env=E, allow=True)]), "FAIL"),
              ("paired bias x1.3", dict(pg, runs=[mkrun("P_bias", base, ok_info, taulog(50, 0.1, math.log10(1.3)), env=E, allow=True)]), "FAIL"),
              ("paired counts only", dict(pg, runs=[mkrun("P_cnt", base, ok_info, "\n".join(cnt_lines) + "\n", env=E, allow=True)]), "FAIL"),
              ("paired NaN", dict(pg, runs=[mkrun("P_nan", base, ok_info, taulog(50, 0.1, 0.0, nan_at=(30, 2, 1)), env=E, allow=True)]), "INPUT_ERROR"),
              ("paired truncated", dict(pg, runs=[mkrun("P_short", base, ok_info, taulog(40, 0.1, 0.0), env=E, allow=True)]), "INPUT_ERROR"),
              ("paired wrong from_pass", dict(pg, from_pass=27, runs=["P_ok"]), "INPUT_ERROR"),
              ("paired undeclared debug var", dict(pg, runs=[mkrun("P_leak", base, ok_info, taulog(50, 0.1, 0.0), env=E, allow=True,
                                                                  debug=["ALF_GPU_TAUDIAG", "ALF_GPU_TAUSWAP"])]), "INPUT_ERROR"),
              ("paired debug without allow", dict(pg, runs=[mkrun("P_noallow", base, ok_info, taulog(50, 0.1, 0.0), env=E, allow=False)]), "INPUT_ERROR"),
              ("paired with 212 env and debug lines", dict(pg, runs=[mkrun("P_envline", base, ok_info, " GPU debug env: ALF_GPU_DEBUG=1, ALF_GPU_TAUDIAG=2\n"
                                                                            " GPU tau pass: DEBUG ALF_GPU_TAUDIAG=1 active (host oracle chains)\n" + taulog(50, 0.1, 0.0), env=E, allow=True)]), "PASS"),
              ("paired malformed data line", dict(pg, runs=[mkrun("P_garbled", base, ok_info, taulog(50, 0.1, 0.0) + " TAUDIAG pt 9 pass x\n", env=E, allow=True)]), "INPUT_ERROR")]
    STL = ACT + ("    G00  3.478E-14 (worst slice 11), GT0  9.805E-14 (40), GTT  4.954E-12 (40), G0T  3.122E-13 (40)\n"
                 " GPU tau self-test CGR2_2 alone: natural calls with D1(1) > D2(1): 4, else: 4, total 8\n"
                 " GPU tau self-test precision reduction vs COMPARE: max rel diff (max and mean)  0.000E+00 (tol  1.0E-12\n"
                 " GPU tau self-test PASSED\n")
    S = ["ALF_GPU_SELFTEST=4"]
    nv = numbers_after("GPU stab self-test CGR at 8 stabilization points x flavours: max rel diff G NVAR=1  1.0E-12, NVAR=2  2.0E-12; phase NVAR=1  3.0E-13, NVAR=2  4.0E-13 (tol  1.0E-10",
                       r"max rel diff G")
    assert nv == [1e-12, 2e-12, 3e-13, 4e-13], nv
    st = {"type": "selftest", "id": "st", "require_lines": ["GPU tau self-test PASSED"],
          "number_checks": [{"name": "pass", "line": r"^\s+G00\s", "after": r"^", "max": 1e-10, "min_count": 4},
                            {"name": "prec", "line": "precision reduction vs COMPARE", "after": r"\(max and mean\)", "max": 1e-12}],
          "int_checks": [{"name": "branches", "regex": r"D1\(1\) > D2\(1\): (\d+), else: (\d+)", "min": 1}]}
    cases += [("selftest clean", dict(st, runs=[mkrun("S_ok", base, ok_info, STL, env=S, allow=True)]), "PASS"),
              ("selftest crash", dict(st, runs=[mkrun("S_crash", base, ok_info, STL.replace(" GPU tau self-test PASSED\n", ""), env=S, allow=True, rc=139)]), "INPUT_ERROR"),
              ("selftest above tol", dict(st, runs=[mkrun("S_big", base, ok_info, STL.replace("4.954E-12", "4.954E-09"), env=S, allow=True)]), "FAIL"),
              ("selftest one branch", dict(st, runs=[mkrun("S_br", base, ok_info, STL.replace("else: 4", "else: 0"), env=S, allow=True)]), "FAIL"),
              ("selftest missing line", dict(st, runs=[mkrun("S_noline", base, ok_info, STL.replace("precision reduction", "xx"), env=S, allow=True)]), "INPUT_ERROR"),
              ("selftest FAILED rc 1", dict(st, runs=[mkrun("S_fail", base, ok_info, ACT + " GPU self-test FAILED: tau: pass\n", env=S, allow=True, rc=1)]), "FAIL"),
              ("selftest setup error rc 1", dict(st, runs=[mkrun("S_setup", base, ok_info, ACT + " gpu_state: device allocation FAILED\n", env=S, allow=True, rc=1)]), "INPUT_ERROR")]
    F = ["ALF_GPU_SELFTEST=4", "ALF_GPU_FAULT_TAU=perm"]
    ft = {"type": "faults", "id": "ft"}
    cases += [("fault detected", dict(ft, cases={mkrun("F_ok", base, ok_info, ACT + " GPU self-test FAILED: tau: pass\n", rc=1, env=F, allow=True): "tau: pass"}), "PASS"),
              ("fault = crash", dict(ft, cases={mkrun("F_crash", base, ok_info, ACT + "Segmentation fault\n", rc=139, env=F, allow=True): "tau: pass"}), "INPUT_ERROR"),
              ("fault not detected", dict(ft, cases={mkrun("F_miss", base, ok_info, STL, rc=0, env=F, allow=True): "tau: pass"}), "FAIL"),
              ("fault wrong test", dict(ft, cases={mkrun("F_other", base, ok_info, ACT + " GPU self-test FAILED: hopping\n", rc=1, env=F, allow=True): "tau: pass"}), "FAIL"),
              ("fault = setup error rc 1", dict(ft, cases={mkrun("F_setup", base, ok_info, ACT + " setup FAILED\n", rc=1, env=F, allow=True): "tau: pass"}), "INPUT_ERROR")]
    os.makedirs(os.path.join(tmp, "g"))
    open(os.path.join(tmp, "g", "stamp"), "w").write(CAND + "\n")
    open(os.path.join(tmp, "g", "stamp_old"), "w").write("deadbeef0000\n")
    open(os.path.join(tmp, "g", "ctest.log"), "w").write("100% tests passed, 0 tests failed out of 52\n")
    open(os.path.join(tmp, "g", "ctest.exit"), "w").write("0\n")
    open(os.path.join(tmp, "g", "bad.exit"), "w").write("99\n")
    lg = {"type": "log_regex", "id": "lg"}
    ck = {"file": "g/ctest.log", "regex": r"(\d+)% tests passed, (\d+) tests failed out of (\d+)", "exit_file": "g/ctest.exit", "stamp_file": "g/stamp"}
    cases += [("log_regex clean", dict(lg, checks=[dict(ck, equals=[100, 0, 52])]), "PASS"),
              ("log_regex count", dict(lg, checks=[dict(ck, equals=[100, 0, 53])]), "FAIL"),
              ("log_regex stale stamp", dict(lg, checks=[dict(ck, equals=[100, 0, 52], stamp_file="g/stamp_old")]), "INPUT_ERROR"),
              ("log_regex exit 99 (invalid run)", dict(lg, checks=[dict(ck, equals=[100, 0, 52], exit_file="g/bad.exit")]), "INPUT_ERROR"),
              ("log_regex stamp lacks marker", dict(lg, checks=[dict(ck, equals=[100, 0, 52], stamp_contains=["plant no"])]), "INPUT_ERROR")]
    bad = 0
    for name, g, exp in cases:
        v, r, why = run_gate(g, tmp, acc)
        bad += v != exp
        print(f"{'ok ' if v == exp else 'BAD'} {name:36s} expected {exp:11s} got {v:11s} {('(' + why + ')') if why else ''}"[:220])
    acc2 = dict(acc, gates=[dict(ss, pairs=[["A", "B", NB]])],
                controls=[{"id": "c_missing", "like": "ss", "replace": {"pairs": [["A_nowhere", "B", NB]]}, "expect": "FAIL"},
                          {"id": "c_real", "like": "ss", "replace": {"pairs": [["A_small", "B_small", NB]]}, "expect": "PASS"}])
    v = evaluate(acc2, tmp)
    sem = [(c["id"], c["outcome"], c["as_expected"]) for c in v["controls"]]
    good = sem == [("c_missing", "INPUT_ERROR", False), ("c_real", "PASS", True)] and v["overall"] == "FAIL"
    acc3 = dict(acc, gates=[dict(sse, pairs=[["A_small", "B_small", NB]])],  # the gate itself FAILs on clean runs
                controls=[{"id": "c_uninf", "like": "ss", "replace": {"pairs": [["A_small", "B_small", NB]]}, "expect": "FAIL"}])
    acc3["runs"] = acc["runs"]
    v3 = evaluate(dict(acc3, controls=acc3["controls"] + [{"id": "c_ok_dummy", "like": "ss", "replace": {"pairs": [["A", "B_small", NB]]}, "expect": "PASS"}]), tmp)
    try:
        resolve_control({"id": "c_bad", "like": "lg", "replace": {"checks": [dict(ck, equals=[-1])]}}, dict(acc, gates=[dict(lg, checks=[dict(ck, equals=[100, 0, 52])])]))
        good = False
    except SystemExit:
        pass
    good &= [(c["outcome"], c["as_expected"]) for c in v3["controls"]][:1] == [("FAIL", False)]
    bad += not good
    print(f"{'ok ' if good else 'BAD'} control semantics (missing runs are not a kill) {sem}")
    # lint: report_only needs an owner approval with a decision id found in DECISIONS.md; controls cannot swap in
    # unrelated runs
    base_acc = {"ticket": "zz9", "builds": builds, "configs": {"x": {}, "y": {}},
                "runs": {"r1": {"build": "gpu", "config": "x"}, "r2": {"build": "cpu", "config": "x"}, "r3": {"build": "gpu", "config": "y"}},
                "gates": [dict(ss, pairs=[["r1", "r2", 1]]), dict(ss, id="ro", report_only=True, pairs=[["r1", "r2", 1]], _reason="decision")],
                "controls": [{"id": "c1", "like": "ss", "replace": {"pairs": [["r3", "r2", 1]]}, "expect": "FAIL"}], "_status": "DRAFT"}
    le = lint_errors(base_acc)
    ok_l = any("report_only needs an owner approval" in e for e in le) and any("matches no run of gate" in e for e in le)
    base_acc["approvals"] = [{"who": "owner", "decision": "99zz", "covers": ["ro"]}]
    ok_l &= any("not found in docs/gpu/DECISIONS.md" in e for e in lint_errors(base_acc))
    bad += not ok_l
    print(f"{'ok ' if ok_l else 'BAD'} lint: report_only without a real owner decision, unrelated control runs")
    # lint: ALF_LA_TEST_* is a debug variable (only in faults gates and controls); no declared selector on GPU=intel
    la_acc = {"ticket": "zz9", "builds": builds, "configs": {"x": {}},
              "runs": {"q1": {"build": "gpu", "config": "x", "env": ["ALF_LA_TEST_X=1"], "allow_debug": True},
                       "q2": {"build": "cpu", "config": "x"},
                       "q3": {"build": "igpu", "config": "x", "env": ["ONEAPI_DEVICE_SELECTOR=level_zero:1"]}},
              "gates": [dict(ss, pairs=[["q1", "q2", 1]])], "controls": [], "_status": "DRAFT"}
    le = lint_errors(la_acc)
    ok_x = any("run q1 carries fault/numerics variables" in e for e in le) and any(e.startswith("run q3:") and "ONEAPI_DEVICE_SELECTOR" in e for e in le)
    la_acc["gates"] = [{"type": "faults", "id": "ft", "cases": {"q1": "x"}}]
    ok_x &= not any("carries fault/numerics" in e for e in lint_errors(la_acc))
    bad += not ok_x
    print(f"{'ok ' if ok_x else 'BAD'} lint: ALF_LA_TEST_ in a gate run, faults gate exempt, selector declared on GPU=intel")
    shutil.rmtree(tmp)
    print("selfcheck:", "PASS" if not bad else f"FAIL ({bad} wrong)")
    return 0 if not bad else 1


GATE_DEBUG_OK = ("ALF_GPU_SELFTEST", "ALF_GPU_EXIT", "ALF_GPU_TAUDIAG", "ALF_GPU_TAUDIAG_FROM", "ALF_GPU_DEBUG")
SELFTEST_ONLY = ("ALF_GPU_UPDTEST",)  # restores no acceptance counters (Wrapgr_mod): only in self-test gates


def lint_errors(acc):
    """Static checks of an acceptance file (enforced by evaluate() and the driver)."""
    errs = []
    runs, builds, configs = acc.get("runs", {}), acc.get("builds", {}), acc.get("configs", {})
    draft = "DRAFT" in acc.get("_status", "")
    used_builds = set()
    for rid, r in runs.items():
        if r.get("build") not in builds:
            errs.append(f"run {rid}: unknown build {r.get('build')}")
        used_builds.add(r.get("build"))
        if r.get("config") not in configs:
            errs.append(f"run {rid}: unknown config {r.get('config')}")
        if "gpu_active" in r:
            errs.append(f"run {rid}: 'gpu_active' is not allowed (every GPU run must show the ACTIVE lines)")
        if "GPU=intel" in builds.get(r.get("build"), {}).get("configure", "").split():
            pin = [kv for kv in r.get("env", []) if kv.split("=")[0] in ("ONEAPI_DEVICE_SELECTOR", "ZE_AFFINITY_MASK")]
            if pin:
                errs.append(f"run {rid}: {pin} declared for a GPU=intel build (acc_run.sh sets the selector and refuses these)")
    for b in builds:
        if b not in used_builds:
            errs.append(f"build {b}: declared but unused")
    if not draft:
        for c, v in configs.items():
            if "parameters_sha256" not in v:
                errs.append(f"config {c}: parameters_sha256 missing (required when frozen)")
    mut_dst = {m["dst"] for m in acc.get("mutations", [])}

    dirs = {t["id"] for t in acc.get("ctests", [])} | {m["dst"] for m in acc.get("log_mutations", [])}
    if acc.get("gpunone"):
        dirs |= {"gpunone", "gpunone_plant"}

    def dir_refs(g):
        out = [x for p in g.get("pairs", []) for x in p[:2]] if g["type"] == "files_identical" else []
        out += [c[k].split("/")[0] for c in g.get("checks", []) for k in ("file", "stamp_file", "exit_file") if k in c]
        return out

    def refs(g):
        if g["type"] in ("log_regex", "files_identical"):
            return []
        cs = g.get("cases", {})
        out = list(g.get("runs", [])) + [x for p in g.get("pairs", []) for x in p[:2]] + (list(cs.keys()) if isinstance(cs, dict) else [])
        out += [g[k] for k in ("a", "b") if k in g]
        for c in cs if isinstance(cs, list) else []:
            out += [c["gpu"]] + c["cpu"]
        return [x for x in out if isinstance(x, str)]
    gates = {g["id"]: g for g in acc.get("gates", [])}
    for g in acc.get("gates", []):
        for x in refs(g):
            if x not in runs:
                errs.append(f"gate {g['id']}: undeclared run {x}")
            elif x in mut_dst:
                errs.append(f"gate {g['id']}: uses the mutated copy {x}")
            else:
                ok_vars = GATE_DEBUG_OK + (SELFTEST_ONLY if g["type"] == "selftest" else ())
                bad = [kv for kv in runs[x].get("env", []) if kv.split("=")[0].startswith(DEBUG_PREFIXES)
                       and kv.split("=")[0] not in ok_vars and g["type"] != "faults"]
                if bad:
                    errs.append(f"gate {g['id']}: run {x} carries fault/numerics variables {bad} (allowed only in controls)")
        for p in g.get("pairs", []):
            if p[0] == p[1]:
                errs.append(f"gate {g['id']}: pair compares {p[0]} with itself")
            if g["type"] == "same_seed" and p[0] in runs and p[1] in runs:
                ga = builds.get(runs[p[0]]["build"], {}).get("gpu", True)
                gb = builds.get(runs[p[1]]["build"], {}).get("gpu", True)
                if not (ga and not gb):
                    errs.append(f"gate {g['id']}: same-seed pair {p[:2]} is not GPU build vs CPU build")
    for g in acc.get("gates", []):
        for x in dir_refs(g):
            if x not in dirs:
                errs.append(f"gate {g['id']}: output directory {x} is produced by no ctest, log mutation or gpunone step")
    fail_ctl = set()
    for c in acc.get("controls", []):
        if c.get("like") not in gates:
            errs.append(f"control {c['id']}: unknown gate {c.get('like')}")
            continue
        if c.get("expect", "FAIL") == "FAIL":
            fail_ctl.add(c["like"])
        try:
            cg = resolve_control(c, acc)
        except SystemExit as e:
            errs.append(str(e))
            continue
        for x in refs(cg):
            if x not in runs:
                errs.append(f"control {c['id']}: undeclared run {x}")
    dec_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "docs", "gpu", "DECISIONS.md")
    decisions = open(dec_path).read() if os.path.isfile(dec_path) else ""
    for gid, g in gates.items():
        if g.get("report_only"):
            ap = [a for a in acc.get("approvals", []) if a.get("who") == "owner" and gid in a.get("covers", []) and a.get("decision")]
            row = [l for l in decisions.splitlines() if ap and l.startswith(f"| {ap[0]['decision']} |")]
            if not ap:
                errs.append(f"gate {gid}: report_only needs an owner approval covering it, with a decision id")
            elif not row or str(acc.get("ticket")) not in row[0]:
                errs.append(f"gate {gid}: decision {ap[0]['decision']} not found in docs/gpu/DECISIONS.md, or it does not name ticket {acc.get('ticket')}")
            elif not ap[0].get("quote") or ap[0]["quote"] not in row[0] or gid not in row[0]:
                errs.append(f"gate {gid}: decision {ap[0]['decision']} must quote the owner verbatim ({ap[0].get('quote')!r}) and name the gate id")
            continue
        if gid not in fail_ctl:
            errs.append(f"gate {gid}: no positive control with expect FAIL")
    null_ctl = set()
    for c in acc.get("controls", []):
        cg = gates.get(c.get("like"))
        if cg and cg["type"] == "compare2" and c.get("expect") == "PASS" and c.get("null"):
            r = c.get("replace", {})
            ra, rb = r.get("a", cg.get("a")), r.get("b", cg.get("b"))
            if ra in runs and rb in runs and not builds.get(runs[ra]["build"], {}).get("gpu", True) \
                    and not builds.get(runs[rb]["build"], {}).get("gpu", True) and runs[ra].get("seeds") != runs[rb].get("seeds"):
                null_ctl.add(c["like"])
    for gid, g in gates.items():
        if g["type"] == "compare2" and gid not in null_ctl:
            errs.append(f"gate {gid}: compare2 needs a null control (\"null\": true, expect PASS, CPU vs CPU on two seed sets)")
    for m in acc.get("mutations", []):
        if m["src"] not in runs or m["dst"] not in runs:
            errs.append(f"mutation {m['dst']}: source or destination not declared")
    # controls may only swap in runs/outputs of the same kind as the gate's (S3, M2)
    mut = {m["dst"]: m["src"] for m in acc.get("mutations", [])} | {m["dst"]: m["src"] for m in acc.get("log_mutations", [])}
    ct = {t["id"]: t for t in acc.get("ctests", [])}

    def same_run(x, y):
        if x == y or mut.get(x) == y:
            return True
        rx, ry = runs.get(x), runs.get(y)
        return bool(rx and ry and rx["config"] == ry["config"] and rx["build"] == ry["build"]
                    and rx.get("seeds", "dev") == ry.get("seeds", "dev"))

    def same_dir(x, y):
        if x == y or mut.get(x) == y:
            return True
        tx, ty = ct.get(x), ct.get(y)
        return bool(tx and ty and all(tx.get(k) == ty.get(k) for k in ("tool", "configure", "regex", "repeat", "base", "gpu")))
    for c in acc.get("controls", []):
        g = gates.get(c.get("like"))
        if not g:
            continue
        r = c.get("replace", {})
        pairs_old, pairs_new = [], []
        for k in ("runs", "a", "b"):
            if k in r:
                o, n = g.get(k), r[k]
                o, n = (o if isinstance(o, list) else [o]), (n if isinstance(n, list) else [n])
                pairs_old += o; pairs_new += n
        for k in ("pairs", "runs"):
            if k in r and len(r[k]) > len(g.get(k, [])):
                errs.append(f"control {c['id']}: more {k} than its gate")
        if isinstance(r.get("cases"), (list, dict)) and len(r["cases"]) > len(g.get("cases", [])):
            errs.append(f"control {c['id']}: more cases than its gate")
        if "pairs" in r:
            for po, pn in zip(g["pairs"], r["pairs"]):  # pairs are positional: each new pair matches the gate's pair
                if not (same_run(pn[0], po[0]) and same_run(pn[1], po[1])) and not (c.get("null") and g["type"] == "compare2"):
                    errs.append(f"control {c['id']}: pair {pn[:2]} does not correspond to the gate's pair {po[:2]}")
                pairs_old += po[:2]; pairs_new += pn[:2]
        if "cases" in r and isinstance(r["cases"], dict):
            pairs_old += list(g["cases"].keys()); pairs_new += list(r["cases"].keys())
        if "cases" in r and isinstance(r["cases"], list):
            for co, cn in zip(g["cases"], r["cases"]):
                pairs_old += [co["gpu"]] + co["cpu"]; pairs_new += [cn["gpu"]] + cn["cpu"]
        for n in pairs_new:  # each swapped-in run must be of the kind of some run of the gate, or a mutation of one
            if not any(same_run(n, o) for o in pairs_old + [x for p in g.get("pairs", []) for x in p[:2]] + list(g.get("runs", []))
                       + ([g[k] for k in ("a", "b") if k in g]) + (list(g["cases"].keys()) if isinstance(g.get("cases"), dict) else [])
                       + ([x for cc in g["cases"] for x in [cc["gpu"]] + cc["cpu"]] if isinstance(g.get("cases"), list) else [])) \
                    and not (c.get("null") and g["type"] == "compare2"):
                errs.append(f"control {c['id']}: {n} matches no run of gate {g['id']} in config and build, and is no mutation of one")
        for chk_o, chk_n in zip(g.get("checks", []), r.get("checks", [])):
            do, dn = chk_o["file"].split("/")[0], chk_n.get("file", chk_o["file"]).split("/")[0]
            if dn in ct and do in ct and not same_dir(dn, do):
                errs.append(f"control {c['id']}: ctest {dn} differs from {do} in more than its environment")
            lm = {m["dst"]: m["src"] for m in acc.get("log_mutations", [])}
            if dn in lm and lm[dn] != do:
                errs.append(f"control {c['id']}: {dn} is a log mutation of {lm[dn]}, not of the gate's {do}")
    for g in acc.get("gates", []):
        for d in dir_refs(g):
            if ct.get(d, {}).get("control") or ct.get(d, {}).get("env") and not ct.get(d, {}).get("control"):
                errs.append(f"gate {g['id']}: uses ctest output {d} that is a control or carries an environment")
    return errs


def lint(acc):
    errs = lint_errors(acc)
    for e in errs:
        print("LINT:", e)
    print("lint:", "OK" if not errs else f"{len(errs)} problems")
    return 1 if errs else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("acceptance", nargs="?")
    ap.add_argument("runroot", nargs="?")
    ap.add_argument("--out", default=None)
    ap.add_argument("--candidate", default=None, help="candidate sha (replaces \"self\" in the acceptance file)")
    ap.add_argument("--seeds-record", default=None, help="the driver's seeds record (sha256 per seed set/config)")
    ap.add_argument("--attempts", default=None, help="the ticket's attempts log (earlier T2 runs are listed)")
    ap.add_argument("--selfcheck", action="store_true")
    ap.add_argument("--lint", action="store_true", help="static checks of the acceptance file only")
    a = ap.parse_args()
    if a.selfcheck:
        return selfcheck()
    if a.lint and a.acceptance:
        return lint(json.load(open(a.acceptance)))
    if not (a.acceptance and a.runroot):
        ap.print_usage()
        return 2
    raw = open(a.acceptance, "rb").read()
    acc = json.loads(raw)
    acc["_sha256"] = hashlib.sha256(raw).hexdigest()
    if acc.get("candidate_sha") in (None, "", "self"):
        if not a.candidate:
            print("the acceptance file says candidate \"self\": give --candidate SHA")
            return 2
        acc["candidate_sha"] = a.candidate
    uses_acc = any(r.get("seeds", "dev").startswith("acc") for r in acc.get("runs", {}).values())
    if a.seeds_record:
        acc["_seeds_record"] = json.load(open(a.seeds_record))
    elif uses_acc:
        print("runs use acceptance seeds: --seeds-record is required")
        return 2
    if a.attempts and os.path.isfile(a.attempts):
        acc["_attempts"] = [l.strip() for l in open(a.attempts) if l.strip()]
    out = a.out or a.runroot
    for f in ("verdict.json", "verdict.md"):
        if os.path.exists(os.path.join(out, f)):
            print(f"refusing to overwrite {os.path.join(out, f)}")
            return 2
    v = evaluate(acc, a.runroot)
    json.dump(v, open(os.path.join(out, "verdict.json"), "w"), indent=1, default=str)
    write_md(v, os.path.join(out, "verdict.md"))
    print(open(os.path.join(out, "verdict.md")).read())
    return 0 if v["overall"] == "PASS" else (3 if v["overall"] == "BLOCKED" else 1)


if __name__ == "__main__":
    sys.exit(main())
