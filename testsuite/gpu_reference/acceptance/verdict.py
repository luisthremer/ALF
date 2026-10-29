#!/usr/bin/env python3
"""Acceptance verdicts from raw run outputs (ALF GPU port, standard of 2026-10-29).

usage: verdict.py ACCEPTANCE.json RUNROOT [--out DIR]     compute every gate and control
       verdict.py --selfcheck                              positive controls of this script (synthetic data, CPU only)

RUNROOT holds one directory per run (as written by acc_run.sh: provenance.json, info, run.log, data.h5, ...).
The acceptance file names the run directories relative to RUNROOT. Writes verdict.json and verdict.md.

Fail-closed rules, applied by every gate (any of them gives FAIL, never a silent skip):
  a missing run directory, file, line or dataset; a nonzero exit code or a timeout (unless the gate expects one);
  any NaN/Inf in a parsed number; a count that differs from the acceptance file; a regex that does not match;
  a provenance mismatch (candidate sha, binary sha256 vs its build record, debug variables in a run that does not
  declare them, another process on the GPU at start).
Exit code: 0 only if every gate PASSes and every control gives its expected outcome; 1 otherwise; 2 on a usage error.
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


class Fail(Exception):
    pass


def num(s):
    v = float(s.replace("D", "E").replace("d", "e"))
    if not math.isfinite(v):
        raise Fail(f"non-finite number {s}")
    return v


def need(path):
    if not os.path.exists(path):
        raise Fail(f"missing {path}")
    return path


def read(path):
    with open(need(path), errors="replace") as f:
        return f.read()


# ------------------------------------------------------------------ provenance
def check_prov(rd, acc, expect_rc=0, allow_debug=False, gpu_exclusive=None):
    """Common checks on a run directory's provenance.json. Returns the provenance dict."""
    p = json.loads(read(os.path.join(rd, "provenance.json")))
    if p.get("schema") != "alf-acc-prov/1":
        raise Fail(f"{rd}: unknown provenance schema")
    rc = p.get("exit_code")
    if p.get("timed_out"):
        raise Fail(f"{rd}: timed out")
    if expect_rc is not None and rc != expect_rc:
        raise Fail(f"{rd}: exit code {rc}, expected {expect_rc}")
    cand = acc.get("candidate_sha")
    got = p.get("build", {}).get("candidate", {})
    src = p.get("build", {}).get("source", "src")
    want = acc.get("base_sha") if src == "base_src" else cand
    have = got.get("base_sha") if src == "base_src" else got.get("sha")
    if p.get("from_refcache"):
        pass  # a cached CPU reference: its key covers the CPU-path sources, build and inputs (refcache.py)
    elif want and (not have or not have.startswith(want)):
        raise Fail(f"{rd}: built from {have}, acceptance names {want}")
    if p.get("mutated") and not acc.get("_allow_mutated"):
        raise Fail(f"{rd}: mutated copy ({p['mutated']}) used outside a control")
    if p.get("binary_sha256") != p.get("build", {}).get("binary_sha256"):
        raise Fail(f"{rd}: binary sha256 differs from its build record")
    if p.get("debug_env_present") and not (allow_debug and p.get("allow_debug")):
        raise Fail(f"{rd}: debug variables in a production run: {p['debug_env_present']}")
    excl = acc.get("gpu_exclusive", True) if gpu_exclusive is None else gpu_exclusive
    if excl and [l for l in p.get("gpu_before", []) if l.strip()]:
        raise Fail(f"{rd}: GPU not idle at start: {p['gpu_before']}")
    return p


# ------------------------------------------------------------------ parsers
def info_precision(rd):
    """Precision Green/tau (mean, max) and Phase max from info (fixed labels; a miss is a FAIL)."""
    t = read(os.path.join(rd, "info")).replace("\n", " ")
    out = {}
    for key, label, n in (("green", "Precision Green", 2), ("tau", "Precision tau", 2), ("phase", "Precision Phase", 1)):
        m = re.search(label + r"[^:]*:\s+" + r"\s+".join([f"({FLOAT})"] * n), t)
        if not m:
            if key == "tau" and "Precision tau" not in t:
                continue  # Ltau = 0 run: no tau line (the gate decides whether that is allowed)
            raise Fail(f"{rd}: cannot parse '{label}' in info")
        out[key] = [num(x) for x in m.groups()]
    return out


def h5_datasets(path, only_obser=False):
    import h5py
    import numpy as np
    out = {}
    with h5py.File(need(path), "r") as f:
        def visit(n, o):
            if isinstance(o, h5py.Dataset) and (not only_obser or n.endswith("obser")):
                out[n] = np.asarray(o[()])
        f.visititems(visit)
    if not out:
        raise Fail(f"{path}: no datasets")
    return out


def taudiag(rd):
    """Rows of the TAUDIAG lines (ALF_GPU_TAUDIAG=2): pass, nst, nfe, br, dev[4], hostfull[4], P[5]."""
    rows = []
    for line in read(os.path.join(rd, "run.log")).splitlines():
        if "TAUDIAG" not in line:
            continue
        t = line.split()
        try:
            i, j, k = t.index("dev"), t.index("hostfull"), t.index("P")
            rows.append({"pass": int(t[t.index("pass") + 1]), "nst": int(t[t.index("nst") + 1]),
                         "dev": [num(x) for x in t[i + 1:i + 5]], "host": [num(x) for x in t[j + 1:j + 5]],
                         "P": [num(x) for x in t[k + 1:k + 6]]})
        except (ValueError, IndexError) as e:
            raise Fail(f"{rd}: malformed TAUDIAG line: {line[:80]} ({e})")
    if not rows:
        raise Fail(f"{rd}: no TAUDIAG lines")
    return rows


# ------------------------------------------------------------------ gates
def g_selftest(g, root, acc):
    """Self-test runs: rc 0, exactly `passed` PASSED lines, no FAILED, all reported numbers finite."""
    res = []
    for r in g["runs"]:
        rd = os.path.join(root, r)
        check_prov(rd, acc, 0, allow_debug=True)
        log = read(os.path.join(rd, "run.log"))
        npass, nfail = len(re.findall(r"self-test PASSED", log)), len(re.findall(r"FAILED", log))
        if npass != g.get("passed", 1) or nfail:
            raise Fail(f"{r}: {npass} PASSED (expected {g.get('passed', 1)}), {nfail} FAILED")
        if re.search(r"\b(?:NaN|Infinity)\b", log, re.I):
            raise Fail(f"{r}: NaN/Inf in the log")
        for pat in g.get("require_lines", []):
            if not re.search(pat, log):
                raise Fail(f"{r}: required line missing: {pat}")
        mx = None
        for pat in g.get("max_patterns", []):
            vals = [num(v) for v in re.findall(pat, log)]
            if not vals:
                raise Fail(f"{r}: pattern {pat} matched nothing")
            mx = max([mx or 0.0] + vals)
        if mx is not None and mx > g["max_rel"]:
            raise Fail(f"{r}: max {mx:.3e} > {g['max_rel']:.1e}")
        res.append({"run": r, "passed": npass, "max": mx})
    return res


def g_faults(g, root, acc):
    """Fault controls: each run must exit 1 and print FAILED naming the expected sub-test (a crash is not a pass)."""
    res = []
    for r, expect in g["cases"].items():
        rd = os.path.join(root, r)
        check_prov(rd, acc, 1, allow_debug=True)
        log = read(os.path.join(rd, "run.log"))
        fl = [l for l in log.splitlines() if "FAILED" in l]
        if not fl or not any(expect in l for l in fl):
            raise Fail(f"{r}: no FAILED line naming '{expect}'")
        for pat in g.get("armed_lines", []):
            if not re.search(pat, log):
                raise Fail(f"{r}: fault not shown as armed ({pat})")
        res.append({"run": r, "failed_line": fl[0].strip()[:120]})
    return res


def _dev_per_bin(a, b, nbin_expect):
    import numpy as np
    if a.shape != b.shape:
        raise Fail(f"shape {a.shape} vs {b.shape}")
    if a.shape[0] != nbin_expect:
        raise Fail(f"{a.shape[0]} bins, expected {nbin_expect}")
    if not (np.isfinite(a).all() and np.isfinite(b).all()):
        raise Fail("non-finite values")
    nb = a.shape[0]
    scale = max(float(np.abs(a).max()), 1e-300)
    d = np.abs(a - b).reshape(nb, -1).max(axis=1) / scale
    floor = 1e-6 * scale  # per-entry relative metric with a floor, reported
    per = (np.abs(a - b) / np.maximum(np.abs(a), floor)).reshape(nb, -1).max(axis=1)
    return d, per


def g_same_seed(g, root, acc):
    """GPU vs CPU oracle, same seeds, one sweep per bin: per-bin deviation of every */obser dataset."""
    import numpy as np
    res = []
    for ra, rb in g["pairs"]:
        A, B = os.path.join(root, ra), os.path.join(root, rb)
        check_prov(A, acc, 0)
        check_prov(B, acc, 0, gpu_exclusive=False)
        for pat in g.get("active_lines", []):
            if pat not in read(os.path.join(A, "run.log")) + read(os.path.join(A, "info")):
                raise Fail(f"{ra}: GPU path not active ({pat})")
        if read(os.path.join(A, "seeds")) != read(os.path.join(B, "seeds")):
            raise Fail(f"{ra}/{rb}: seeds differ (same-seed gate)")
        da = h5_datasets(os.path.join(A, "data.h5"), only_obser=True)
        db = h5_datasets(os.path.join(B, "data.h5"), only_obser=True)
        if set(da) != set(db):
            raise Fail(f"dataset sets differ: {sorted(set(da) ^ set(db))}")
        if len(da) != g["n_datasets"]:
            raise Fail(f"{len(da)} datasets, expected {g['n_datasets']}")
        tau = [n for n in da if "_tau/" in n]
        if len(tau) != g["n_tau_datasets"]:
            raise Fail(f"{len(tau)} tau datasets, expected {g['n_tau_datasets']}")
        dev_all = dev_tau = per_all = None
        for n in sorted(da):
            d, per = _dev_per_bin(da[n], db[n], g["nbin"])
            dev_all = d if dev_all is None else np.maximum(dev_all, d)
            per_all = per if per_all is None else np.maximum(per_all, per)
            if n in tau:
                dev_tau = d if dev_tau is None else np.maximum(dev_tau, d)
        out = {"pair": [ra, rb], "overall_max": float(dev_all.max()), "per_entry_max": float(per_all.max()),
               "tau_max": float(dev_tau.max()) if dev_tau is not None else None,
               "first_bin_above_1e-9": int(np.argmax(dev_all > 1e-9)) + 1 if (dev_all > 1e-9).any() else None}
        if (dev_all > g["no_bin_above"]).any():
            raise Fail(f"{ra}: a bin above {g['no_bin_above']:.0e} (divergent trajectory)")
        if out["overall_max"] > g["overall_max"]:
            raise Fail(f"{ra}: overall max {out['overall_max']:.3e} > {g['overall_max']:.1e}")
        if dev_tau is not None and out["tau_max"] > g["tau_max"]:
            raise Fail(f"{ra}: tau max {out['tau_max']:.3e} > {g['tau_max']:.1e}")
        if "per_entry_max" in g and out["per_entry_max"] > g["per_entry_max"]:
            raise Fail(f"{ra}: per-entry max {out['per_entry_max']:.3e} > {g['per_entry_max']:.1e}")
        res.append(out)
    return res


def g_replay(g, root, acc):
    """Bit identity of two runs (same binary re-run, or two builds of one commit): every dataset and the
    full-precision Precision lines. With "datasets": "equal_time" only non-tau datasets are compared
    (non-interference gate)."""
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
            raise Fail(f"dataset sets differ: {sorted(set(da) ^ set(db))}")
        if len(da) != g["n_datasets"]:
            raise Fail(f"{len(da)} datasets, expected {g['n_datasets']}")
        diff = [n for n in sorted(da) if da[n].shape != db[n].shape or da[n].tobytes() != db[n].tobytes()]
        for n in da:
            if np.issubdtype(da[n].dtype, np.number) and not np.isfinite(da[n]).all():
                raise Fail(f"{ra}: non-finite values in {n}")
        if g.get("datasets") != "equal_time":
            pa = [l.strip() for l in read(os.path.join(A, "info")).splitlines() if "Precision" in l]
            pb = [l.strip() for l in read(os.path.join(B, "info")).splitlines() if "Precision" in l]
            if not pa or pa != pb:
                diff.append("info:Precision")
        if diff:
            raise Fail(f"{ra} vs {rb}: {len(diff)} differ, first {diff[0]}")
        res.append({"pair": [ra, rb], "identical": True, "n": len(da)})
    return res


def paired_stats(rows, frm, ref_row="b/a"):
    """GTT pair (index 1 of dev/host), passes >= frm: log10(dev/host) and the 1-ulp reference.
    P = [P0 dev host-propagated, P1 host+1ulp(a), P2 host+1ulp(b), P3 dev+1ulp(a), P4 dev+1ulp(b)]."""
    sel = [r for r in rows if r["pass"] >= frm]
    if not sel:
        raise Fail(f"no TAUDIAG rows at pass >= {frm}")
    for r in sel:
        for v in [r["dev"][1], r["host"][1], r["P"][1], r["P"][2]]:
            if not v > 0:
                raise Fail(f"non-positive value in a used row (pass {r['pass']}): {v}")
    lr = [math.log10(r["dev"][1] / r["host"][1]) for r in sel]
    if ref_row == "b/a":
        ref = [math.log10(r["P"][2] / r["P"][1]) for r in sel]
    else:
        raise Fail(f"unknown reference row {ref_row}")
    return sel, lr, ref


def g_paired(g, root, acc):
    """Paired same-trajectory test (ticket 207 acceptance 3, as decided by the owner; thresholds from the file)."""
    res = []
    for r in g["runs"]:
        rd = os.path.join(root, r)
        check_prov(rd, acc, 0, allow_debug=True)
        prov = json.loads(read(os.path.join(rd, "provenance.json")))
        for kv in g.get("env", []):
            if kv not in prov.get("extra_env", []):
                raise Fail(f"{r}: required setting {kv} not used")
        rows = taudiag(rd)
        if len(rows) != g["n_points"]:
            raise Fail(f"{r}: {len(rows)} TAUDIAG lines, expected {g['n_points']} (incomplete run?)")
        sel, lr, ref = paired_stats(rows, g["from_pass"], g.get("ref_row", "b/a"))
        if len(sel) != g["n_used"]:
            raise Fail(f"{r}: {len(sel)} rows used, expected {g['n_used']}")
        sd, sref, mean = statistics.pstdev(lr), statistics.pstdev(ref), statistics.mean(lr)
        cnt = {}
        for t in g["count_thresholds"]:
            nd = sum(r_["dev"][1] > t for r_ in rows)  # all points, GTT pair
            nh = sum(r_["host"][1] > t for r_ in rows)
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
            raise Fail(f"{r}: " + "; ".join(bad))
    return res


def g_precision(g, root, acc):
    """Precision mean ratio GPU / larger CPU mean <= limit (Green and tau); max reported only."""
    res = []
    for case in g["cases"]:
        gp = info_precision(os.path.join(root, case["gpu"]))
        check_prov(os.path.join(root, case["gpu"]), acc, 0, allow_debug=case.get("allow_debug", False))
        cps = []
        for c in case["cpu"]:
            check_prov(os.path.join(root, c), acc, 0, gpu_exclusive=False)
            cps.append(info_precision(os.path.join(root, c)))
        pg = json.loads(read(os.path.join(root, case["gpu"], "provenance.json")))["inputs"]
        for c in case["cpu"]:
            pc = json.loads(read(os.path.join(root, c, "provenance.json")))["inputs"]
            if pc != pg:
                raise Fail(f"{case['gpu']} vs {c}: parameters or seeds differ")
        out = {"gpu": case["gpu"]}
        for key in g.get("quantities", ["green", "tau"]):
            if key not in gp or any(key not in c for c in cps):
                raise Fail(f"{case['gpu']}: Precision {key} missing")
            ratio = gp[key][0] / max(c[key][0] for c in cps)
            out[key] = {"gpu_mean": gp[key][0], "gpu_max": gp[key][1], "cpu_mean_max": max(c[key][0] for c in cps),
                        "mean_ratio": ratio, "max_ratio_info": gp[key][1] / max(c[key][1] for c in cps)}
            if ratio > g["mean_ratio_max"]:
                raise Fail(f"{case['gpu']}: {key} mean ratio {ratio:.3f} > {g['mean_ratio_max']}")
        res.append(out)
    return res


def g_compare2(g, root, acc):
    """compare2.py (fail-closed version) on two independent-seed runs; PASS required."""
    A, B = os.path.join(root, g["a"]), os.path.join(root, g["b"])
    check_prov(A, acc, 0)
    check_prov(B, acc, 0, gpu_exclusive=False)
    here = os.path.dirname(os.path.abspath(__file__))
    js = os.path.join(tempfile.mkdtemp(), "c2.json")
    cmd = [sys.executable, os.path.join(here, "..", "compare2.py"), A, B, "--json", js,
           "--expect", ";".join(g.get("expect", [])), "--expect-m", str(g.get("m", 1))]
    p = subprocess.run(cmd, capture_output=True, text=True)
    out = {"exit": p.returncode, "tail": p.stdout.strip().splitlines()[-1:] if p.stdout else []}
    if p.returncode != 0:
        raise Fail(f"compare2 exit {p.returncode}: {out['tail']} {p.stderr.strip()[-200:]}")
    c = json.loads(read(js))
    out.update({"m": c["m"], "sum_t2": c["sum_t2"], "max_abs_t": max(abs(r["t"]) for r in c["rows"])})
    return [out]


def g_log_regex(g, root, acc):
    """Generic check on a log file: every regex must match; captured integers must equal the expected values;
    'exit_file' (if given) must contain 0. Used for ctest and ppcmp logs."""
    res = []
    for chk in g["checks"]:
        t = read(os.path.join(root, chk["file"]))
        m = re.search(chk["regex"], t, re.M)
        if not m:
            raise Fail(f"{chk['file']}: no match for {chk['regex']}")
        if "equals" in chk and [int(x) for x in m.groups()] != chk["equals"]:
            raise Fail(f"{chk['file']}: {m.groups()} != {chk['equals']}")
        if "exit_file" in chk and read(os.path.join(root, chk["exit_file"])).strip() != "0":
            raise Fail(f"{chk['exit_file']}: nonzero exit")
        res.append({"file": chk["file"], "groups": list(m.groups())})
    return res


GATES = {"selftest": g_selftest, "faults": g_faults, "same_seed": g_same_seed, "replay": g_replay,
         "paired": g_paired, "precision": g_precision, "compare2": g_compare2, "log_regex": g_log_regex}


def run_gate(g, root, acc):
    try:
        return "PASS", GATES[g["type"]](g, root, acc), None
    except Fail as e:
        return "FAIL", None, str(e)
    except Exception as e:  # anything unexpected is a FAIL, never a pass
        return "FAIL", None, f"{type(e).__name__}: {e}"


def evaluate(acc, root):
    out = {"schema": "alf-acc-verdict/1", "ticket": acc.get("ticket"), "candidate_sha": acc.get("candidate_sha"),
           "acceptance_sha256": acc.get("_sha256"), "gates": [], "controls": []}
    ok = True
    for g in acc["gates"]:
        v, r, why = run_gate(g, root, acc)
        out["gates"].append({"id": g["id"], "type": g["type"], "verdict": v, "result": r, "reason": why})
        ok &= v == "PASS"
    for c in acc.get("controls", []):
        v, r, why = run_gate(c, root, dict(acc, _allow_mutated=True))
        killed = v == "FAIL"
        exp = c.get("expect", "FAIL")
        good = (exp == "FAIL" and killed) or (exp == "PASS" and not killed)
        out["controls"].append({"id": c["id"], "gate": c["type"], "expect": exp, "outcome": v, "as_expected": good,
                                "reason": why})
        ok &= good
    out["overall"] = "PASS" if ok else "FAIL"
    return out


def write_md(v, path):
    L = [f"# Verdict ticket {v['ticket']} candidate {v['candidate_sha']}", "",
         f"acceptance file sha256 `{v['acceptance_sha256']}`", "", f"**Overall: {v['overall']}**", "",
         "| gate | type | verdict | detail |", "|---|---|---|---|"]
    for g in v["gates"]:
        L.append(f"| {g['id']} | {g['type']} | {g['verdict']} | {(g['reason'] or json.dumps(g['result'], default=str))[:300]} |")
    L += ["", "| control | gate | expected | outcome | as expected | reason |", "|---|---|---|---|---|---|"]
    for c in v["controls"]:
        L.append(f"| {c['id']} | {c['gate']} | {c['expect']} | {c['outcome']} | {c['as_expected']} | {(c['reason'] or '')[:200]} |")
    open(path, "w").write("\n".join(L) + "\n")


# ------------------------------------------------------------------ self-check (positive controls of this script)
def selfcheck():
    """Synthetic run directories through every gate type: a clean case must PASS, each mutation must FAIL."""
    import h5py
    import numpy as np
    tmp = tempfile.mkdtemp(prefix="verdict_selfcheck_")
    rng = np.random.default_rng(12345)
    acc = {"candidate_sha": "c0ffee", "gpu_exclusive": True}
    build = {"binary_sha256": "b" * 64, "candidate": {"sha": "c0ffee" + "0" * 34}, "source": "src"}

    def mkrun(name, data, info, log="", extra_env=(), debug=(), rc=0, seeds="1 2 3\n"):
        d = os.path.join(tmp, name)
        os.makedirs(d)
        with h5py.File(os.path.join(d, "data.h5"), "w") as f:
            for k, v in data.items():
                f[k] = v
        open(os.path.join(d, "info"), "w").write(info)
        open(os.path.join(d, "run.log"), "w").write(log)
        open(os.path.join(d, "seeds"), "w").write(seeds)
        open(os.path.join(d, "parameters"), "w").write("&VAR_x\n/\n")
        prov = {"schema": "alf-acc-prov/1", "binary_sha256": "b" * 64, "build": build, "exit_code": rc,
                "timed_out": False, "debug_env_present": list(debug), "allow_debug": bool(debug),
                "extra_env": list(extra_env), "gpu_before": [],
                "inputs": {"parameters_sha256": "p", "seeds_sha256": "s"}}
        json.dump(prov, open(os.path.join(d, "provenance.json"), "w"))
        return name

    def info(gm, gx, tm, tx):
        return (f"  Precision Green  Mean, Max :    {gm:.16E} \n   {gx:.16E}\n  Precision Phase, Max       :     0.0     \n"
                f"  Precision tau    Mean, Max :    {tm:.16E} \n   {tx:.16E}\n")

    base = {f"Ener_scal/obser": rng.normal(size=(20, 1)), "Green_tau/obser": rng.normal(size=(20, 4, 3)) * 1e-3,
            "Green_eq/obser": rng.normal(size=(20, 3))}
    ok_info = info(5e-11, 1e-6, 6e-12, 5e-7)
    act = "GPU tau pass (207): ACTIVE"
    mkrun("A", base, ok_info, log=act)
    mkrun("B", {k: v * (1 + 1e-12) for k, v in base.items()}, ok_info)
    mut = dict(base)
    g = mut["Green_tau/obser"].copy()
    g[7, 3, 2] *= 1 + 1e-6  # a small entry, relative 1e-6 (shows as < 1e-9 overall: per-entry metric catches it)
    mut["Green_tau/obser"] = g
    mkrun("B_small", mut, ok_info)
    nanb = {k: v.copy() for k, v in base.items()}
    nanb["Ener_scal/obser"][3, 0] = float("nan")
    mkrun("B_nan", nanb, ok_info)
    mkrun("B_missing", {k: v for k, v in base.items() if k != "Green_tau/obser"}, ok_info)
    mkrun("B_seed", base, ok_info, seeds="9 9 9\n")
    mkrun("A_inactive", base, ok_info, log="GPU tau pass (207): INACTIVE")
    mkrun("A_debug", base, ok_info, log=act, debug=["ALF_GPU_NS"])
    mkrun("A_rc", base, ok_info, log=act, rc=137)
    mkrun("A_twin", base, ok_info, log=act)
    nan_info = info(float("nan"), 1e-6, 6e-12, 5e-7)
    assert nan_info != ok_info and "NAN" in nan_info.upper()  # the mutation must really change the input
    mkrun("A_nanprec", base, nan_info, log=act)

    ss = {"type": "same_seed", "id": "ss", "n_datasets": 3, "n_tau_datasets": 1, "nbin": 20, "no_bin_above": 1e-4,
          "overall_max": 3e-9, "tau_max": 1e-9, "per_entry_max": 1e-8, "active_lines": [act]}
    rp = {"type": "replay", "id": "rp", "n_datasets": 3}
    expect = [("same_seed clean", dict(ss, pairs=[["A", "B"]]), "PASS"),
              ("same_seed small entry 1e-6", dict(ss, pairs=[["A", "B_small"]]), "FAIL"),
              ("same_seed NaN", dict(ss, pairs=[["A", "B_nan"]]), "FAIL"),
              ("same_seed missing dataset", dict(ss, pairs=[["A", "B_missing"]]), "FAIL"),
              ("same_seed other seeds", dict(ss, pairs=[["A", "B_seed"]]), "FAIL"),
              ("same_seed GPU inactive", dict(ss, pairs=[["A_inactive", "B"]]), "FAIL"),
              ("same_seed debug env", dict(ss, pairs=[["A_debug", "B"]]), "FAIL"),
              ("same_seed crash rc 137", dict(ss, pairs=[["A_rc", "B"]]), "FAIL"),
              ("replay identical", dict(rp, pairs=[["A", "A_twin"]]), "PASS"),
              ("replay 1e-12", dict(rp, pairs=[["A", "B"]]), "FAIL"),
              ("precision NaN", {"type": "precision", "id": "pr", "mean_ratio_max": 1.5,
                                 "cases": [{"gpu": "A_nanprec", "cpu": ["B"]}]}, "FAIL"),
              ("precision clean", {"type": "precision", "id": "pr", "mean_ratio_max": 1.5,
                                   "cases": [{"gpu": "A", "cpu": ["B"]}]}, "PASS"),
              ("precision 2x mean", {"type": "precision", "id": "pr", "mean_ratio_max": 1.5,
                                     "cases": [{"gpu": mkrun("A_2x", base, info(1e-10, 1e-6, 6e-12, 5e-7), log=act),
                                                "cpu": ["B"]}]}, "FAIL")]

    # paired: synthetic TAUDIAG with known statistics
    def taulog(n_pass, sd_dev, bias, nan_at=None):
        lines = []
        for ps in range(1, n_pass + 1):
            for nst in range(1, 5):
                for nfe in (1, 2):
                    h = 10 ** rng.normal(-8, 0.5)
                    pa, pb = h * 10 ** rng.normal(0, 0.1), h * 10 ** rng.normal(0, 0.1)
                    dv = h * 10 ** (rng.normal(0, sd_dev) + bias)
                    if nan_at == (ps, nst, nfe):
                        dv = float("nan")
                    lines.append(f" TAUDIAG pt 0 pass {ps} nst {nst} nfe {nfe} br 1 dev 1E-12 {dv:.3E} 1E-12 1E-12 | "
                                 f"hostfull 1E-12 {h:.3E} 1E-12 1E-12 | P {h:.3E} {pa:.3E} {pb:.3E} {pa:.3E} {pb:.3E} |")
        return "\n".join(lines)
    pg = {"type": "paired", "id": "pd", "n_points": 400, "n_used": 400, "from_pass": 1, "abs_mean_max": 0.05,
          "sd_ratio_max": 1.25, "count_thresholds": [1e-7], "count_k": 2.0}
    expect += [("paired clean", dict(pg, runs=[mkrun("P_ok", base, ok_info, taulog(50, 0.1, 0.0), debug=["ALF_GPU_TAUDIAG"])]), "PASS"),
               ("paired variance x1.5", dict(pg, runs=[mkrun("P_var", base, ok_info, taulog(50, 0.15 * 1.5 / 1.0, 0.0), debug=["ALF_GPU_TAUDIAG"])]), "FAIL"),
               ("paired bias x1.3", dict(pg, runs=[mkrun("P_bias", base, ok_info, taulog(50, 0.1, math.log10(1.3)), debug=["ALF_GPU_TAUDIAG"])]), "FAIL"),
               ("paired NaN", dict(pg, runs=[mkrun("P_nan", base, ok_info, taulog(50, 0.1, 0.0, nan_at=(3, 2, 1)), debug=["ALF_GPU_TAUDIAG"])]), "FAIL"),
               ("paired truncated", dict(pg, runs=[mkrun("P_short", base, ok_info, taulog(40, 0.1, 0.0), debug=["ALF_GPU_TAUDIAG"])]), "FAIL")]
    # self-test / faults
    st = {"type": "selftest", "id": "st", "passed": 1, "max_rel": 1e-10, "max_patterns": [r"G00\s+(" + FLOAT + r")"]}
    expect += [("selftest clean", dict(st, runs=[mkrun("S_ok", base, ok_info, "G00  3.5E-14\n GPU tau self-test PASSED", debug=["ALF_GPU_SELFTEST"])]), "PASS"),
               ("selftest crash", dict(st, runs=[mkrun("S_crash", base, ok_info, "G00  3.5E-14\n", rc=139, debug=["ALF_GPU_SELFTEST"])]), "FAIL"),
               ("selftest above tol", dict(st, runs=[mkrun("S_big", base, ok_info, "G00  3.5E-09\n GPU tau self-test PASSED", debug=["ALF_GPU_SELFTEST"])]), "FAIL"),
               ("fault detected", {"type": "faults", "id": "ft", "cases": {mkrun("F_ok", base, ok_info, " GPU self-test FAILED: tau: pass", rc=1, debug=["ALF_GPU_FAULT_TAU"]): "tau: pass"}}, "PASS"),
               ("fault = crash", {"type": "faults", "id": "ft", "cases": {mkrun("F_crash", base, ok_info, "Segmentation fault", rc=139, debug=["ALF_GPU_FAULT_TAU"]): "tau: pass"}}, "FAIL"),
               ("fault wrong test", {"type": "faults", "id": "ft", "cases": {mkrun("F_other", base, ok_info, " GPU self-test FAILED: hopping", rc=1, debug=["ALF_GPU_FAULT_TAU"]): "tau: pass"}}, "FAIL")]
    bad = 0
    for name, g, exp in expect:
        v, r, why = run_gate(g, tmp, acc)
        flag = "ok " if v == exp else "BAD"
        bad += v != exp
        print(f"{flag} {name:32s} expected {exp:4s} got {v:4s} {('(' + why + ')') if why else ''}"[:200])
    shutil.rmtree(tmp)
    print("selfcheck:", "PASS" if not bad else f"FAIL ({bad} wrong)")
    return 0 if not bad else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("acceptance", nargs="?")
    ap.add_argument("runroot", nargs="?")
    ap.add_argument("--out", default=None)
    ap.add_argument("--candidate", default=None, help="candidate sha (replaces \"self\" in the acceptance file)")
    ap.add_argument("--selfcheck", action="store_true")
    a = ap.parse_args()
    if a.selfcheck:
        return selfcheck()
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
    v = evaluate(acc, a.runroot)
    out = a.out or a.runroot
    json.dump(v, open(os.path.join(out, "verdict.json"), "w"), indent=1, default=str)
    write_md(v, os.path.join(out, "verdict.md"))
    print(open(os.path.join(out, "verdict.md")).read())
    return 0 if v["overall"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
