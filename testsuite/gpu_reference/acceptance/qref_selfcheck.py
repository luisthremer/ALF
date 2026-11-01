#!/usr/bin/env python3
"""Self-check of the QREF gates of verdict.py (ticket 207-5), called by `verdict.py --selfcheck`.

Synthetic QREF runs (qref_synth.py) go through the real 207-2 scripts (parser, aggregation, decision) and the gates of qref_gate.py: the
clean input passes, and every planted defect gives its outcome (and, where named, its reason).  CPU only.
"""
import copy
import hashlib
import json
import os
import re
import shutil
import time

import qref_gate
import qref_synth as S


def _edit_json(q, name, **kw):
    p = os.path.join(q, name)
    d = json.load(open(p))
    d.update(kw)
    json.dump(d, open(p, "w"), indent=1)


def selfcheck(V):
    """returns the number of wrong cases"""
    import h5py
    import numpy as np
    t0 = time.time()
    tmp = V.tmpdir()
    root, bindir = os.path.join(tmp, "runs"), os.path.join(tmp, "bin")
    os.makedirs(root)
    qbs = {"def": S.make_qbuild(bindir, "def")}
    acc = {"ticket": "zz9", "candidate_sha": S.CAND[:12], "gpu_exclusive": True, "active_lines": ["GPU tau pass (207): ACTIVE"],
           "builds": {"gpu": {"tool": "nvhpc", "configure": "PGI noMPI HDF5 GPU=nvidia", "gpu": True}}, "configs": {"synth": {"parameters_sha256": "x"}},
           "qbuilds": {"def": {"configure": "GNU noMPI", "qflags": "-DQR_PREPIVOT", "flags_ext": "-DQR_PREPIVOT"}}, "runs": {}}
    nul = ",".join(str(c) for c in S.NUL["def"])
    bad = 0

    def mk(rid, **kw):
        acc["runs"][rid] = S.make_run(root, rid, **kw)
        S.seal(root, rid, qbs[kw.get("path", "def")])
        return rid

    def clone(new, src, fn=None):
        shutil.copytree(os.path.join(root, src), os.path.join(root, new))
        acc["runs"][new] = copy.deepcopy(acc["runs"][src])
        if fn:
            fn(os.path.join(root, new, "qproc"))
        return new

    def edit_q(fn):
        def f(q):
            p = os.path.join(q, "qref.txt")
            ls = open(p).read().split("\n")[:-1]
            open(p, "w").write("\n".join(fn(ls[:-1]) + [ls[-1]]) + "\n")
        return f

    def swap_blocks(a, b):
        def one(line):
            segs = line.split(" | ")
            ia = next(i for i, s in enumerate(segs) if s.split()[0] == a)
            ib = next(i for i, s in enumerate(segs) if s.split()[0] == b)
            na, nb = segs[ia].split(None, 1)[1], segs[ib].split(None, 1)[1]
            segs[ia], segs[ib] = "%s %s" % (a, nb), "%s %s" % (b, na)
            return " | ".join(segs)
        return lambda ls: [one(x) for x in ls]

    mk("m1", seed=1, scale={"var8": 10.0})
    mk("m2", seed=2, scale={"var8": 10.0})
    mk("c1", seed=3, role="ctl", scale={"dev": 10.0})
    G0 = {"id": "qref_def", "type": "qref", "path": "def", "stride": 3, "qs": [20], "keep": 39, "nul": nul, "runs": ["m1", "m2"], "ctl_runs": ["c1"],
          "refs": "def", "bound": 0.15, "expect_points": 20, "required_controls": ["var8", "ctl_dev"],
          "pointwise": {"measure": "emax", "bound": 0.15, "expect_points": 40}, "tools_sha256": S.tool_hashes()}

    def case(name, expect, gate=None, why=None, a=None):
        nonlocal bad
        v, r, w = V.run_gate(gate or G0, root, a or acc)
        ok = v == expect and (why is None or (w is not None and why in w))
        bad += not ok
        extra = json.dumps({k: r[k]["verdict"] for k in ("pTT", "pointwise")}) if (v == "PASS" and r and "pTT" in r) else ""
        print(f"{'ok ' if ok else 'BAD'} qref {name:48s} expected {expect}{' (' + why + ')' if why else ''} got {v}"
              f"{'' if ok else ': ' + str(w)[:300]} {extra}", flush=True)
        return v, r, w

    v, r, w = case("clean input", "PASS")
    if r:
        n = r["pTT"]["numbers"]
        print(f"    clean: S_ens(pTT) dev {n['dev']['S_ens']:+.3f} (n {n['dev']['n']}), var8 {n['var8']['S_ens']:+.3f}, ctl_dev {n['ctl_dev']['S_ens']:+.3f}; "
              f"pointwise dev {r['pointwise']['numbers']['dev']['S_ens']:+.3f}, ctl_dev {r['pointwise']['numbers']['ctl_dev']['S_ens']:+.3f}")
        dm = n["dev"]["S_ens"]
        case("bound just above the measured device value", "PASS", dict(G0, bound=dm + 0.01))
        case("bound just below the measured device value", "FAIL", dict(G0, bound=dm - 0.01), why="device S_ens")
    # the harness control: the judged device replaced by the stand-in defect (var8) must FAIL, through evaluate() with lint
    acc_c = dict(acc, gates=[G0], controls=[{"id": "c_device_is_var8", "like": "qref_def", "replace": {"device": "var8"}, "expect": "FAIL"}], _status="")
    acc_c["runs"] = {k: acc["runs"][k] for k in ("m1", "m2", "c1")}
    ev = V.evaluate(acc_c, root)
    ok = ev["gates"][0]["verdict"] == "PASS" and [c["outcome"] for c in ev["controls"]] == ["FAIL"] and ev["controls"][0]["as_expected"] and not ev["lint"]
    bad += not ok
    print(f"{'ok ' if ok else 'BAD'} qref {'control device := var8 is killed (evaluate)':48s} gate {ev['gates'][0]['verdict']}, control "
          f"{[(c['outcome'], c['reason']) for c in ev['controls']]} lint {ev['lint']}"[:340])

    def with_m1(name, fn, expect, why=None):
        clone(name, "m1", edit_q(fn))
        case(name, expect, dict(G0, runs=[name, "m2"]), why)

    def tok_edit(i, repl):
        def f(ls):
            t = ls[i].split(" | ")
            t[1] = t[1].replace(t[1].split()[3], repl, 1)       # third number of the dev member
            ls[i] = " | ".join(t)
            return ls
        return f
    with_m1("m1_nan", tok_edit(5, "NaN"), "INPUT_ERROR", "ES12.4E3")
    with_m1("m1_infinity", tok_edit(5, "Infinity"), "INPUT_ERROR", "ES12.4E3")
    with_m1("m1_noE", tok_edit(5, "1.2000+300"), "INPUT_ERROR", "ES12.4E3")
    with_m1("m1_missing_member", lambda ls: ls[:6] + [" | ".join(s for s in ls[6].split(" | ") if not s.startswith("n5 "))] + ls[7:],
            "INPUT_ERROR", "members, expected")
    with_m1("m1_missing_point", lambda ls: ls[:8] + ls[9:], "INPUT_ERROR")

    def shift_nob(ls):
        i = next(j for j, x in enumerate(ls) if " hasb 0 " in x)
        segs = ls[i].split(" | ")
        t = segs[1].split()
        segs[1] = " ".join(t[:1] + t[-1:] + t[1:-1])      # the 8 numbers of dev rotated by one place: the -1 is no longer in its slot
        ls[i] = " | ".join(segs)
        return ls
    def pos_in_nob(ls):
        i = next(j for j, x in enumerate(ls) if " hasb 0 " in x)
        ls[i] = ls[i].replace("-1.0000E+000", "1.0000E-009", 1)       # a propagated number where the point has no B
        return ls
    with_m1("m1_shifted_fields_no_B", shift_nob, "INPUT_ERROR", "negative error number")
    with_m1("m1_propagated_number_without_B", pos_in_nob, "INPUT_ERROR", "propagated numbers")

    def hdr(ls):
        ls[3] = ls[3].replace(" nfe 1 ", " nfe 2 ", 1) if " nfe 1 " in ls[3] else ls[3].replace(" nfe 2 ", " nfe 1 ", 1)
        return ls
    with_m1("m1_header_nfe", hdr, "INPUT_ERROR", "header nfe")
    clone("m1_swap", "m1", edit_q(swap_blocks("dev", "var8")))
    clone("m2_swap", "m2", edit_q(swap_blocks("dev", "var8")))
    case("device and var8 swapped (both main runs)", "FAIL", dict(G0, runs=["m1_swap", "m2_swap"]), why="device S_ens")
    mk("c1_weak", seed=3, role="ctl")           # a production control that does not show its defect
    case("control run not killed (ctl_dev ~ 0)", "FAIL", dict(G0, ctl_runs=["c1_weak"]), why="NOT KILLED")
    mk("m1_pw", seed=1, scale={"var8": 10.0, "dev": (10.0, 1.0)})      # errors x10, propagated numbers unchanged
    mk("m1_pt", seed=1, scale={"var8": 10.0, "dev": (1.0, 10.0)})      # propagated numbers x10, errors unchanged
    v, r, w = case("pointwise-only defect (emax x10)", "FAIL", dict(G0, runs=["m1_pw", "m2"]), why="pointwise")
    bad += bool(w and "S_ens(pTT)" in w)
    v, r, w = case("pTT-only defect (propagated x10)", "FAIL", dict(G0, runs=["m1_pt", "m2"]), why="S_ens(pTT)")
    bad += bool(w and "pointwise" in w)
    case("pointwise bound null (DRAFT): not judged", "PASS", dict(G0, pointwise=dict(G0["pointwise"], bound=None)))
    case("pointwise control not killed", "FAIL", dict(G0, pointwise=dict(G0["pointwise"], ctl_runs=["c1_weak"])), why="NOT KILLED")
    mk("m1_env", seed=1, scale={"var8": 10.0}, extra_env=["ALF_GPU_TAUSWAP=2"])
    case("main run with an extra ALF_GPU_TAUSWAP", "INPUT_ERROR", dict(G0, runs=["m1_env", "m2"]), why="allow-list")
    case("expect_points off by one", "INPUT_ERROR", dict(G0, expect_points=19))
    case("pointwise expect_points off by one", "INPUT_ERROR", dict(G0, pointwise=dict(G0["pointwise"], expect_points=39)))
    case("wrong reference list", "INPUT_ERROR", dict(G0, refs="cpu,nvh,n4"))
    case("required control member absent", "INPUT_ERROR", dict(G0, required_controls=["var8", "ctl_dev", "ctl_graded"]), why="required control")
    case("unfrozen script (tools_sha256)", "INPUT_ERROR", dict(G0, tools_sha256=dict(G0["tools_sha256"], **{"t207_d5_gate.py": "0" * 64})),
         why="pre-registered script")
    case("same run twice", "INPUT_ERROR", dict(G0, runs=["m1", "m1"]))
    clone("m1_copy", "m1")
    case("a copy of a run (same qref.txt)", "INPUT_ERROR", dict(G0, runs=["m1", "m1_copy"]), why="same run")
    clone("m1_prog", "m1", lambda q: _edit_json(q, "qref_proc.json", prog="1c" * 32))
    case("program hash of the processing differs", "INPUT_ERROR", dict(G0, runs=["m1_prog", "m2"]), why="build record")
    clone("m1_info", "m1", lambda q: open(os.path.join(q, "info"), "a").write(" extra line\n"))
    case("qproc records differ from the harness provenance", "INPUT_ERROR", dict(G0, runs=["m1_info", "m2"]), why="harness provenance")
    clone("m1_noq", "m1", lambda q: shutil.rmtree(q))
    case("no qproc/ directory", "INPUT_ERROR", dict(G0, runs=["m1_noq", "m2"]), why="no qproc")
    clone("m1_dump", "m1", lambda q: os.remove(os.path.join(q, "qref", "q000040.bin")))
    case("a kept dump deleted after the processing", "INPUT_ERROR", dict(G0, runs=["m1_dump", "m2"]), why="check --dumps")
    clone("m1_qb", "m1", lambda q: _edit_json(q, "acc_qproc.json", qbuild_sha256="2d" * 32))
    case("build record not the one used", "INPUT_ERROR", dict(G0, runs=["m1_qb", "m2"]), why="build record")
    clone("m1_rc", "m1", lambda q: _edit_json(q, "acc_qproc.json", proc_rc=3))
    case("processing exited with an error", "INPUT_ERROR", dict(G0, runs=["m1_rc", "m2"]), why="processing status")
    case("offline program of another candidate", "INPUT_ERROR", G0, why="acceptance names", a=dict(acc, candidate_sha="deadbeef0000"))
    # lint
    le = qref_gate.lint(V, dict(acc, gates=[G0], controls=[]))
    ok = not [e for e in le if "no qref gate uses" not in e]
    bad += not ok
    print(f"{'ok ' if ok else 'BAD'} qref lint: clean file {le[:2]}")
    lcases = [("stride 2", dict(G0, stride=2), "stride must be 3"), ("no control run", dict(G0, ctl_runs=[]), "no production control"),
              ("pointwise bound null in a frozen file", dict(G0, pointwise=dict(G0["pointwise"], bound=None)), "DRAFT"),
              ("main run listed as control", dict(G0, ctl_runs=["m2"]), "declared in the ctl"),
              ("tools list incomplete", dict(G0, tools_sha256={}), "tools_sha256"), ("bound NaN", dict(G0, bound=float("nan")), "finite")]
    for nm, gg, sub in lcases:
        e = qref_gate.lint(V, dict(acc, gates=[gg], controls=[]))
        ok = any(sub in x for x in e)
        bad += not ok
        print(f"{'ok ' if ok else 'BAD'} qref lint: {nm:38s} expects '{sub}' {'found' if ok else e[:3]}")

    # h5_identical
    def mkh5(rid, env, data):
        d = os.path.join(root, rid)
        os.makedirs(d)
        with h5py.File(os.path.join(d, "data.h5"), "w") as f:
            for k, x in data.items():
                f[k] = x
        for nm, tx in (("info", "GPU tau pass (207): ACTIVE\n"), ("run.log", "x\n")):
            open(os.path.join(d, nm), "w").write(tx)

        def hh(f):
            return hashlib.sha256(open(os.path.join(d, f), "rb").read()).hexdigest()
        b = {"name": "gpu", "configure": "PGI noMPI HDF5 GPU=nvidia", "ALF_FLAGS_EXT": "", "binary_sha256": "b" * 64, "source": "src",
             "tar_sha256_built": S.TAR, "candidate": {"sha": S.CAND, "tar_sha256": S.TAR}}
        json.dump({"schema": "alf-acc-prov/1", "binary_sha256": "b" * 64, "build": b, "exit_code": 0, "timed_out": False, "extra_env": list(env),
                   "debug_env_present": [t.split("=")[0] for t in env], "env": {"OMP_TARGET_OFFLOAD": "MANDATORY"}, "gpu_check_rc": 0, "load_ok": True,
                   "gpu_before": [], "inputs": {"parameters_sha256": "p", "seeds_sha256": "s"},
                   "outputs": {f: hh(f) for f in ("info", "data.h5", "run.log")}}, open(os.path.join(d, "provenance.json"), "w"))
        acc["runs"][rid] = {"build": "gpu", "config": "synth", "env": list(env), "allow_debug": True, "active_lines": ["GPU tau pass (207): ACTIVE"]}
    rng = np.random.default_rng(5)
    data = {"Ener_scal/obser": rng.normal(size=(8, 2)), "Green_eq/obser": rng.normal(size=(8, 3))}
    mkh5("h_on", S.env_tokens("def", "main"), data)
    mkh5("h_off", [], data)
    GI = {"id": "qref_identity", "type": "h5_identical", "pairs": [["h_on", "h_off"]]}
    case("identity: data.h5 QREF on == off", "PASS", GI)
    d2 = {k: x.copy() for k, x in data.items()}
    d2["Green_eq/obser"][3, 1] = np.nextafter(d2["Green_eq/obser"][3, 1], 1e9)       # one ulp
    mkh5("h_off_ulp", [], d2)
    case("identity: one entry one ulp different", "FAIL", dict(GI, pairs=[["h_on", "h_off_ulp"]]), why="differ")
    mkh5("h_off_dbg", ["ALF_GPU_DEBUG=1"], data)
    case("identity: the 'off' run has a debug variable", "INPUT_ERROR", dict(GI, pairs=[["h_on", "h_off_dbg"]]), why="'off' run")
    case("identity: the 'on' run has no QREF", "INPUT_ERROR", dict(GI, pairs=[["h_off", "h_off"]]), why="'on' run")
    shutil.rmtree(tmp)
    print(f"qref self-check: {bad} wrong ({time.time() - t0:.0f} s)")
    return bad
