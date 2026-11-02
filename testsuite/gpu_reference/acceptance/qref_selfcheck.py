#!/usr/bin/env python3
"""Self-check of the QREF gates of verdict.py (ticket 207-5, second round), called by `verdict.py --selfcheck [--qref-fixture CANDROOT]`.

Part A (always, synthetic, seconds): the lint of qref gates (qref_gate.lint, which also runs the 207-2 gate's own spec check) on a synthetic
acceptance file, and the h5_identical gate on synthetic run directories.
Part B (needs --qref-fixture CANDROOT: a candidate root of a finished QREF dry run of acc_driver.py on THIS machine: runs/, bin/qref_*.qbuild.json,
qw/bin/ with the offline program, configs/seeds_record.json, src/ with the acceptance file): every case copies the fixture's QREF run
directories into its own directory, plants one defect, writes its own stamps (as the driver does; "not restamped" = the stamps of the clean
copy, i.e. the defect planted after the stamp), and runs verdict.run_gate -> qref_gate.g_qref -> the real t207_d5_gate.py (with the real
re-processing of the kept dumps).  Expected outcomes are fixed per case; the decision cases use spec variants whose bounds are derived from
the clean report (the self-check tests the mapping of the gate's verdicts, it is no calibration).
Returns (number of wrong cases, part B run).
"""
import copy
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time

import qref_gate as QG

CAND = "c0ffee000000" + "0" * 28
TAR = "7a" * 32
PAR = "5e" * 32  # a concrete parameters_sha256 (required outside a DRAFT, fix round 3)


def _res(V, name, expect, v, w, why, out):
    """why: a substring of the reason, or a compiled regular expression it must match"""
    hit = why is None or (w is not None and (why.search(w) is not None if hasattr(why, "search") else why in w))
    ok = v == expect and hit
    why = getattr(why, "pattern", why)
    print(f"{'ok ' if ok else 'BAD'} qref {name:58s} expected {expect}{' (' + why + ')' if why else ''} got {v}"
          + (f" | {(w or '')[:160]}" if ok else f": {(w or '')[:300]}"), flush=True)
    out.append((name, expect, v, ok))
    return ok


def _mem(nul, var=True):
    return ["dev"] + (["var8"] if var else []) + ["nvh", "cpu"] + ["n%d" % c for c in nul] + ["u1", "u2"]


NUL_DEF = [0, 4, 100004, 1, 100001, 2, 1608, 4808, 5, 100005]
# the environment names of a GPU=intel run on `intel` (read-only, 2026-11-01: env -i HOME USER LANG LC_ALL PATH TMPDIR bash --noprofile --norc -c
# '. ~/alf-gpu/env/intel.sh; export OMP_TARGET_OFFLOAD=MANDATORY ONEAPI_DEVICE_SELECTOR=level_zero:0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# MKL_NUM_THREADS=1; env')
INTEL_ENV = ("CLASSPATH CMAKE_PREFIX_PATH CMPLR_ROOT CPATH CPLUS_INCLUDE_PATH C_INCLUDE_PATH DIAGUTIL_PATH HOME I_MPI_ROOT LANG LC_ALL LD_LIBRARY_PATH "
             "LIBRARY_PATH MANPATH MKLROOT MKL_NUM_THREADS NEO_CACHE_PERSISTENT NLSPATH OCL_ICD_FILENAMES OMP_NUM_THREADS OMP_TARGET_OFFLOAD "
             "ONEAPI_DEVICE_SELECTOR OPENBLAS_NUM_THREADS PATH PKG_CONFIG_PATH PWD SHLVL TBBROOT TCMROOT TCM_ROOT TMPDIR UMF_ROOT USER _").split()


def synth_spec():
    src = QG.sha(os.path.join(QG.TOOLS, "t207_d5_qref.F90"))
    return {"spec_version": 1, "path": "def", "host_member": "nvh", "members": {"main": _mem(NUL_DEF), "ctl_ns0": _mem(NUL_DEF, False)},
            "refs": ["cpu", "nvh", "n4", "n100004", "n1", "n100001", "n2", "n1608", "n4808", "n5", "n100005"], "qs": 20, "nul": NUL_DEF,
            "cover": 1, "keep": 39, "hostin": 0, "min_reprocessed": 2, "expect_points": {"pTT": 60}, "prog_sha": "@harness", "src_sha": src,
            "binary_sha256": "@harness", "parameters_sha256": PAR, "seeds": "@harness",
            "env_roles": {"main": {"ALF_GPU_DEBUG": "1", "ALF_GPU_TAUDIAG_QHOST": "1", "ALF_GPU_TAUDIAG_QVAR": "8", "ALF_GPU_TAUDIAG_QREF": "3",
                                   "ALF_GPU_TAUDIAG_FROM": "11"},
                          "ctl_ns0": {"extra_env": {"ALF_GPU_TAU_NS": "0"}, "omit_env": ["ALF_GPU_TAUDIAG_QVAR"]}},
            "n_runs": {"main": 2, "ctl_ns0": 1}, "gates": [{"name": "acc", "measure": "pTT", "statistic": "S_ens", "pool": "pooled", "bound": 0.15}],
            "controls": [{"name": "var8", "kind": "variant", "member": "var8", "trips": "acc"},
                         {"name": "ns0", "kind": "run", "role": "ctl_ns0", "trips": "acc"}]}


def with_sha(g):
    g = copy.deepcopy(g)
    g["spec_template_sha256"] = hashlib.sha256(QG.canon(g["spec"])).hexdigest()
    return g


def part_a(V, tmp, out):
    import h5py
    import numpy as np
    main_env = ["ALF_GPU_DEBUG=1", "ALF_GPU_TAUDIAG_QHOST=1", "ALF_GPU_TAUDIAG_QVAR=8", "ALF_GPU_TAUDIAG_QREF=3", "ALF_GPU_TAUDIAG_FROM=11"]
    ctl_env = ["ALF_GPU_DEBUG=1", "ALF_GPU_TAUDIAG_QHOST=1", "ALF_GPU_TAUDIAG_QREF=3", "ALF_GPU_TAUDIAG_FROM=11", "ALF_GPU_TAU_NS=0"]
    qe = {"path": "def", "qs": 20, "nul": ",".join(str(c) for c in NUL_DEF), "hostin": 0, "keep": 39}
    acc = {"_status": "", "builds": {"gpu": {"tool": "nvhpc", "configure": "PGI noMPI HDF5 GPU=nvidia", "gpu": True}},
           "qbuilds": {"def": {"configure": "GNU noMPI", "qflags": "-DQR_PREPIVOT", "flags_ext": "-DQR_PREPIVOT"}},
           "runs": {"m1": {"build": "gpu", "config": "c", "seeds": "acc1", "env": main_env, "allow_debug": True, "qref": qe},
                    "m2": {"build": "gpu", "config": "c", "seeds": "acc2", "env": main_env, "allow_debug": True, "qref": qe},
                    "c1": {"build": "gpu", "config": "c", "seeds": "acc3", "env": ctl_env, "allow_debug": True, "qref": qe}},
           "controls": []}
    tools = {f: QG.sha(os.path.join(QG.TOOLS, f)) for f in QG.TOOL_FILES}
    G0 = with_sha({"id": "qref_def", "type": "qref", "spec": synth_spec(), "runs": ["m1", "m2"], "ctl_runs": {"ctl_ns0": ["c1"]}, "tools_sha256": tools})

    def lint_case(name, gate, sub, a=None):
        e = QG.lint(V, dict(a or acc, gates=[gate]))
        if sub is None:
            ok = not e
        else:
            ok = any(sub in x for x in e)
        print(f"{'ok ' if ok else 'BAD'} qref lint: {name:52s} expects {repr(sub) if sub else 'no error'} {'' if ok else e[:3]}", flush=True)
        out.append(("lint " + name, sub, e[:2], ok))

    lint_case("clean synthetic file", G0, None)
    g = copy.deepcopy(G0); g["spec"]["gates"][0]["bound"] = 0.2
    lint_case("spec edited, template sha not updated", g, "spec_template_sha256")
    g = with_sha(dict(G0, spec=dict(G0["spec"], test_only="x")))
    lint_case("test_only spec in a frozen file", g, "test_only")
    lint_case("run keep differs from the spec", G0, "differs from the spec's",
              a=dict(acc, runs=dict(acc["runs"], m2=dict(acc["runs"]["m2"], qref=dict(qe, keep=40)))))
    lint_case("main run with a control env", with_sha(dict(G0, runs=["m1", "c1"], ctl_runs={"ctl_ns0": ["m2"]})), "listed for role")
    lint_case("QREF run declaring ALF_LA_TEST_X", G0, "only ALF_GPU_",
              a=dict(acc, runs=dict(acc["runs"], m1=dict(acc["runs"]["m1"], env=main_env + ["ALF_LA_TEST_X=1"]))))
    sp = copy.deepcopy(G0["spec"]); sp["env_roles"]["main"]["ALF_GPU_TAUDIAG_QREF"] = "2"
    lint_case("stride 2", with_sha(dict(G0, spec=sp)), "decision 35")
    sp = copy.deepcopy(G0["spec"]); sp["n_runs"]["main"] = 3
    lint_case("n_runs differs from the runs", with_sha(dict(G0, spec=sp)), "n_runs")
    lint_case("qref_env_allowed with LIBOMPTARGET_DEBUG", G0, "may not be allowed", a=dict(acc, qref_env_allowed=["NVHPC", "LIBOMPTARGET_DEBUG"]))
    sp = copy.deepcopy(G0["spec"]); sp["refs"] = sp["refs"][:-1]
    lint_case("unfrozen reference list", with_sha(dict(G0, spec=sp)), "207-2 gate refuses")
    lint_case("control replacing 'device'", G0, "may replace only runs",
              a=dict(acc, controls=[{"id": "c", "like": "qref_def", "replace": {"device": "var8"}, "expect": "FAIL"}]))
    lint_case("tools_sha256 incomplete", dict(G0, tools_sha256={k: v for k, v in tools.items() if k != "t207_d5_gate.py"}), "tools_sha256")
    sp = copy.deepcopy(G0["spec"]); del sp["prog_sha"]
    lint_case("instance field missing", with_sha(dict(G0, spec=sp)), "spec.prog_sha missing")
    # fix round 3
    sp = copy.deepcopy(G0["spec"]); sp["parameters_sha256"] = "@harness"
    lint_case("parameters_sha256 @harness in a frozen file", with_sha(dict(G0, spec=sp)), "parameters_sha256 must be concrete")
    lint_case("parameters_sha256 @harness in a DRAFT", with_sha(dict(G0, spec=sp)), None, a=dict(acc, _status="DRAFT x"))
    devr = {r: dict(x, seeds="dev%d" % (i + 1)) for i, (r, x) in enumerate(acc["runs"].items())}
    lint_case("@harness seeds with dev seeds in a frozen file", G0, "acceptance (acc*) seeds", a=dict(acc, runs=devr))
    lint_case("@harness seeds with dev seeds in a DRAFT", G0, None, a=dict(acc, runs=devr, _status="DRAFT x"))
    sp = copy.deepcopy(G0["spec"])
    sp["seeds"] = {"main": ["1" * 64, "2" * 64], "ctl_ns0": ["3" * 64]}
    lint_case("concrete seeds with dev seeds in a frozen file", with_sha(dict(G0, spec=sp)), None, a=dict(acc, runs=devr))
    swr = dict(acc["runs"], m1s=dict(acc["runs"]["m1"]), m2s=dict(acc["runs"]["m2"]))
    swm = [{"src": "m1", "dst": "m1s", "kind": "qref_swap:dev:var8"}, {"src": "m2", "dst": "m2s", "kind": "qref_swap:dev:var8"}]
    swc = [{"id": "c", "like": "qref_def", "replace": {"runs": ["m1s", "m2s"]}, "expect": "FAIL"}]
    lint_case("harness control: qref_swap copies of both main runs", G0, None, a=dict(acc, runs=swr, mutations=swm, controls=swc))
    lint_case("harness control: a copy of another kind", G0, "must be a mutations entry of kind",
              a=dict(acc, runs=swr, mutations=[swm[0], dict(swm[1], kind="seeds")], controls=swc))
    lint_case("harness control: both copies of one main run", G0, "exactly the gate's main runs",
              a=dict(acc, runs=swr, mutations=[swm[0], dict(swm[1], src="m1")], controls=swc))
    lint_case("harness control copy with dev seeds, @harness, frozen", G0, "acceptance (acc*) seeds",
              a=dict(acc, runs=dict(swr, m2s=dict(swr["m2s"], seeds="dev2")), mutations=swm, controls=swc))
    lint_case("QREF run on a base build", G0, "base build",
              a=dict(acc, builds=dict(acc["builds"], base_gpu=acc["builds"]["gpu"]), runs=dict(acc["runs"], m2=dict(acc["runs"]["m2"], build="base_gpu"))))
    pinr = "intel.sh sets it (no kernel cache outside ~/alf-gpu)"
    lint_case("qref_env_pinned NEO_CACHE_PERSISTENT=0 with reason", G0, None,
              a=dict(acc, qref_env_pinned={"NEO_CACHE_PERSISTENT": "0"}, _qref_env_pinned_reason=pinr))
    lint_case("qref_env_pinned without reason", G0, "_qref_env_pinned_reason", a=dict(acc, qref_env_pinned={"NEO_CACHE_PERSISTENT": "0"}))
    for n in ("OMP_TARGET_OFFLOAD", "ONEAPI_DEVICE_SELECTOR", "ZE_AFFINITY_MASK", "LIBOMPTARGET_DEBUG", "ALF_GPU_TAU_NS", "LD_PRELOAD", "OMP_NUM_THREADS"):
        lint_case(f"qref_env_pinned {n}", G0, "may not be pinned", a=dict(acc, qref_env_pinned={n: "0"}, _qref_env_pinned_reason=pinr))
    lint_case("qref_env_allowed NEO_CACHE_PERSISTENT (by name)", G0, "may not be allowed", a=dict(acc, qref_env_allowed=["NEO_CACHE_PERSISTENT"]))
    lint_case("qref_gate_timeout_s 0", G0, "qref_gate_timeout_s", a=dict(acc, qref_gate_timeout_s=0))

    # check_env with the complete environment of a GPU=intel run on `intel` (names read on intel 2026-11-01: env -i + ~/alf-gpu/env/intel.sh
    # + the harness exports, as acc_run.sh builds it); only the names matter, values are placeholders except the checked ones
    ienv = {n: "x" for n in INTEL_ENV}
    ienv.update(OMP_TARGET_OFFLOAD="MANDATORY", ONEAPI_DEVICE_SELECTOR="level_zero:0", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                MKL_NUM_THREADS="1", NEO_CACHE_PERSISTENT="0", **dict(kv.split("=", 1) for kv in main_env))
    iacc = dict(acc, qref_env_allowed=sorted(n for n in INTEL_ENV if n not in QG.ENV_BASE and n != "NEO_CACHE_PERSISTENT"),
                qref_env_pinned={"NEO_CACHE_PERSISTENT": "0"}, _qref_env_pinned_reason=pinr,
                builds={"gpu": {"tool": "intel", "configure": "IntelLLVM noMPI HDF5 GPU=intel", "gpu": True}})
    lint_case("intel: qref_env_allowed + pinned for the full intel env", G0, None, a=iacc)

    def env_case(name, env, a, why):
        try:
            QG.check_env(V, "m1", {"env": env}, acc["runs"]["m1"], a, a["builds"]["gpu"])
            v, w = "PASS", None
        except V.InputFail as e:
            v, w = "INPUT_ERROR", str(e)
        _res(V, name, "PASS" if why is None else "INPUT_ERROR", v, w, why, out)
    env_case("check_env: full intel env, NEO_CACHE_PERSISTENT=0 pinned", ienv, iacc, None)
    env_case("check_env: full intel env without the pin", ienv, dict(iacc, qref_env_pinned={}), "cannot see")
    env_case("check_env: NEO_CACHE_PERSISTENT=1 (wrong value)", dict(ienv, NEO_CACHE_PERSISTENT="1"), iacc, "pre-registered value")
    env_case("check_env: NEO_CACHE_PERSISTENT missing", {k: v for k, v in ienv.items() if k != "NEO_CACHE_PERSISTENT"}, iacc, "pre-registered value")
    env_case("check_env: intel env + SYCL_PI_TRACE", dict(ienv, SYCL_PI_TRACE="1"), iacc, "cannot see")

    # h5_identical on synthetic run directories (unchanged gate of the first round)
    root = os.path.join(tmp, "runs_h5")
    os.makedirs(root)
    hacc = {"candidate_sha": CAND[:12], "gpu_exclusive": True, "active_lines": ["GPU tau pass (207): ACTIVE"],
            "builds": {"gpu": {"tool": "nvhpc", "configure": "PGI noMPI HDF5 GPU=nvidia", "gpu": True}}, "runs": {}}

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
             "tar_sha256_built": TAR, "candidate": {"sha": CAND, "tar_sha256": TAR}}
        json.dump({"schema": "alf-acc-prov/1", "binary_sha256": "b" * 64, "build": b, "exit_code": 0, "timed_out": False, "extra_env": list(env),
                   "debug_env_present": [t.split("=")[0] for t in env], "env": {"OMP_TARGET_OFFLOAD": "MANDATORY"}, "gpu_check_rc": 0, "load_ok": True,
                   "gpu_before": [], "inputs": {"parameters_sha256": "p", "seeds_sha256": "s"},
                   "outputs": {f: hh(f) for f in ("info", "data.h5", "run.log")}}, open(os.path.join(d, "provenance.json"), "w"))
        hacc["runs"][rid] = {"build": "gpu", "config": "synth", "env": list(env), "allow_debug": True, "active_lines": ["GPU tau pass (207): ACTIVE"]}

    def hcase(name, expect, gate, why=None):
        v, r, w = V.run_gate(gate, root, hacc)
        _res(V, name, expect, v, w, why, out)

    rng = np.random.default_rng(5)
    data = {"Ener_scal/obser": rng.normal(size=(8, 2)), "Green_eq/obser": rng.normal(size=(8, 3))}
    mkh5("h_on", main_env, data)
    mkh5("h_off", [], data)
    GI = {"id": "qref_identity", "type": "h5_identical", "pairs": [["h_on", "h_off"]]}
    hcase("identity: data.h5 QREF on == off", "PASS", GI)
    d2 = {k: x.copy() for k, x in data.items()}
    d2["Green_eq/obser"][3, 1] = np.nextafter(d2["Green_eq/obser"][3, 1], 1e9)
    mkh5("h_off_ulp", [], d2)
    hcase("identity: one entry one ulp different", "FAIL", dict(GI, pairs=[["h_on", "h_off_ulp"]]), why="differ")
    mkh5("h_off_dbg", ["ALF_GPU_DEBUG=1"], data)
    hcase("identity: the 'off' run has a debug variable", "INPUT_ERROR", dict(GI, pairs=[["h_on", "h_off_dbg"]]), why="'off' run")
    hcase("identity: the 'on' run has no QREF", "INPUT_ERROR", dict(GI, pairs=[["h_off", "h_off"]]), why="'on' run")


# ------------------------------------------------------------------ part B: against the real 207-2 gate
def _fixture(V, C):
    cand = dict(l.split(None, 1) for l in open(os.path.join(C, ".alf_candidate")).read().splitlines() if " " in l)
    parts = [l.split()[1] for l in open(os.path.join(C, ".alf_candidate")).read().splitlines() if l.startswith("acceptance_part ")]
    for name in parts:
        acc = json.load(open(os.path.join(C, "src", "testsuite", "gpu_reference", "acceptance", name)))
        g = next((x for x in acc.get("gates", []) if x.get("type") == "qref"), None)
        if g is not None:
            acc["candidate_sha"] = cand["sha"]
            acc["_seeds_record"] = json.load(open(os.path.join(C, "configs", "seeds_record.json")))
            return acc, g, cand
    raise SystemExit(f"fixture {C}: no acceptance part with a qref gate")


def _lines(p):
    return open(p).read().split("\n")[:-1]


def _write_lines(p, ls):
    open(p, "w").write("\n".join(ls) + "\n")


def _fix_qsha(q):
    pj = os.path.join(q, "qref_proc.json")
    pr = json.load(open(pj))
    pr["qref_sha256"] = QG.sha(os.path.join(q, "qref.txt"))
    pr["manifest_sha256"] = QG.sha(os.path.join(q, "qref_manifest.txt"))
    json.dump(pr, open(pj, "w"), indent=1)


def _kept(q):
    return sorted(int(f[1:7]) for f in os.listdir(os.path.join(q, "qref")) if f.endswith(".bin"))


def edit_line(rid, pick, fn):
    """qref.txt of run rid: the line of a kept dump (pick 'kept') or of a deleted one (pick 'deleted') edited by fn; qref_sha256 follows"""
    def m(runs):
        q = os.path.join(runs, rid, "qproc")
        ls = _lines(os.path.join(q, "qref.txt"))
        kept = _kept(q)
        idx = [i for i, l in enumerate(ls) if l.startswith("QREF file ") and ((int(l.split()[2]) in kept) == (pick == "kept"))]
        i = idx[len(idx) // 2]
        new = fn(ls[i])
        assert new != ls[i], "edit had no effect"
        if new is None:
            del ls[i]
        else:
            ls[i] = new
        _write_lines(os.path.join(q, "qref.txt"), ls)
        _fix_qsha(q)
    return m


def tok_edit(member, k, val):
    def f(line):
        segs = line.split(" | ")
        j = next(i for i, s in enumerate(segs) if s.split()[:1] == [member])
        t = segs[j].split()
        t[1 + k] = val(t[1 + k]) if callable(val) else val
        segs[j] = " ".join(t)
        return " | ".join(segs)
    return f


def part_b(V, C, tmp, out):
    acc, g0, cand = _fixture(V, C)
    spec0 = g0["spec"]
    rids = list(g0["runs"]) + [r for v in g0["ctl_runs"].values() for r in v]
    mains = list(g0["runs"])
    swaps = {}
    for c in acc.get("controls", []):
        if c.get("like") == g0["id"]:
            swaps = {r: r2 for r, r2 in zip(mains, c["replace"]["runs"])}
    allr = rids + list(swaps.values())
    n = [0]

    def setup(mutate=None, restamp=True, extra_lines=None, drop_stamp=None, stamp_cand=None):
        n[0] += 1
        T = os.path.join(tmp, "b%03d" % n[0])
        runs = os.path.join(T, "runs")
        os.makedirs(os.path.join(T, "bin"))
        os.makedirs(os.path.join(T, "out"))
        for r in allr:
            shutil.copytree(os.path.join(C, "runs", r), os.path.join(runs, r), symlinks=True)
        for f in os.listdir(os.path.join(C, "bin")):
            if f.startswith("qref_") and f.endswith(".qbuild.json"):
                shutil.copy2(os.path.join(C, "bin", f), os.path.join(T, "bin", f))

        def stamps():
            ls = []
            for r in allr:
                if r == drop_stamp:
                    continue
                ls.append("2026-11-01T00:00:00 " + QG.stamp_line("agent", stamp_cand or cand["sha"], r, QG.stamp_object(os.path.join(runs, r))))
            return ls
        lines = None if restamp else stamps()
        if mutate:
            mutate(runs)
        if restamp:
            lines = stamps()
        return T, runs, lines + (extra_lines(runs) if extra_lines else [])

    def case(name, expect, gate=None, why=None, mutate=None, restamp=True, extra_lines=None, drop_stamp=None, a=None, stamp_cand=None):
        T, runs, lines = setup(mutate, restamp, extra_lines, drop_stamp, stamp_cand)
        a2 = dict(acc, _attempt_lines=lines, _out=os.path.join(T, "out"), **dict({"_allow_mutated": True}, **(a or {})))
        v, r, w = V.run_gate(gate or g0, runs, a2)
        _res(V, name, expect, v, w, why, out)
        shutil.rmtree(T)
        return v, r, w

    t0 = time.time()
    T, runs, lines = setup()
    v0, r0, w0 = V.run_gate(g0, runs, dict(acc, _attempt_lines=lines, _out=os.path.join(T, "out")))
    rep = json.load(open(r0["report"])) if r0 else None
    if rep is None:  # FAIL / INPUT_ERROR: the report path is in the reason
        m = re.search(r"report (\S+\.report\.json)", w0 or "")
        rep = json.load(open(m.group(1))) if m else None
    ok = v0 in ("PASS", "FAIL") and rep is not None
    print(f"{'ok ' if ok else 'BAD'} qref {'clean fixture, pre-registered spec':58s} PASS or FAIL (informational) got {v0}: {(w0 or '')[:200]}", flush=True)
    out.append(("clean fixture", "PASS|FAIL", v0, ok))
    if rep is None:
        return
    # points per run for every measure, from the fixture's qref.txt (pTT: dumps with B; emax / eF: all)
    q1 = os.path.join(runs, mains[0], "qproc", "qref.txt")
    allp = sum(1 for l in _lines(q1) if l.startswith("QREF file "))
    withb = sum(1 for l in _lines(q1) if l.startswith("QREF file ") and " hasb 1 " in l)
    npts = {"pTT": withb, "emax": allp, "eF": allp}
    shutil.rmtree(T)

    def variant(gates, controls=None, runs_=None, **kw):
        sp = copy.deepcopy(spec0)
        sp["gates"] = gates
        sp["expect_points"] = {x["measure"]: npts[x["measure"]] for x in gates}
        sp["controls"] = controls if controls is not None else [dict({k: v for k, v in c.items() if k != "must_not_trip"}, trips=gates[0]["name"])
                                                                  for c in spec0["controls"]]
        sp.update(kw)
        return with_sha(dict(g0, spec=sp, runs=runs_ or g0["runs"]))

    # a separable (measure, statistic) on the original main runs for the PASS path: device below every control.  Probe gates on every measure
    # (bounds +50: the numbers only, from the report's localized section); the qref_swap copies (device and var8 exchanged: the copies' device is
    # the originals' var8) are the planted counterpart of the FAIL cases (fix round 3: only the original runs can be the clean gate of the
    # harness control, whose runs must be qref_swap copies).
    ctl_who = [("main", c["member"]) if c["kind"] == "variant" else (c["role"], "dev") for c in spec0["controls"]]
    sets = [("original main runs", mains)]
    sep = None
    for label, mset in sets:
        probe = variant([{"name": "p_" + mm, "measure": mm, "statistic": "S_ens", "pool": "pooled", "bound": 50.0} for mm in ("pTT", "emax", "eF")],
                        runs_=mset)
        T, runs, lines = setup()
        vp, rp_, wp = V.run_gate(probe, runs, dict(acc, _attempt_lines=lines, _out=os.path.join(T, "out"), _allow_mutated=True))
        m_ = re.search(r"report (\S+\.report\.json)", wp or "")
        rpp = rp_["report"] if rp_ else (m_.group(1) if m_ else None)
        loc = json.load(open(rpp)).get("localized", {}) if rpp else {}
        shutil.rmtree(T)
        print(f"    probe on the {label}: gate {vp}; " + "; ".join(
            f"{mm}: " + " ".join(f"{e['role']}/{e['member']} {e['pooled']['mean']:+.3f}" for e in loc.get(mm, [])) for mm in sorted(loc)), flush=True)

        def val(mm, role, member, stat):
            e = next(x for x in loc[mm] if x["role"] == role and x["member"] == member)
            return e["pooled"]["mean" if stat == "S_ens" else "q95" if stat == "S_ens_q95" else "max"]
        for mm in sorted(loc):
            for st in ("S_ens", "S_ens_max", "S_ens_q95"):
                d = val(mm, "main", "dev", st)
                cmin = min(val(mm, r, x, st) for r, x in ctl_who)
                if cmin > d + 1e-3:
                    sep = (label, mset, mm, st, d, cmin)
                    break
            if sep:
                break
        if sep:
            break
    if sep is None:
        print("BAD qref no separable (measure, statistic) on this fixture: the PASS path cannot be shown (not a code defect; another fixture is needed)")
        out.append(("separable fixture", "exists", None, False))
        return
    label, mset, m, st, d, cmin = sep
    other = [x for x in (mains + [swaps[r] for r in mains] if swaps else mains) if x not in mset]
    MU = {"_allow_mutated": True}          # the qref_swap copies are mutated copies (allowed in controls only; here the self-check's cases)
    b = (d + cmin) / 2
    gP = variant([{"name": "sc", "measure": m, "statistic": st, "pool": "pooled", "bound": b}], runs_=mset)
    print(f"    PASS-path spec on the {label}: {st}({m}) device {d:+.4f}, smallest control {cmin:+.4f}, bound {b:+.4f} (self-check only, no calibration)")
    vP, rP, _w = case("separable spec (bound between device and controls)", "PASS", gP, a=MU)
    DEVFAIL = re.compile(r"device [+-][0-9]+\.[0-9]+ > bound [+-][0-9]+\.[0-9]+")
    case("bound below the device value", "FAIL", variant([{"name": "sc", "measure": m, "statistic": st, "pool": "pooled", "bound": d - 0.05}], runs_=mset),
         why=DEVFAIL, a=MU)
    case("bound above every control (controls not killed)", "FAIL",
         variant([{"name": "sc", "measure": m, "statistic": st, "pool": "pooled", "bound": cmin + 0.05}], runs_=mset), why="not killed", a=MU)
    m2 = "eF" if m != "eF" else "emax"
    g2 = variant([{"name": "sc", "measure": m, "statistic": st, "pool": "pooled", "bound": b},
                  {"name": "pw", "measure": m2, "statistic": "S_ens_max", "pool": "pooled", "bound": -50.0}], runs_=mset)
    case("defect in one of two statistics (second bound below its device value)", "FAIL", g2, why="gate pw", a=MU)
    g3 = variant([{"name": "sc", "measure": m, "statistic": st, "pool": "pooled", "bound": b},
                  {"name": "nc", "measure": m, "statistic": st, "pool": "pooled", "bound": 50.0}], runs_=mset)
    case("a gate no control trips: DEV-ONLY", "INPUT_ERROR", g3, why="DEV-ONLY", a=MU)
    case("spec marked test_only: TEST-ONLY", "INPUT_ERROR", variant(gP["spec"]["gates"], runs_=mset, test_only="self-check"), why="TEST-ONLY", a=MU)
    if swaps and rP:
        # the harness control (fix round 3): a FAIL counts only through the device (> bound), with the copies' device = the clean var8
        GR = {"_gate_results": {gP["id"]: rP}}
        ctl = dict(gP, runs=other, _control_of=gP["id"])
        CTLDET = re.compile(r"^control detection: gate sc: device [+-][0-9]+\.[0-9]+ > bound [+-][0-9]+\.[0-9]+ \(= the clean report's var8 value\)")
        case("harness control: qref_swap copies vs the clean gate (detection)", "FAIL", ctl, why=CTLDET, a=dict(MU, **GR))
        case("the same without the clean result (swap gate alone)", "INPUT_ERROR", ctl, why="uninformative control", a=MU)
        gH = variant([{"name": "sc", "measure": m, "statistic": st, "pool": "pooled", "bound": 50.0}], runs_=other)
        case("harness control failing only by the internal var8 (bound 50)", "INPUT_ERROR", dict(gH, _control_of=gP["id"]),
             why="without the device failing", a=dict(MU, **GR))
        bad_rP = copy.deepcopy(rP)
        for c_ in bad_rP["controls"]:
            for ch in c_["checks"]:
                if c_["member"] == "var8":
                    ch["pooled"] = ch["pooled"] + 1e-3
        case("harness control: copies' device is not the clean var8", "INPUT_ERROR", ctl, why="is not the clean report's var8", a=dict(MU, _gate_results={gP["id"]: bad_rP}))
        gB = variant([{"name": "sc", "measure": m, "statistic": st, "pool": "pooled", "bound": d - 0.05}], runs_=mset)
        case("harness control on the original (unswapped) runs", "INPUT_ERROR", dict(gB, _control_of=gP["id"]), why="not qref_swap:dev:var8 copies",
             a=dict(MU, **GR))
        ev = V.evaluate
        T, runs, lines = setup()
        accE = dict(acc, gates=[gP], controls=[{"id": "c_swap", "like": gP["id"], "replace": {"runs": other}, "expect": "FAIL"}],
                    _attempt_lines=lines, _out=os.path.join(T, "out"), _status="DRAFT selfcheck", blocked=[], **MU)
        e = ev(accE, runs)
        okE = e["gates"][0]["verdict"] == "PASS" and e["controls"][0]["outcome"] == "FAIL" and e["controls"][0]["as_expected"] \
            and CTLDET.search(e["controls"][0]["reason"] or "") is not None
        print(f"{'ok ' if okE else 'BAD'} qref {'evaluate(): gate PASS, harness control (swapped copies) FAIL':58s} got {e['gates'][0]['verdict']} / "
              f"{e['controls'][0]['outcome']} | {(e['controls'][0]['reason'] or '')[:160]}"
              + ("" if okE else f": {e['gates'][0]['reason']} | {e['controls'][0]['reason']}"), flush=True)
        out.append(("evaluate with harness control", "PASS/FAIL", None, okE))
        shutil.rmtree(T)
    else:
        print("BAD qref the fixture has no qref_swap control: the swapped-copy cases cannot run")
        out.append(("swap control present", "yes", None, False))
    r1 = mset[0]

    def swap_text(line):
        segs = line.split(" | ")
        i = next(j for j, s in enumerate(segs) if s.split()[:1] == ["dev"])
        k = next(j for j, s in enumerate(segs) if s.split()[:1] == ["var8"])
        segs[i], segs[k] = "dev" + segs[k][4:], "var8" + segs[i][3:]
        return " | ".join(segs)

    def all_lines(rid, fn):
        def mm(runs):
            q = os.path.join(runs, rid, "qproc")
            _write_lines(os.path.join(q, "qref.txt"), [fn(l) if l.startswith("QREF file ") else l for l in _lines(os.path.join(q, "qref.txt"))])
            _fix_qsha(q)
        return mm

    case("device and var8 swapped in qref.txt only (restamped)", "INPUT_ERROR", gP, mutate=all_lines(r1, swap_text), why="re-comput")
    case("NaN in a deleted dump's line, not restamped", "INPUT_ERROR", gP, mutate=edit_line(r1, "deleted", tok_edit("dev", 0, "NaN")), restamp=False,
         why="changed after the stamp")
    case("NaN in a deleted dump's line (restamped)", "INPUT_ERROR", gP, mutate=edit_line(r1, "deleted", tok_edit("dev", 0, "NaN")), why="ES12.4E3")
    case("Infinity in a deleted dump's line (restamped)", "INPUT_ERROR", gP, mutate=edit_line(r1, "deleted", tok_edit("cpu", 1, "Infinity")), why="ES12.4E3")
    case("E-less number (restamped)", "INPUT_ERROR", gP, mutate=edit_line(r1, "deleted", tok_edit("n4", 2, lambda x: x.split("E")[0])), why="ES12.4E3")
    case("missing member n4 in one line (restamped)", "INPUT_ERROR", gP,
         mutate=edit_line(r1, "deleted", lambda l: " | ".join(s for s in l.split(" | ") if s.split()[:1] != ["n4"])), why="members")
    case("missing point (one line deleted, restamped)", "INPUT_ERROR", gP, mutate=edit_line(r1, "deleted", lambda l: None), why="QREF lines")

    def rot(member):
        def f(line):
            segs = line.split(" | ")
            j = next(i for i, s in enumerate(segs) if s.split()[:1] == [member])
            t = segs[j].split()
            t = [t[0]] + t[2:6] + [t[1]] + t[6:]
            segs[j] = " ".join(t)
            return " | ".join(segs)
        return f
    case("shifted fields in a kept dump's line (restamped)", "INPUT_ERROR", gP, mutate=edit_line(r1, "kept", rot("dev")))
    case("shifted fields in a deleted dump's line, not restamped", "INPUT_ERROR", gP, mutate=edit_line(r1, "deleted", rot("dev")), restamp=False,
         why="changed after the stamp")
    case("changed header (nfe) in a line (restamped)", "INPUT_ERROR", gP,
         mutate=edit_line(r1, "deleted", lambda l: re.sub(r" nfe ([0-9]+) ", lambda x: " nfe %d " % (3 - int(x.group(1))), l, count=1)), why="header")

    def prov_env(rid, k, v):
        def mm(runs):
            pj = os.path.join(runs, rid, "provenance.json")
            p = json.load(open(pj))
            p["env"][k] = v
            json.dump(p, open(pj, "w"), indent=1)
        return mm
    case("extra ALF_GPU_TAUSWAP in the recorded env (restamped)", "INPUT_ERROR", gP, mutate=prov_env(r1, "ALF_GPU_TAUSWAP", "8"), why="ALF_GPU")
    case("ALF_LA_TEST_X in the recorded env (restamped)", "INPUT_ERROR", gP, mutate=prov_env(r1, "ALF_LA_TEST_X", "1"), why="undeclared debug")
    case("ALFSH_X in the recorded env (restamped)", "INPUT_ERROR", gP, mutate=prov_env(r1, "ALFSH_X", "1"), why="undeclared debug")
    case("OMP_TARGET_OFFLOAD=DISABLED in the recorded env (restamped)", "INPUT_ERROR", gP, mutate=prov_env(r1, "OMP_TARGET_OFFLOAD", "DISABLED"),
         why="MANDATORY")
    case("LIBOMPTARGET_DEBUG in the recorded env (restamped)", "INPUT_ERROR", gP, mutate=prov_env(r1, "LIBOMPTARGET_DEBUG", "1"), why="cannot see")
    case("LD_PRELOAD in the recorded env (restamped)", "INPUT_ERROR", gP, mutate=prov_env(r1, "LD_PRELOAD", "/x.so"), why="cannot see")
    case("unknown variable FOO in the recorded env (restamped)", "INPUT_ERROR", gP, mutate=prov_env(r1, "FOO", "1"), why="cannot see")
    case("changed script hash (tools_sha256)", "INPUT_ERROR", dict(gP, tools_sha256=dict(gP["tools_sha256"], **{"t207_d5_gate.py": "0" * 64})),
         why="pre-registered script")
    sp = copy.deepcopy(gP["spec"]); sp["prog_sha"] = "1" * 64
    case("pre-registered prog_sha differs from the offline build", "INPUT_ERROR", with_sha(dict(gP, spec=sp)), why="prog_sha")

    def qbuild_edit(runs):
        f = [x for x in os.listdir(os.path.join(runs, "..", "bin")) if x.endswith(".qbuild.json")][0]
        p = os.path.join(runs, "..", "bin", f)
        d = json.load(open(p))
        d["prog_sha256"] = "2" * 64
        json.dump(d, open(p, "w"), indent=1)
    case("offline build record edited (program hash)", "INPUT_ERROR", gP, mutate=qbuild_edit, why="build record")

    def rec_edit(runs):
        p = os.path.join(runs, r1, "qproc", "info")
        open(p, "a").write(" \n")
    case("altered record qproc/info (restamped)", "INPUT_ERROR", gP, mutate=rec_edit, why="harness provenance")

    def del_dump(runs):
        q = os.path.join(runs, r1, "qproc", "qref")
        os.remove(os.path.join(q, sorted(os.listdir(q))[-1]))
    case("a kept dump deleted, not restamped", "INPUT_ERROR", gP, mutate=del_dump, restamp=False, why="files differ from the stamp")
    case("a kept dump deleted (restamped)", "INPUT_ERROR", gP, mutate=del_dump, why="qref/ holds")

    def sym_dump(runs):
        q = os.path.join(runs, r1, "qproc", "qref")
        f = sorted(os.listdir(q))[0]
        dst = os.path.join(runs, "..", "elsewhere_" + f)
        shutil.move(os.path.join(q, f), dst)
        os.symlink(os.path.abspath(dst), os.path.join(q, f))
    case("a kept dump replaced by a symbolic link (same bytes)", "INPUT_ERROR", gP, mutate=sym_dump, restamp=False, why="symbolic link")

    def mach_edit(runs):
        p = os.path.join(runs, r1, "qproc", "acc_qproc.json")
        d = json.load(open(p))
        d["machine"]["host"] = "elsewhere"
        json.dump(d, open(p, "w"), indent=1)
    case("processed on another machine (restamped)", "INPUT_ERROR", gP, mutate=mach_edit, why="same machine")
    case("no stamp for a run", "INPUT_ERROR", gP, drop_stamp=r1, why="0 qref_stamp")
    case("two stamps for a run", "INPUT_ERROR", gP,
         extra_lines=lambda runs: ["2026-11-01T00:00:01 " + QG.stamp_line("agent", cand["sha"], r1, QG.stamp_object(os.path.join(runs, r1)))],
         why="2 qref_stamp")
    case("stamp of another candidate", "INPUT_ERROR", gP, a={"candidate_sha": "deadbeef0000" + "0" * 28}, why="built from")
    g5 = copy.deepcopy(gP); g5["spec"]["gates"][0]["bound"] = b + 0.01
    case("spec differs from spec_template_sha256", "INPUT_ERROR", g5, why="spec_template_sha256")
    ctl0 = next(iter(g0["ctl_runs"].values()))[0]
    case("control run listed as a main run", "INPUT_ERROR", dict(gP, runs=[mset[0], ctl0], ctl_runs={k: [mset[1]] for k in g0["ctl_runs"]}),
         why="role")
    case("a run listed twice", "INPUT_ERROR", dict(gP, runs=[mset[0], mset[0]]), why="listed twice")

    # ---- fix round 3
    def jedit(rel, fn):
        def mm(runs):
            pj = os.path.join(runs, r1, rel)
            d = json.load(open(pj))
            fn(d)
            json.dump(d, open(pj, "w"), indent=1)
        return mm

    def base_build(p):
        p["build"]["source"] = "base_src"
        p["build"]["candidate"] = dict(p["build"]["candidate"], base_sha=cand["sha"], base_tar_sha256=p["build"]["tar_sha256_built"])
    case("binary from a base build record (source base_src, consistent)", "INPUT_ERROR", gP, mutate=jedit("provenance.json", base_build),
         why="never the base", a={"base_sha": cand["sha"]})
    case("run's env file differs from its build record's", "INPUT_ERROR", gP, mutate=jedit("provenance.json", lambda p: p.update(env_file_sha256="3" * 64)),
         why="env file")
    case("acc_qproc.json proc_rc 1 (restamped)", "INPUT_ERROR", gP, mutate=jedit("qproc/acc_qproc.json", lambda d: d.update(proc_rc=1)),
         why="processing status 1")
    case("qproc/ missing (restamped)", "INPUT_ERROR", gP, mutate=lambda runs: shutil.rmtree(os.path.join(runs, r1, "qproc")), why="no qproc/")
    case("bridge qproc/provenance.json edited (seeds, restamped)", "INPUT_ERROR", gP,
         mutate=jedit("qproc/provenance.json", lambda d: d.update(seeds_sha256="4" * 64)), why="qproc/provenance.json is not")

    def stamp_txt(runs):
        p = os.path.join(runs, r1, "qproc", "stamp.txt")
        ls = open(p).read().splitlines()
        ls[1] = ls[1] + " ALF_GPU_TAU_NS=0"
        open(p, "w").write("\n".join(ls) + "\n")
    case("bridge qproc/stamp.txt edited (env line, restamped)", "INPUT_ERROR", gP, mutate=stamp_txt, why="stamp.txt is not generated")
    case("qref_proc.json with another python (restamped)", "INPUT_ERROR", gP, mutate=jedit("qproc/qref_proc.json", lambda d: d.update(python="/usr/bin/python3.x")),
         why="interpreter")

    def lib_edit(d):
        k = sorted(d["machine"]["libs"])[0]
        d["machine"]["libs"][k] = dict(d["machine"]["libs"][k], sha256="6" * 64)
    case("acc_qproc.json: one library of the program differs (restamped)", "INPUT_ERROR", gP, mutate=jedit("qproc/acc_qproc.json", lib_edit),
         why="same machine")
    case("stamps written for another candidate (find_stamp)", "INPUT_ERROR", gP, stamp_cand="deadbeef0000" + "0" * 28, why="stamp of candidate deadbeef0000")
    case("@harness seeds without the driver's seeds record", "INPUT_ERROR", gP, a={"_seeds_record": None}, why="no seeds record")
    case("@harness seeds, frozen file, a main run on dev seeds", "INPUT_ERROR", gP,
         a={"_status": "", "runs": dict(acc["runs"], **{r1: dict(acc["runs"][r1], seeds="dev1")})}, why="outside a DRAFT")
    sp = copy.deepcopy(gP["spec"]); sp["parameters_sha256"] = rP["filled"]["parameters_sha256"] if rP else "0" * 64
    case("concrete parameters_sha256 = the runs' (frozen-file form)", "PASS", with_sha(dict(gP, spec=sp)), a=MU)
    sp2 = copy.deepcopy(sp); sp2["parameters_sha256"] = "1" * 64
    case("concrete parameters_sha256 differs from the runs'", "INPUT_ERROR", with_sha(dict(gP, spec=sp2)), why="parameters_sha256", a=MU)
    PIN = {"qref_env_pinned": {"NEO_CACHE_PERSISTENT": "0"}, "_qref_env_pinned_reason": "self-check"}
    case("pinned NEO_CACHE_PERSISTENT=0 in every run's recorded env (restamped)", "PASS", gP,
         mutate=lambda runs: [prov_env(r, "NEO_CACHE_PERSISTENT", "0")(runs) for r in rids], a=dict(MU, **PIN))
    case("pinned NEO_CACHE_PERSISTENT=0 in one run only (restamped)", "INPUT_ERROR", gP, mutate=prov_env(r1, "NEO_CACHE_PERSISTENT", "0"),
         why="NEO_CACHE_PERSISTENT=None", a=dict(MU, **PIN))
    case("pinned NEO_CACHE_PERSISTENT=1 in the recorded env (restamped)", "INPUT_ERROR", gP, mutate=prov_env(r1, "NEO_CACHE_PERSISTENT", "1"),
         why="pre-registered value", a=dict(MU, **PIN))
    case("NEO_CACHE_PERSISTENT=0 recorded, not pinned (restamped)", "INPUT_ERROR", gP, mutate=prov_env(r1, "NEO_CACHE_PERSISTENT", "0"), why="cannot see",
         a=MU)
    print(f"    part B: {n[0]} case directories, {time.time() - t0:.0f} s")


def selfcheck(V, fixture=None):
    t0 = time.time()
    tmp = V.tmpdir()
    out = []
    part_a(V, tmp, out)
    ran = False
    if fixture:
        part_b(V, os.path.abspath(fixture), tmp, out)
        ran = True
    else:
        print("    qref part B (cases against the real 207-2 gate) NOT RUN: no --qref-fixture")
    shutil.rmtree(tmp)
    bad = sum(1 for x in out if not x[3])
    print(f"qref self-check: {len(out)} cases, {bad} wrong ({time.time() - t0:.0f} s)")
    return bad, ran
