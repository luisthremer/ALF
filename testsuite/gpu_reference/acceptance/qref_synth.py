#!/usr/bin/env python3
"""Synthetic QREF run directories for verdict.py --selfcheck and the dry-run tests of ticket 207-5 (no GPU, no dumps of real data).

Writes, for one run, exactly what acc_run.sh --qref + acc_qproc.sh + t207_d5_proc.sh leave in runs/<id>/: the harness provenance and
records, and under qproc/ the records of the 207-2 layout, qref.txt in the format of the offline program (ES12.4E3), qref_spec.json and
qref_proc.json, a manifest, empty dump files (the parser checks names and hashes, not their content), thinned by the real
t207_d5_qparse.py `thin`.  Everything is read back through the real parser; this module only generates text.
Tiny geometry: Beta 1.0, Dtau 0.1 (Ltrot 10), Nwrap 5 (NSTM 2), 2 flavours, N = 16, 30 passes, stride 3 -> 40 dumps, 20 with B.
"""
import hashlib
import importlib.util
import json
import math
import os
import random
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.join(HERE, "..", "..", "..", "docs", "gpu", "scripts", "nividia")
_sp = importlib.util.spec_from_file_location("qref_synth_qp", os.path.join(TOOLS, "t207_d5_qparse.py"))
QP = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(QP)

NUL = {"def": [0, 4, 100004, 1, 100001, 2, 1608, 4808, 5, 100005], "log": [0, 4, 100004, 1, 100001, 2, 5, 100005]}
NSTM, NFL, NPASS, N = 2, 2, 30, 16
CAND = "c0ffee000000" + "0" * 28
TAR = "7a" * 32
PROG = "9b" * 32


def h(b):
    return hashlib.sha256(b).hexdigest()


def hf(p):
    with open(p, "rb") as f:
        return h(f.read())


def fmt(x):
    m, e = ("%.4E" % x).split("E")
    return "%sE%+04d" % (m, int(e))


def tool_hashes():
    return {n: hf(os.path.join(TOOLS, n)) for n in sorted(os.listdir(TOOLS)) if n.startswith("t207_d5_")}


def env_tokens(path, role, stride=3, frm=1, qvar=True):
    d = {"ALF_GPU_DEBUG": "1", "ALF_GPU_TAUDIAG_QHOST": "1", "ALF_GPU_TAUDIAG_QREF": str(stride), "ALF_GPU_TAUDIAG_FROM": str(frm)}
    if qvar:
        d["ALF_GPU_TAUDIAG_QVAR"] = "8"
    if role == "ctl":
        d.update([QP.CTL_EXTRA[path]])
    return ["%s=%s" % kv for kv in d.items()]


def qref_lines(path, role, rng, stride, qs, scale):
    """(lines of qref.txt, list of the expected br)"""
    nul = NUL[path]
    names = ["dev"] + (["var8"] if role == "main" else []) + ["nvh", "cpu"] + ["n%d" % c for c in nul] + ["u1", "u2"]
    out, brs = [], []
    cnt = k = 0
    for p in range(1, NPASS + 1):
        for nst in range(1, NSTM + 1):
            for nfe in range(1, NFL + 1):
                cnt += 1
                if (cnt - 1) % stride:
                    continue
                k += 1
                br = 1 + (k % 2)
                brs.append(br)
                hasb = 1 if nst < NSTM else 0
                s = 10 ** (-12 + 0.3 * rng.gauss(0, 1))
                segs = {}
                for nm in names:
                    if nm == "n0":
                        continue
                    f = scale.get(nm, 1.0)      # a factor, or (factor of the 5 errors, factor of the 3 propagated numbers)
                    f, fp = f if isinstance(f, tuple) else (f, f)
                    e = [f * s * 10 ** (0.3 * rng.gauss(0, 1)) for _ in range(5)]
                    pr = [fp * 1e-9 * 10 ** (0.3 * rng.gauss(0, 1)) for _ in range(3)] if hasb else [-1.0] * 3
                    segs[nm] = " ".join(fmt(x) for x in e + pr)
                segs["n0"] = segs["cpu"]
                ms = []
                for nm in names:
                    ms.append("%s %s" % (nm, segs[nm]))
                if (k - 1) % qs == 0:
                    ms += ["qself " + " ".join([fmt(1e-30)] * 5), "qbr " + " ".join([fmt(1e-30)] * 5)]
                hd = "QREF file %d pass %d nst %d nfe %d br %d hbr %d n %d ver 2 hasb %d dg 0 tq %9.2f" % (k, p, nst, nfe, br, br, N, hasb, 1.5)
                out.append(hd + " | " + " | ".join(ms))
    return out, brs


def make_run(root, rid, path="def", role="main", seed=1, stride=3, qs=20, keep=39, frm=1, scale=None, qvar=None, extra_env=(), cfg="synth",
             prog=PROG, qbuild=True):
    """runs/<rid> of the harness + qproc/.  scale: member name -> factor on the error numbers (the planted defects).  Returns the run's
    declaration for acc['runs']."""
    rng = random.Random(1000 * seed + (7 if role == "ctl" else 0) + (3 if path == "log" else 0))
    d = os.path.join(root, rid)
    q = os.path.join(d, "qproc")
    os.makedirs(os.path.join(q, "qref"))
    qv = (role == "main") if qvar is None else qvar
    tok = env_tokens(path, role, stride, frm, qv) + list(extra_env)
    nul = NUL[path]
    par = ("&VAR_Lattice\nL1 = 4\n/\n&VAR_Model_Generic\nDtau = 0.1d0\nBeta = 1.0d0\n/\n&VAR_QMC\nNbin = 3\nNsweep = 10\nNwrap = 5\n/\n")
    envinfo = ", ".join(sorted(tok))
    info = ("  Sweeps : 10\n  Bins : 3\n  Nstm : %d\n  N_FL : %d\n  dtau,Ltrot_eff: 0.1 10\n"
            " GPU tau pass (207): ACTIVE (synthetic) N = %d, Nstm = %d, flavours %d, CGR2_2 QR: %s\n GPU debug env: %s\n"
            % (NSTM, NFL, N, NSTM, NFL, "pre-pivoted (default)" if path == "def" else "host ZGEQP3 (STABLOG)", envinfo))
    lines, brs = qref_lines(path, role, rng, stride, qs, scale or {})
    qvl = [8] if qv else []
    log = " GPU tau pass: DEBUG ALF_GPU_TAUDIAG_QREF=%d active (dumps to qref/, from pass %d, host CGR2_2 outputs T, device variant masks:%s)\n" % (
        stride, frm, "".join(" %d" % m for m in qvl))
    for k, br in enumerate(brs, 1):
        # pass / nst / nfe of dump k from the header of its qref line
        t = lines[k - 1].split()
        log += " TAUQREF file %d pass %s nst %s nfe %s br %d qref/q%06d.bin\n" % (k, t[4], t[6], t[8], br, k)
    log += " synthetic run, seed %d, role %s, path %s\n" % (seed, role, path)       # real run.logs differ between seeds
    seeds = "%d %d\n" % (seed, seed + 1)
    for dd in (d, q):
        os.makedirs(dd, exist_ok=True)
        for nm, tx in (("parameters", par), ("info", info), ("run.log", log)):
            open(os.path.join(dd, nm), "w").write(tx)
    open(os.path.join(d, "seeds"), "w").write(seeds)
    ph, sh = hf(os.path.join(d, "parameters")), hf(os.path.join(d, "seeds"))
    ih, lh = hf(os.path.join(d, "info")), hf(os.path.join(d, "run.log"))
    prov = {"schema": "alf-acc-prov/1", "run_dir": d, "binary_sha256": "b" * 64, "exit_code": 0, "timed_out": False,
            "build": {"name": "gpu", "configure": "PGI noMPI HDF5 GPU=nvidia", "ALF_FLAGS_EXT": "", "binary_sha256": "b" * 64, "source": "src",
                      "tar_sha256_built": TAR, "candidate": {"sha": CAND, "tar_sha256": TAR}},
            "inputs": {"parameters_sha256": ph, "seeds_sha256": sh}, "extra_env": list(tok), "allow_debug": True,
            "debug_env_present": [t.split("=")[0] for t in tok], "env": {"OMP_TARGET_OFFLOAD": "MANDATORY"},
            "gpu_check_rc": 0, "load_ok": True, "gpu_before": [], "outputs": {"info": ih, "data.h5": None, "run.log": lh}}
    json.dump(prov, open(os.path.join(d, "provenance.json"), "w"), indent=1)
    alf = [t for t in tok if t.startswith("ALF_GPU_")]
    open(os.path.join(q, "stamp.txt"), "w").write("binary %s\nenv %s\nseeds %s parameters %s\n" % ("b" * 64, " ".join(alf), sh[:16], ph[:16]))
    json.dump({"rc": 0, "env": alf, "seeds_sha256": sh, "parameters_sha256": ph, "run_log_sha256": lh, "info_sha256": ih,
               "binary_sha256": "b" * 64}, open(os.path.join(q, "provenance.json"), "w"), indent=1)
    open(os.path.join(q, "qref.txt"), "w").write("\n".join(lines) + "\nQREFEND dumps %d ok %d err 0\n" % (len(lines), len(lines)))
    for k in range(1, len(lines) + 1):
        open(os.path.join(q, "qref", "q%06d.bin" % k), "wb").close()
    names = sorted("qref/" + f for f in os.listdir(os.path.join(q, "qref")))
    open(os.path.join(q, "qref_manifest.txt"), "w").write("".join("%s  %s\n" % (hf(os.path.join(q, n)), n) for n in names))
    sp = QP.make_spec(q, path, qs, nul, 0, 1, keep)
    json.dump(sp, open(os.path.join(q, "qref_spec.json"), "w"), indent=1)
    json.dump({"prog": prog, "args": {"path": path, "qs": qs, "nul": ",".join(str(c) for c in nul), "hostin": 0, "cover": 1, "keep": keep},
               "rc": 0, "check_rc": 0, "check": "QPARSE check ok (synthetic)"}, open(os.path.join(q, "qref_proc.json"), "w"), indent=1)
    if keep > 0:
        QP.thin(q)
    json.dump({"schema": "alf-acc-qproc/1", "args": {"path": path, "qs": qs, "nul": ",".join(str(c) for c in nul), "hostin": 0, "keep": keep},
               "proc_rc": 0, "qbuild_sha256": None, "tools_sha256": tool_hashes()}, open(os.path.join(q, "acc_qproc.json"), "w"), indent=1)
    return {"build": "gpu", "config": cfg, "env": list(tok), "allow_debug": True, "active_lines": ["GPU tau pass (207): ACTIVE"],
            "qref": {"path": path, "qs": qs, "nul": ",".join(str(c) for c in nul), "hostin": 0, "keep": keep}}


def make_qbuild(bindir, path, prog=PROG):
    os.makedirs(bindir, exist_ok=True)
    t = tool_hashes()
    rec = {"schema": "alf-acc-qbuild/1", "path": path, "prog_sha256": prog, "tar_sha256_built": TAR,
           "candidate": {"sha": CAND, "tar_sha256": TAR}, "tools_sha256": t, "qref_src_sha256": t["t207_d5_qref.F90"]}
    p = os.path.join(bindir, "qref_%s.qbuild.json" % path)
    json.dump(rec, open(p, "w"), indent=1)
    return hf(p)


def seal(root, rid, qbuild_sha):
    """record the sha256 of the qbuild record in acc_qproc.json (as acc_qproc.sh does)"""
    p = os.path.join(root, rid, "qproc", "acc_qproc.json")
    r = json.load(open(p))
    r["qbuild_sha256"] = qbuild_sha
    json.dump(r, open(p, "w"), indent=1)
