#!/usr/bin/env python3
"""Script-level positive controls: copy a run directory and plant one defect, so that a gate can be shown to FAIL.

usage: mutate.py RUN_IN RUN_OUT KIND
  KIND  scale:<dataset>:<rel>     multiply one entry of the middle bin of <dataset> by (1 + rel): the smallest entry
                                  whose magnitude is >= 1e-5 of the dataset maximum (above the per-entry floor)
        shift:<dataset>:<k>       add k * (error of the mean from 8-fold rebinned bins, i.e. robust against
                                  autocorrelation) to every bin of <dataset>
        nan:<dataset>             put a NaN into one bin of <dataset>
        drop:<group>              delete an HDF5 group or dataset (e.g. SpinZ_tau)
        seeds                     rotate the seeds file by 37 lines (must break a same-seed or identity gate)
        taudiag_inflate:<f>       multiply each GTT dev value in TAUDIAG lines by 10**N(0, f) (variance inflation)
        taudiag_bias:<x>          multiply each GTT dev value in TAUDIAG lines by x
        qref_swap:dev:var8        QREF run (207-5): device output and stand-in var8 exchanged in the kept dumps and in qproc/qref.txt
                                  (manifest and qref_proc.json follow; the driver stamps the copy in the attempt log)
The copy keeps provenance.json (marked "mutated") so the gate's provenance checks still run. RUN_OUT must not exist.
Exit 0 when the mutation was applied and verified to change the copy, 1 otherwise.
"""
import hashlib
import json
import os
import random
import shutil
import sys

import h5py
import numpy as np


def qref_swap(q):
    """QREF run (ticket 207-5): the device output and the stand-in var8 exchanged consistently, so that the processed run is what the pinned
    offline program computes from the exchanged dumps: in every kept dump the dev block and the var8 block (64 N^2 bytes each, layout of
    t207_d5_qref.F90 / Prog/gpu_tau_mod.F90 qref_post: header, QVAR masks, two states, dev, var<m>..., host, ...) are swapped; in qref.txt the
    numbers of the members dev and var8 are swapped in every line; the manifest lines of the kept dumps and qref_sha256 / manifest_sha256 of
    qref_proc.json follow.  The gate then judges var8's numbers as the device (planted defect) and dev's as the stand-in control."""
    import struct
    qd = os.path.join(q, "qref")
    kept = sorted(f for f in os.listdir(qd) if f.endswith(".bin"))
    if not kept:
        sys.exit("no kept dump to swap")
    for f in kept:
        p = os.path.join(qd, f)
        if os.path.islink(p):
            sys.exit(f"{p} is a symbolic link: refusing")
        with open(p, "r+b") as fh:
            h = struct.unpack("<13i", fh.read(52))
            if h[0] != 20791 or h[1] != 2:
                sys.exit(f"{p}: not a version 2 dump")
            n, logp, nv = h[7], h[8], h[9]
            masks = list(struct.unpack(f"<{nv}i", fh.read(4 * nv))) if nv else []
            if masks[:1] != [8]:
                sys.exit(f"{p}: the first QVAR mask is {masks[:1]}, not 8")
            sb = 32 * n * n + (8 if logp else 16) * n
            dev = 52 + 4 * nv + 2 * sb
            fh.seek(dev)
            a = fh.read(64 * n * n)
            v = fh.read(64 * n * n)
            if len(a) != 64 * n * n or len(v) != 64 * n * n or a == v:
                sys.exit(f"{p}: short dump or identical dev/var8 blocks")
            fh.seek(dev)
            fh.write(v + a)
    qt = os.path.join(q, "qref.txt")
    out = []
    for ln in open(qt).read().split("\n")[:-1]:
        if ln.startswith("QREF file "):
            segs = ln.split(" | ")
            i = [j for j, s in enumerate(segs) if s.split()[:1] == ["dev"]]
            k = [j for j, s in enumerate(segs) if s.split()[:1] == ["var8"]]
            if len(i) != 1 or len(k) != 1:
                sys.exit(f"qref.txt line without exactly one dev and one var8: {ln[:80]}")
            i, k = i[0], k[0]
            nd, nv8 = segs[i][len("dev"):], segs[k][len("var8"):]
            segs[i], segs[k] = "dev" + nv8, "var8" + nd
            ln = " | ".join(segs)
        out.append(ln)
    open(qt, "w").write("\n".join(out) + "\n")
    mp = os.path.join(q, "qref_manifest.txt")
    ml = []
    for ln in open(mp).read().split("\n")[:-1]:
        hsh, name = ln.split("  ", 1)
        if os.path.basename(name) in kept:
            hsh = hashlib.sha256(open(os.path.join(q, name), "rb").read()).hexdigest()
        ml.append(f"{hsh}  {name}")
    open(mp, "w").write("\n".join(ml) + "\n")
    pj = os.path.join(q, "qref_proc.json")
    pr = json.load(open(pj))
    pr["qref_sha256"] = hashlib.sha256(open(qt, "rb").read()).hexdigest()
    pr["manifest_sha256"] = hashlib.sha256(open(mp, "rb").read()).hexdigest()
    pr["mutated"] = "qref_swap:dev:var8 (mutate.py)"
    json.dump(pr, open(pj, "w"), indent=1)


def main(a, b, kind):
    if os.path.exists(b):
        sys.exit(f"refusing: {b} exists")
    shutil.copytree(a, b)
    k, *arg = kind.split(":")
    before = {f: hashlib.sha256(open(os.path.join(b, f), "rb").read()).hexdigest()
              for f in ("data.h5", "run.log", "seeds") if os.path.isfile(os.path.join(b, f))}
    if k in ("scale", "shift", "nan", "drop"):
        with h5py.File(os.path.join(b, "data.h5"), "r+") as f:
            if k == "drop":
                del f[arg[0]]
            else:
                d = f[arg[0]]
                x = d[()]
                if k == "scale":
                    flat = x[x.shape[0] // 2].reshape(-1)
                    nz = np.nonzero(np.abs(flat) >= 1e-5 * np.abs(x).max())[0]
                    if nz.size == 0:
                        sys.exit("no entry above the floor")
                    i = nz[np.argmin(np.abs(flat[nz]))]
                    flat[i] = flat[i] * (1 + float(arg[1]))
                    x[x.shape[0] // 2] = flat.reshape(x[x.shape[0] // 2].shape)
                elif k == "shift":
                    nb = (x.shape[0] // 8) * 8
                    if nb < 16:
                        sys.exit("too few bins for an 8-fold rebinned error")
                    xb = x[x.shape[0] - nb:].reshape((nb // 8, 8) + x.shape[1:]).mean(axis=1)
                    x = x + float(arg[1]) * xb.std(axis=0, ddof=1) / np.sqrt(xb.shape[0])
                else:
                    x.reshape(x.shape[0], -1)[x.shape[0] // 2, 0] = np.nan
                d[...] = x
    elif k == "seeds":
        s = open(os.path.join(b, "seeds")).read().splitlines()
        open(os.path.join(b, "seeds"), "w").write("\n".join(s[37:] + s[:37]) + "\n")
    elif k == "qref_swap":
        if arg != ["dev", "var8"]:
            sys.exit("qref_swap supports only dev:var8")
        before["qproc/qref.txt"] = hashlib.sha256(open(os.path.join(b, "qproc", "qref.txt"), "rb").read()).hexdigest()
        qref_swap(os.path.join(b, "qproc"))
    elif k in ("taudiag_inflate", "taudiag_bias"):
        rng = random.Random(4711)
        out = []
        for line in open(os.path.join(b, "run.log")):
            if line.lstrip().startswith("TAUDIAG pt "):  # data lines only, as verdict.py
                t = line.split()  # same tokenisation as verdict.py
                i = t.index("dev") + 2  # GTT value
                v = float(t[i])
                v = v * (10 ** rng.gauss(0, float(arg[0])) if k == "taudiag_inflate" else float(arg[0]))
                t[i] = f"{v:.2E}"
                line = " " + " ".join(t) + "\n"
            out.append(line)
        open(os.path.join(b, "run.log"), "w").write("".join(out))
    else:
        sys.exit(f"unknown kind {kind}")
    after = {f: hashlib.sha256(open(os.path.join(b, f), "rb").read()).hexdigest() for f in before}
    if after == before:
        print("mutation did not change anything")
        return 1
    pp = os.path.join(b, "provenance.json")
    if not os.path.isfile(pp):
        print("no provenance.json in the source run")
        return 1
    p = json.load(open(pp))
    p["mutated"] = kind
    p["outputs_before_mutation"] = p.get("outputs", {})
    p["outputs"] = {f: hashlib.sha256(open(os.path.join(b, f), "rb").read()).hexdigest() if os.path.isfile(os.path.join(b, f)) else None
                    for f in ("info", "data.h5", "run.log")}
    if k == "seeds":
        p.setdefault("inputs", {})["seeds_sha256"] = hashlib.sha256(open(os.path.join(b, "seeds"), "rb").read()).hexdigest()
    json.dump(p, open(pp, "w"), indent=1)
    print("mutated", b, kind)
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    sys.exit(main(*sys.argv[1:]))
