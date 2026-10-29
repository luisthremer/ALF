#!/usr/bin/env python3
"""Script-level positive controls: copy a run directory and plant one defect, so that a gate can be shown to FAIL.

usage: mutate.py RUN_IN RUN_OUT KIND
  KIND  scale:<dataset>:<rel>     multiply the smallest-magnitude nonzero entry of one bin of <dataset> by (1 + rel)
        shift:<dataset>:<nsigma>  add nsigma * (per-bin standard deviation) to every bin of <dataset>
        nan:<dataset>             put a NaN into one bin of <dataset>
        drop:<group>              delete an HDF5 group or dataset (e.g. SpinZ_tau)
        seeds                     rotate the seeds file by 37 lines (must break a same-seed or identity gate)
        taudiag_inflate:<f>       multiply each GTT dev value in TAUDIAG lines by 10**N(0, f) (variance inflation)
        taudiag_bias:<x>          multiply each GTT dev value in TAUDIAG lines by x
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
                    nz = np.nonzero(np.abs(flat) > 0)[0]
                    i = nz[np.argmin(np.abs(flat[nz]))]
                    flat[i] = flat[i] * (1 + float(arg[1]))
                    x[x.shape[0] // 2] = flat.reshape(x[x.shape[0] // 2].shape)
                elif k == "shift":
                    x = x + float(arg[1]) * x.std(axis=0, ddof=1)
                else:
                    x.reshape(x.shape[0], -1)[x.shape[0] // 2, 0] = np.nan
                d[...] = x
    elif k == "seeds":
        s = open(os.path.join(b, "seeds")).read().splitlines()
        open(os.path.join(b, "seeds"), "w").write("\n".join(s[37:] + s[:37]) + "\n")
    elif k in ("taudiag_inflate", "taudiag_bias"):
        rng = random.Random(4711)
        out = []
        for line in open(os.path.join(b, "run.log")):
            if "TAUDIAG" in line:
                t = line.split(" ")
                i = t.index("dev") + 2  # GTT value
                while t[i] == "":
                    i += 1
                v = float(t[i])
                v = v * (10 ** rng.gauss(0, float(arg[0])) if k == "taudiag_inflate" else float(arg[0]))
                t[i] = f"{v:.2E}"
                line = " ".join(t)
            out.append(line)
        open(os.path.join(b, "run.log"), "w").write("".join(out))
    else:
        sys.exit(f"unknown kind {kind}")
    after = {f: hashlib.sha256(open(os.path.join(b, f), "rb").read()).hexdigest() for f in before}
    if after == before:
        print("mutation did not change anything")
        return 1
    pp = os.path.join(b, "provenance.json")
    if os.path.isfile(pp):
        p = json.load(open(pp))
        p["mutated"] = kind
        json.dump(p, open(pp, "w"), indent=1)
    print("mutated", b, kind)
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    sys.exit(main(*sys.argv[1:]))
