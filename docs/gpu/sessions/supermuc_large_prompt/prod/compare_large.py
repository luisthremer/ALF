#!/usr/bin/env python3
"""Same-seed GPU vs CPU acceptance for L >= 16 (owner decision 2026-11-11, see docs/gpu/DECISIONS.md).

usage: compare_large.py GPU_DIR CPU_DIR CPU2_DIR [--factor 10] [--plant FILE:INDEX:REL]
  CPU2_DIR = the same CPU run with another MKL thread count (different summation order): its deviation from CPU_DIR is
  the rounding spread of correct code. Criterion per observable file:
    (1) confout_* bit-identical GPU vs CPU, acceptance identical;
    (2) real parts: element-wise relative deviation (entries with |x| > 1e-3 x file max) <= factor x CPU-vs-CPU spread,
        and absolute deviation over all entries (small ones too) <= factor x CPU-vs-CPU absolute spread;
    (3) imaginary parts: absolute deviation <= factor x CPU-vs-CPU absolute spread of the imaginary parts.
  Imaginary parts are judged in absolute terms: in the scalar observables they are small sums with heavy cancellation and
  their value depends on the summation order (GPU reductions sum in another order than the CPU).
  A spread of exactly 0 is replaced by a floor of 1e-15 (relative) / 1e-15 x file max (absolute).
  --plant FILE:INDEX:REL adds REL x (file max) to real entry INDEX of FILE in the GPU data: negative control.
exit 0 = PASS, 1 = FAIL."""
import sys, os, re, glob, filecmp
cplx = re.compile(r"\(\s*([-+0-9.EeDd]+)\s*,\s*([-+0-9.EeDd]+)\s*\)")
num = re.compile(r"[-+]?\d+\.\d*(?:[EeDd][-+]?\d+)?")

def parse(path):
    txt = open(path).read()
    re_, im = [], []
    for a, b in cplx.findall(txt):
        re_.append(float(a.replace("D", "E"))); im.append(float(b.replace("D", "E")))
    rest = cplx.sub(" ", txt)
    re_ += [float(x.replace("D", "E")) for x in num.findall(rest)]
    return re_, im

def acc(d):
    for l in open(os.path.join(d, "info")):
        if l.strip().startswith("Acceptance  "): return l.split(":")[1].strip()
    return None

def dev(a, b):
    """(element-wise relative over significant real entries, absolute over imaginary parts, file max)"""
    (ra, ia), (rb, ib) = a, b
    if len(ra) != len(rb) or len(ia) != len(ib): return None
    fmax = max([abs(x) for x in rb] + [abs(x) for x in ib] + [0.0]) or 1.0
    rel = max([abs(x - y) / abs(y) for x, y in zip(ra, rb) if abs(y) > 1e-3 * fmax] + [0.0])
    ab = max([abs(x - y) for x, y in zip(ia, ib)] + [0.0])
    rab = max([abs(x - y) for x, y in zip(ra, rb)] + [0.0])        # real parts, absolute, all entries (small ones too)
    return rel, ab, fmax, rab

def main():
    g, c, c2 = sys.argv[1:4]
    f = float(sys.argv[sys.argv.index("--factor") + 1]) if "--factor" in sys.argv else 10.0
    plant = sys.argv[sys.argv.index("--plant") + 1] if "--plant" in sys.argv else None
    ok = True
    for t in sorted(glob.glob(os.path.join(c, "Temp_*"))) or [c]:
        sub = os.path.relpath(t, c) if t != c else ""
        G, C, C2 = (os.path.join(x, sub) for x in (g, c, c2))
        print(f"== {sub or '.'}")
        for cf in sorted(glob.glob(os.path.join(C, "confout_*"))):
            same = filecmp.cmp(cf, os.path.join(G, os.path.basename(cf)), shallow=False)
            ok &= same
            if not same: print(f"FAIL {os.path.basename(cf)} differs")
        a1, a2 = acc(G), acc(C)
        ok &= (a1 == a2)
        print(f"acceptance GPU {a1} CPU {a2} {'ok' if a1 == a2 else 'FAIL'}")
        names = sorted(os.path.basename(p) for p in glob.glob(os.path.join(C, "*"))
                       if not p.endswith(("_info", "parameters", "seeds", "info", "RUNNING", "TREE")) and "confout" not in p)
        print(f"{'file':16s} {'re GPU-CPU':>11s} {'re spread':>10s} {'im GPU-CPU':>11s} {'im spread':>10s} verdict")
        for n in names:
            pg, pc, pc2 = parse(os.path.join(G, n)), parse(os.path.join(C, n)), parse(os.path.join(C2, n))
            if plant and plant.split(":")[0] == n:
                i, r = int(plant.split(":")[1]), float(plant.split(":")[2]); pg[0][i] += r * max(abs(x) for x in pc[0] + pc[1])
            d, s = dev(pg, pc), dev(pc2, pc)
            if d is None or s is None:
                ok = False; print(f"{n:16s} SHAPE MISMATCH"); continue
            lim_re = f * max(s[0], 1e-15); lim_im = f * max(s[1], 1e-15 * s[2])
            # imaginary parts of the scalar observables (*_scal): report only (decision 43). The GPU sums them with device
            # reductions (gpu_obs_mod.F90 ~284-303), the CPU with host loops; the CPU-vs-CPU control (MKL threads) does not
            # vary that summation order, so its spread is no yardstick for these cancellation sums.
            im_gate = not n.endswith("_scal")
            lim_rab = f * max(s[3], 1e-15 * s[2])
            good = d[0] <= lim_re and d[3] <= lim_rab and (d[1] <= lim_im or not im_gate)
            ok &= good
            note = "" if im_gate else f"  (im report-only, {d[1] / max(s[1], 1e-300):.0f}x spread)"
            print(f"{n:16s} {d[0]:11.2e} {s[0]:10.2e} {d[1]:11.2e} {s[1]:10.2e} {'PASS' if good else 'FAIL'}{note}")
    print("RESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

main()
