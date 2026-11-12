#!/usr/bin/env python3
"""Statistical error of SpinZ_localtau at a given tau slice and site, from a live or finished data.h5.
usage: localtau_err.py RUNDIR [WARMUP_BINS] [TAU_SLICE] [SITE]   (RUNDIR contains Temp_*/data.h5)
SpinZ_localtau/obser has shape (nbin, 1, Nt, N, 2): [bin, 0, tau slice k (tau = k*dtau, k = 0..Ltrot), site, re/im].
Site 0 = lattice point (0,0). The pinning site is N/2 - L/2 (1-based) -> index N/2 - L/2 - 1 (0-based).
Mean = sum(obs*sign)/sum(sign) per bin set; error = jackknife over bins, repeated for rebinned bin sizes (1,2,4,...)
to expose autocorrelation: a growing error with the bin size means neighbouring bins are correlated."""
import sys, glob, os, numpy as np, h5py
run = sys.argv[1]; n0 = int(sys.argv[2]) if len(sys.argv) > 2 else 0
tau = int(sys.argv[3]) if len(sys.argv) > 3 else 0; site = int(sys.argv[4]) if len(sys.argv) > 4 else 0
def jack(x, s):                                   # x, s: (nbin,) ; ratio estimator sum(x s)/sum(s), jackknife error
    n = len(x); xs = x*s; tot_x, tot_s = xs.sum(), s.sum()
    est = tot_x/tot_s; loo = (tot_x - xs)/(tot_s - s)
    return est, np.sqrt((n-1)/n*((loo - loo.mean())**2).sum())
for d in sorted(glob.glob(os.path.join(run, "Temp_*"))):
    with h5py.File(os.path.join(d, "data.h5"), "r") as f:
        g = f["SpinZ_localtau"]; Nt = g["obser"].shape[2]; N = g["obser"].shape[3]; L = int(round(N**0.5))
        pin = N//2 - L//2 - 1
        sg = np.array(g["sign"]).ravel()
        U = [l for l in open(os.path.join(d, "parameters")) if l.strip().startswith("ham_U ")][0].split("=")[1].split()[0]
        print(f"== {os.path.basename(d)}  U={U}  bins={len(sg)} (discard first {n0}), Nt={Nt}, tau slice {tau}")
        for name, st in (("site 0 (0,0)", site), ("pin site", pin)):
            x = np.array(g["obser"][:, 0, tau, st, 0])[n0:]; s = sg[n0:]
            print(f"  {name} (index {st}):")
            b = 1
            while len(x)//b >= 8:
                m = len(x)//b; xb = x[:m*b].reshape(m, b).mean(1); sb = s[:m*b].reshape(m, b).mean(1)
                est, err = jack(xb, sb)
                print(f"    bin size {b:3d} ({m:3d} bins): mean {est: .6e}  err {err:.2e}  rel {abs(err/est)*100:.3f} %")
                b *= 2
