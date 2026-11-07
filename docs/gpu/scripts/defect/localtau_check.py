#!/usr/bin/env python3
"""Physics sanity of SpinZ_localtau: symmetry about beta/2 and position of the maximum.
usage: localtau_check.py RUNDIR L   (reads RUNDIR/SpinZ_localtau, all bins; error = std of bin means / sqrt(nbin))"""
import sys, re, numpy as np
d, L = sys.argv[1], int(sys.argv[2]); N = L*L
lines = open(d + "/SpinZ_localtau").read().split("\n")
bins = []; i = 0
while i < len(lines) and lines[i].strip():
    h = lines[i].split(); Nt = int(h[3]); i += 1
    arr = np.zeros((N, Nt))
    for s in range(N):
        i += 1                                   # coordinate line
        for t in range(Nt):
            arr[s, t] = float(re.findall(r"[-+]?\d\.\d+E[-+]\d+", lines[i])[0]); i += 1
    bins.append(arr)
B = np.array(bins); nb = B.shape[0]; Ltrot = B.shape[2]
pin = N//2 - L//2 - 1                            # pin_site = N/2 - L/2 (1-based)
mean = B.mean(0); err = B.std(0, ddof=1)/np.sqrt(nb)
f = mean[pin]; e = err[pin]
print(f"bins {nb}, Nt {Ltrot}, pin site index {pin+1}")
print("max over sites and slices at (site, slice):", np.unravel_index(np.abs(mean).argmax(), mean.shape), "pin site:", pin, "slice Ltrot/2 =", Ltrot//2)
# slice index nt = 1..Ltrot  (python index nt-1); symmetric partner of nt is Ltrot - nt for 1 <= nt <= Ltrot-1
print("central slices (slice, value, err):", [(nt, round(f[nt-1], 4), round(e[nt-1], 4)) for nt in range(Ltrot//2-3, Ltrot//2+5)])
def sym_z(center2):          # partner of slice nt is center2 - nt (center2 = 2 * centre slice)
    zz = []
    for nt in range(1, center2//2 + 1):
        p = center2 - nt
        if p < 1 or p > Ltrot or p == nt: continue
        dd = B[:, pin, nt-1] - B[:, pin, p-1]
        zz.append(dd.mean()/(dd.std(ddof=1)/np.sqrt(nb)))
    return np.array(zz)
for c2, lab in ((Ltrot, "slice Ltrot/2 (partner nt <-> Ltrot-nt)"), (Ltrot-1, "Ltrot/2-1/2 (partner nt <-> Ltrot-1-nt)"), (Ltrot+1, "Ltrot/2+1/2 (partner nt <-> Ltrot+1-nt)")):
    zz = sym_z(c2); print(f"symmetry about {lab}: max|z| = {np.abs(zz).max():.2f}, rms z = {np.sqrt((zz**2).mean()):.2f} over {len(zz)} pairs")
z = sym_z(Ltrot)
print("pin-site profile: tau/beta, value, err (every 5th slice)")
for nt in range(1, Ltrot+1, max(1, Ltrot//12)): print(f"  {nt/Ltrot:5.3f} {f[nt-1]: .5f} {e[nt-1]:.5f}")
print(f"symmetry f(t) - f(beta-t): max |z| = {np.abs(z).max():.2f}, rms z = {np.sqrt((z**2).mean()):.2f} over {len(z)} pairs")
print("argmax of |profile| on pin site: slice", int(np.abs(f).argmax())+1, "(expected", Ltrot//2, "or", Ltrot//2+1, ")")
