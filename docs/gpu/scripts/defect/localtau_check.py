#!/usr/bin/env python3
"""Physics sanity of SpinZ_localtau: symmetry f(tau) = f(beta - tau) and position of the maximum.
Grid (JAKOB_BACKUP 5563aecd/7e40b650): Nt = Ltrot + 1 columns, column k = slice k = 0..Ltrot (tau = k*dtau, 0..beta);
the pin is on g_t(1), i.e. at tau = 0 (= beta), so the response peaks at the endpoints k = 0 and k = Ltrot.
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
B = np.array(bins); nb = B.shape[0]; Nt = B.shape[2]; Ltrot = Nt - 1
pin = N//2 - L//2 - 1                            # pin_site = N/2 - L/2 (1-based)
mean = B.mean(0); err = B.std(0, ddof=1)/np.sqrt(nb)
f = mean[pin]; e = err[pin]
print(f"bins {nb}, Nt {Nt} (Ltrot {Ltrot}), pin site index {pin+1}")
print("max over sites and slices at (site, slice k):", np.unravel_index(np.abs(mean).argmax(), mean.shape), "pin site:", pin)
print("endpoint slices (k, value, err):", [(k, round(f[k], 4), round(e[k], 4)) for k in list(range(0, 4)) + list(range(Ltrot-3, Ltrot+1))])
zz = []
for k in range(0, Ltrot//2):                     # partner of k is Ltrot - k
    dd = B[:, pin, k] - B[:, pin, Ltrot-k]
    zz.append(dd.mean()/(dd.std(ddof=1)/np.sqrt(nb)))
zz = np.array(zz)
print("pin-site profile: tau/beta, value, err")
for k in range(0, Ltrot+1, max(1, Ltrot//12)): print(f"  {k/Ltrot:5.3f} {f[k]: .5f} {e[k]:.5f}")
print(f"symmetry f(tau) - f(beta-tau): max |z| = {np.abs(zz).max():.2f}, rms z = {np.sqrt((zz**2).mean()):.2f} over {len(zz)} pairs")
print("argmax of |profile| on pin site: slice", int(np.abs(f).argmax()), "(expected 0 or", Ltrot, ")")
