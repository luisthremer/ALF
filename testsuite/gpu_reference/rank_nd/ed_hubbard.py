#!/usr/bin/env python3
"""Ticket 401: exact finite-temperature ED of the Hubbard ring (N sites, periodic), grand canonical, for the Plain_Vanilla check.
usage: ed_hubbard.py [N=4] [U=4] [beta=2] [mu=0] [t=1]
Ensemble (ALF convention, HS decoupling of the shifted interaction):
    H = -t sum_<ij>,s (c+_is c_js + h.c.) - mu N_tot + U sum_i (n_iu - 1/2)(n_id - 1/2).
Observables as in Hamiltonian_Hubbard_Plain_Vanilla Obser (scal): Kin = <H_hop>, Pot = U <sum_i n_iu n_id> (unshifted), Part = <N_tot>,
Ener = Kin + Pot.  For N = 4, U = 4, beta = 2, mu = 0 this gives Ener = -1.47261997 (the value of test_vs_ed/test_specs.yaml).
numpy only; 4**N states (N <= 6)."""
import sys
import numpy as np

def parts(N, t):
    nm = 2 * N                      # modes: site i spin s -> i + N*s
    dim = 1 << nm
    Hhop = np.zeros((dim, dim)); D = np.zeros(dim); Nt = np.zeros(dim)
    for st in range(dim):
        occ = [(st >> m) & 1 for m in range(nm)]
        D[st] = sum(occ[i] * occ[i + N] for i in range(N)); Nt[st] = sum(occ)
        for s in range(2):
            for i in range(N):
                j = (i + 1) % N
                if N == 2 and i == 1:
                    continue
                a, b = i + N * s, j + N * s
                for (p, q) in ((a, b), (b, a)):       # c+_p c_q
                    if occ[q] == 1 and occ[p] == 0:
                        new = st ^ (1 << q)
                        sign = (-1) ** bin(st & ((1 << q) - 1)).count("1")
                        sign *= (-1) ** bin(new & ((1 << p) - 1)).count("1")
                        new ^= (1 << p)
                        Hhop[new, st] += -t * sign
    return Hhop, D, Nt

def observables(N=4, U=4.0, beta=2.0, mu=0.0, t=1.0):
    Hhop, D, Nt = parts(N, t)
    assert np.allclose(Hhop, Hhop.T)
    H = Hhop + np.diag(-mu * Nt + U * (D - 0.5 * Nt + N / 4.0))
    w, v = np.linalg.eigh(H)
    p = np.exp(-beta * (w - w.min())); p /= p.sum()
    rho_diag = lambda A: float(sum(p[k] * (v[:, k] @ (A @ v[:, k])) for k in range(len(w))))
    kin = rho_diag(Hhop)
    pot = U * rho_diag(np.diag(D))
    part = rho_diag(np.diag(Nt))
    return {"Kin_scal": kin, "Pot_scal": pot, "Part_scal": part, "Ener_scal": kin + pot}

if __name__ == "__main__":
    a = [float(x) for x in sys.argv[1:]]
    N = int(a[0]) if len(a) > 0 else 4
    U = a[1] if len(a) > 1 else 4.0
    beta = a[2] if len(a) > 2 else 2.0
    mu = a[3] if len(a) > 3 else 0.0
    t = a[4] if len(a) > 4 else 1.0
    for k, x in observables(N, U, beta, mu, t).items():
        print("%-10s %.8f" % (k, x))
