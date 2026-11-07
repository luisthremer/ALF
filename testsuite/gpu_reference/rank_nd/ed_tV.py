#!/usr/bin/env python3
"""Exact diagonalisation of the ALF tV model (Prog/Hamiltonians/Hamiltonian_tV_smod.F90), N_SUN = 2, N_FL = 1.

Operator ALF samples (Ham_Hop + Ham_Vint, lattice 'Square', Phi = 0, Symm / Checkerboard irrelevant for H):
    H = sum_{I, d} [ -t sum_s (c^+_{I s} c_{I+d, s} + h.c.) ]  - mu sum_{I,s} n_{I s}
        - (V/N_SUN) sum_{I, d} ( sum_s (c^+_{I s} c_{I+d, s} + h.c.) )^2
  d runs over a_1 (and a_2 if L2 > 1) for EVERY cell I, i.e. for L = 2 the bond I -> I+a1 exists twice
  (I=0 -> 1 and I=1 -> 0): hopping 2t and two interaction vertices.  Each vertex is Op_V with alpha = 0 and
  g = sqrt(dtau*V/N_SUN) type 2, i.e. exp(-dtau w A^2) with w = -g^2/dtau = -V/N_SUN (Obser: weight = -g**2/dtau).
Measured (Obser): Kin = <H_T> (hopping and -mu n), Pot = <H_V>, Part = <N>, Ener = Kin + Pot.

DUPLICATED BONDS (L1 = 2 or L2 = 2, found with this ED against the CPU build, see README.md):
  * Checkerboard = .true.: every (cell, direction) bond is its own two-site operator -> hopping 2t on a duplicated pair, and
    Kin (sum over the same bond list) is consistent: Kin = <H_T>.
  * Checkerboard = .false.: Predefined_Hoppings_set_OPT fills Op_T%O(I1,J1) with '=' (not '+='), so a duplicated pair hops
    with t only in the PROPAGATION, while Compute_Kin still sums over every (cell, direction) bond: Kin = 2 * <H_T,propagated>
    on a duplicated pair.  The ED mirrors this (argument checkerboard).  Pot (Op_V per bond) is doubled in both cases.
"""
import sys
import json
import numpy as np
import scipy.sparse as sp

from edcommon import Fock, thermal, lattice_sites, shift, eq_lattice_ops


def build(L1, L2, t, mu, V, N_SUN=2, checkerboard=True):
    """returns F, cells, idx, HT (propagated hopping + chemical potential), HT_kin (what Compute_Kin measures), HV"""
    assert N_SUN == 2
    cells, idx = lattice_sites(L1, L2)
    F = Fock(len(cells))
    dirs = [(1, 0)] + ([(0, 1)] if L2 > 1 else [])
    Z = sp.csr_matrix((F.dim, F.dim))
    HT, HTk, HV = Z, Z, Z
    seen = set()
    for I in cells:
        for d in dirs:
            J = shift(I, d, L1, L2)
            A = F.hopsym(idx[I], idx[J])
            HTk = HTk - t * A                                   # Compute_Kin: every (cell, direction) bond
            key = tuple(sorted((idx[I], idx[J])))
            if checkerboard or key not in seen:                 # dense Op_T: '=' overwrites a duplicated pair
                HT = HT - t * A
            seen.add(key)
            HV = HV - (V / N_SUN) * (A @ A)                     # one Op_V per (cell, direction) bond
    chem = mu * sum(F.n(i) for i in range(len(cells)))
    return F, cells, idx, (HT - chem).tocsr(), (HTk - chem).tocsr(), HV.tocsr()


def compute(L1, L2, t=1.0, mu=0.0, V=4.0, beta=2.0, N_SUN=2, eq=True, checkerboard=True):
    F, cells, idx, HT, HTk, HV = build(L1, L2, t, mu, V, N_SUN, checkerboard)
    ops = {"Kin": HTk, "Pot": HV, "Part": sum(F.n(i) for i in range(len(cells))).tocsr()}
    if eq:
        ops.update(eq_lattice_ops(F, L1, L2, 1, lambda c, a: idx[c]))
    r = thermal(HT + HV, ops, beta, F.sector)   # propagated Hamiltonian = HT + HV
    r["Ener"] = r["Kin"] + r["Pot"]
    return r


if __name__ == "__main__":
    L1, L2 = int(sys.argv[1]), int(sys.argv[2])
    kw = dict(a.split("=") for a in sys.argv[3:])
    kw = {k: (bool(int(v)) if k == "checkerboard" else float(v)) for k, v in kw.items()}
    r = compute(L1, L2, **kw)
    print(json.dumps({k: v for k, v in r.items() if not k.startswith("_") and "_eq[" not in k}, indent=1))
