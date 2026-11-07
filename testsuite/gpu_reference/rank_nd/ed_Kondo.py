#!/usr/bin/env python3
"""Exact diagonalisation of the ALF Kondo-lattice model (Prog/Hamiltonians/Hamiltonian_Kondo_smod.F90),
Lattice 'Bilayer_square', N_SUN = 2, N_FL = 1, Phi = 0.

Orbital 1 of each cell = c (conduction), orbital 2 = f (Invlist(I,2)).
Operator ALF samples (Ham_Hop, Ham_V; Uc, Uf, JK > 0 as set):
  H = -t sum_{I,d} sum_s (c^+_{I s} c_{I+d s} + h.c.)  - mu sum n_c          (only c-layer hops, T_Loc(f) = 0 !)
      + (Uc/N_SUN) sum_I ( sum_s (n^c_{I s} - 1/2) )^2  + (Uf/N_SUN) sum_I ( sum_s (n^f_{I s} - 1/2) )^2
      - (JK/(2 N_SUN)) sum_I ( sum_s (c^+_{I s} f_{I s} + h.c.) )^2
  d = a_1 for every cell (L2 = 1) or a_1, a_2 (L2 > 1); for L = 2 the bond exists twice (hopping 2t).
  JK: Predefined_Int_V_SUN(.., JK/2): g = sqrt(dtau JK/(2 N_SUN)), alpha = 0, weight -g^2/dtau = -JK/(2 N_SUN) = -JK/4.
  Uc/Uf: Predefined_Int_U_SUN: alpha = -1/2, weight U/N_SUN.   (U/N)(n-1)^2 = U (n_up-1/2)(n_dn-1/2) + U/4.
Measured by Obser (mirrored exactly, including its constant):
  Kin        = < H_T >                               (hopping and -mu n_c)
  Pot        = < Uc-term > + < -(JK/4) A^2 > + JK * ((Latt%N*Norb)/8 INTEGER division)
               (the Uf term is NOT part of Pot / Ener; the constant is JK*N_cells/4 only if N_cells*2 is a multiple of 8)
  Ener       = Kin + Pot
  Constraint = (1/N_SUN) sum_I < (n^f_I - 1)^2 >
  Part       = < N >
Duplicated bonds (L1 = 2 or L2 = 2): see ed_tV.py -- Checkerboard = .false. propagates a duplicated pair with t but Compute_Kin
  counts it twice; the argument checkerboard mirrors both cases.
"""
import sys
import json
import numpy as np
import scipy.sparse as sp

from edcommon import Fock, thermal, lattice_sites, shift, eq_lattice_ops


def build(L1, L2, t, mu, Uc, Uf, JK, N_SUN=2, checkerboard=True):
    assert N_SUN == 2
    cells, idx = lattice_sites(L1, L2)
    Ns = len(cells)
    F = Fock(2 * Ns)
    c = lambda I: idx[I]
    f = lambda I: Ns + idx[I]
    dirs = [(1, 0)] + ([(0, 1)] if L2 > 1 else [])
    Z = sp.csr_matrix((F.dim, F.dim))
    HT, HTk, HUc, HJ, HJe, HUf, Con = Z, Z, Z, Z, Z, Z, Z
    seen = set()
    for I in cells:
        for d in dirs:
            J = shift(I, d, L1, L2)
            A = F.hopsym(c(I), c(J))
            HTk = HTk - t * A
            key = tuple(sorted((c(I), c(J))))
            if checkerboard or key not in seen:
                HT = HT - t * A
            seen.add(key)
        HT = HT - mu * F.n(c(I))
        HTk = HTk - mu * F.n(c(I))
        mc = F.n(c(I)) - F.identity()
        mf = F.n(f(I)) - F.identity()
        HUc = HUc + (Uc / N_SUN) * (mc @ mc)
        HUf = HUf + (Uf / N_SUN) * (mf @ mf)
        Con = Con + (1.0 / N_SUN) * (mf @ mf)
        A = F.hopsym(c(I), f(I))
        HJ = HJ - (JK / (2.0 * N_SUN)) * (A @ A)
        # Obser's JK estimator is the Wick form of A^2 WITHOUT the pair terms sum_{s!=s'} (c+_s f_s c+_s' f_s' + f+_s c_s f+_s' c_s'):
        P = sp.csr_matrix((F.dim, F.dim))
        for s1 in (0, 1):
            s2 = 1 - s1
            cf1, cf2 = F.cdc(F.mode(c(I), s1), F.mode(f(I), s1)), F.cdc(F.mode(c(I), s2), F.mode(f(I), s2))
            P = P + cf1 @ cf2 + cf1.T @ cf2.T
        HJe = HJe - (JK / (2.0 * N_SUN)) * (A @ A - P)
    const = JK * float((Ns * 2) // 8)   # Obser: Real(Latt%N* Latt_unit%Norb/8) -> integer division
    return F, cells, idx, dict(HT=HT.tocsr(), HTk=HTk.tocsr(), HUc=HUc.tocsr(), HUf=HUf.tocsr(), HJ=HJ.tocsr(), HJe=HJe.tocsr(), Con=Con.tocsr()), const


def compute(L1, L2, t=1.0, mu=0.0, Uc=0.0, Uf=2.0, JK=2.0, beta=2.0, N_SUN=2, eq=True, checkerboard=True):
    F, cells, idx, H, const = build(L1, L2, t, mu, Uc, Uf, JK, N_SUN, checkerboard)
    Ns = len(cells)
    Htot = H["HT"] + H["HUc"] + H["HUf"] + H["HJ"]
    ops = {"Kin": H["HTk"], "Hubc": H["HUc"], "HJ": H["HJ"], "HJe": H["HJe"], "Constraint": H["Con"],
           "Part": sum(F.n(i) for i in range(2 * Ns)).tocsr()}
    if eq:
        ops.update(eq_lattice_ops(F, L1, L2, 2, lambda cell, a: idx[cell] + a * Ns))
    r = thermal(Htot, ops, beta, F.sector)
    r["Pot_exact_op"] = r["Hubc"] + r["HJ"] + const   # <Uc term> + <H_JK> + const (true operator expectation)
    r["Pot"] = r["Hubc"] + r["HJe"] + const             # what Obser measures (estimator without the pair terms)
    r["Ener"] = r["Kin"] + r["Pot"]
    r["_const"] = const
    r["_Ns"] = Ns
    return r


if __name__ == "__main__":
    L1, L2 = int(sys.argv[1]), int(sys.argv[2])
    kw = {k: (bool(int(v)) if k == "checkerboard" else float(v)) for k, v in (a.split("=") for a in sys.argv[3:])}
    r = compute(L1, L2, **kw)
    print(json.dumps({k: v for k, v in r.items() if "_eq[" not in k}, indent=1))
