#!/usr/bin/env python3
"""Exact diagonalisation of the ALF Spin-Peierls model (Prog/Hamiltonians/Hamiltonian_Spin_Peierls_smod.F90),
SU(2) symmetric branch (Ham_h = 0, N_SUN = 2, N_FL = 1), Lattice 'Square' with L2 = 1 (the smallest the code allows: L1 = 2).

What ALF samples (Ham_Hop: Op_T with g = 0, i.e. NO fermion hopping; Ham_V; S0):
  per site I    :  type-2 vertex, alpha = -1/2, weight U/N_SUN:       H_U = (U/N_SUN) (sum_s (n_s - 1/2))^2 = (U/2)(n_I - 1)^2
  per bond b    :  type-4 vertex  exp( phi g K_b ),  K_b = sum_s (c^+_{I s} c_{J s} + h.c.),  g = sqrt(dtau J/4),
                   phi = eta(l) sqrt(1 + y_b)   =>  fermionic operator  exp(-dtau * (-(J/4)(1 + y_b)) K_b^2)  i.e. H_b = -(J/4)(1+Q_b) K_b^2
  S0(y)         :  M (y(t+1) - y(t))^2 /(2 dtau) + k dtau (y + J/(4k))^2 / 2   per time slice and bond,
                   with k = 1/(2 Lambda), M = k/Omega0^2 (Ham_k, Ham_M)  =>  oscillator  P^2/(2M) + k Q^2/2 + (J/4) Q  (+const)
  so   H = sum_I (U/2)(n_I-1)^2 + sum_b [ P_b^2/(2M) + k Q_b^2/2 + (J/4) Q_b - (J/4)(1 + Q_b) K_b^2 ]   (constants dropped)
  (for U -> inf: K^2 = 1 - 4 S_I.S_J on singly occupied states, H_b = P^2/2M + kQ^2/2 + J(1+Q) S_I.S_J + const: the Spin-Peierls
   model with spin-phonon coupling J Q S.S.)
  Bonds: for every cell I one x-bond I -> I+a1 (and, for L2 > 1, one y-bond); for L1 = 2 the bond exists twice (two phonons).
  Each bond has its own oscillator, y_b = Q_b = (a + a^+)/sqrt(2 M omega), omega = Omega0 (harmonic basis, truncated at nmax).
ALF writes no energy for this model; the observables are (Obser):
  Part   = <N>,   Pot = sum_I <n_up n_dn> (= sum_I GRC_II^2: double occupancy, NOT an energy),
  PhiXY  = ( (1/Ns) sum_I <y_{x-bond(I)}> , (1/Ns) sum_I <y_{y-bond(I)}> ),
  SpinZ_eq, Phi_eq (<y y>, Latt0 = <y>), Dimer_eq (<D D>, Latt0 = <D>, D_{i m} = sum_{ab} S^a_b(i) S^b_a(m),
  S^a_b(i) = c^+_a c_b - delta_ab/2),  all per site/cell normalised, analysis subtracts Latt0 products for Phi and Dimer.
Equal-time y observables of the time-discretised path integral equal the quantum expectation values up to Trotter error.
"""
import sys
import json
import numpy as np
import scipy.sparse as sp

from edcommon import Fock, thermal, lattice_sites, shift, latt_displacements


def boson_ops(nmax):
    a = sp.diags(np.sqrt(np.arange(1, nmax)), 1, format="csr")
    return a, a.T.tocsr()


def build(L1, L2, J, U, lam, om0, nmax):
    cells, idx = lattice_sites(L1, L2)
    Ns = len(cells)
    F = Fock(Ns)
    dirs = [(1, 0)] + ([(0, 1)] if L2 > 1 else [])
    bonds = [(I, d) for d in dirs for I in cells]       # one phonon per (cell, direction)
    nb = len(bonds)
    k = 1.0 / (2.0 * lam)
    M = k / om0 ** 2
    om = om0
    a, ad = boson_ops(nmax)
    x = (a + ad) / np.sqrt(2.0 * M * om)            # Q
    nn = (ad @ a)
    Ib = sp.identity(nmax, format="csr")

    def bop(m, ib):
        """operator m acting on oscillator ib, identity on the others, as the LAST factors of F (x) b_1 (x) ... (x) b_nb"""
        out = None
        for j in range(nb):
            f = m if j == ib else Ib
            out = f if out is None else sp.kron(out, f, format="csr")
        return out

    dimF, dimB = F.dim, nmax ** nb
    IF = sp.identity(dimF, format="csr")
    IB = sp.identity(dimB, format="csr")
    kronF = lambda A: sp.kron(A, IB, format="csr")       # fermion operator on full space
    kronB = lambda B: sp.kron(IF, B, format="csr")       # boson operator on full space
    H = sp.csr_matrix((dimF * dimB, dimF * dimB))
    Qs = []
    for ib, (I, d) in enumerate(bonds):
        J_ = shift(I, d, L1, L2)
        K = F.hopsym(idx[I], idx[J_])
        K2 = kronF(K @ K)
        Q = kronB(bop(x, ib))
        Qs.append(Q)
        H = H + kronB(bop(om * (nn + 0.5 * Ib), ib)) + (J / 4.0) * Q - (J / 4.0) * (K2 + Q @ K2)
    for I in cells:
        m = F.n(idx[I]) - F.identity()
        H = H + kronF((U / 2.0) * (m @ m))
    sector = np.repeat(F.sector, dimB)
    return dict(F=F, cells=cells, idx=idx, bonds=bonds, H=H.tocsr(), Qs=Qs, kronF=kronF, sector=sector, dimB=dimB, Ns=Ns)


def dimer_op(F, i, m):
    """D_{i m} = sum_{a,b} S^a_b(i) S^b_a(m),  S^a_b(i) = c^+_{i a} c_{i b} - delta_{ab}/2"""
    D = sp.csr_matrix((F.dim, F.dim))
    for al in (0, 1):
        for be in (0, 1):
            Si = F.cdc(F.mode(i, al), F.mode(i, be)) - (0.5 * F.identity() if al == be else 0 * F.identity())
            Sm = F.cdc(F.mode(m, be), F.mode(m, al)) - (0.5 * F.identity() if al == be else 0 * F.identity())
            D = D + Si @ Sm
    return D.tocsr()


def compute(L1, L2=1, J=1.0, U=4.0, lam=0.1, om0=2.0, beta=2.0, nmax=12, eq=True):
    S = build(L1, L2, J, U, lam, om0, nmax)
    F, cells, idx, bonds, Ns = S["F"], S["cells"], S["idx"], S["bonds"], S["Ns"]
    kronF = S["kronF"]
    nbond = len(bonds)
    ops = {}
    ops["Part"] = kronF(sum(F.n(i) for i in range(Ns)).tocsr())
    ops["Pot"] = kronF(sum(F.nspin(i, 0) @ F.nspin(i, 1) for i in range(Ns)).tocsr())
    # Phi: bond (I, d=a1) -> no = 1, (I, d=a2) -> no = 2
    nco = 1 if L2 == 1 else 2
    bidx = {(I, no): S["bonds"].index((I, (1, 0) if no == 0 else (0, 1))) for I in cells for no in range(nco)}
    ops["PhiXY[0]"] = sum(S["Qs"][bidx[(I, 0)]] for I in cells) / Ns
    ops["PhiXY[1]"] = (sum(S["Qs"][bidx[(I, 1)]] for I in cells) / Ns) if nco == 2 else 0 * ops["PhiXY[0]"]
    if eq:
        for (dx, dy) in latt_displacements(L1, L2):
            acc_s = None
            acc_p = {}
            acc_d = {}
            for I in cells:
                J_ = ((I[0] - dx) % L1, (I[1] - dy) % L2)
                t = (F.sz2(idx[I]) @ F.sz2(idx[J_])).tocsr()
                acc_s = t if acc_s is None else acc_s + t
            ops[f"SpinZ_eq[{dx},{dy},0,0]"] = kronF(acc_s / Ns)
            for no in range(nco):
                for no2 in range(nco):
                    sp_ = None
                    sd_ = None
                    for I in cells:
                        J_ = ((I[0] - dx) % L1, (I[1] - dy) % L2)
                        t = S["Qs"][bidx[(I, no)]] @ S["Qs"][bidx[(J_, no2)]]
                        sp_ = t if sp_ is None else sp_ + t
                        I1 = shift(I, (1, 0) if no == 0 else (0, 1), L1, L2)
                        J1 = shift(J_, (1, 0) if no2 == 0 else (0, 1), L1, L2)
                        d = (dimer_op(F, idx[I], idx[I1]) @ dimer_op(F, idx[J_], idx[J1])).tocsr()
                        sd_ = d if sd_ is None else sd_ + d
                    ops[f"Phi_eq[{dx},{dy},{no},{no2}]"] = (sp_ / Ns).tocsr()
                    ops[f"Dimer_eq[{dx},{dy},{no},{no2}]"] = kronF(sd_ / Ns)
        for no in range(nco):
            ops[f"Phi_eq0[{no}]"] = sum(S["Qs"][bidx[(I, no)]] for I in cells) / Ns
            ops[f"Dimer_eq0[{no}]"] = kronF(sum(dimer_op(F, idx[I], idx[shift(I, (1, 0) if no == 0 else (0, 1), L1, L2)]) for I in cells) / Ns)
    r = thermal(S["H"], ops, beta, S["sector"])
    r["_nmax"] = nmax
    r["_dim"] = S["H"].shape[0]
    return r


def alf_arrays(r, L1, L2):
    """ALF/pyALF-style arrays [a][b][iR] in Lattice%list order; Phi and Dimer connected (Latt0 products subtracted)"""
    from edcommon import alf_lattice_listr
    listr = alf_lattice_listr(L1, L2)
    nco = 1 if L2 == 1 else 2
    out = {}
    for k, bg in (("SpinZ", None), ("Phi", "Phi_eq0"), ("Dimer", "Dimer_eq0")):
        arr = np.zeros((nco, nco, len(listr)))
        for iR, (i1, i2) in enumerate(listr):
            for a in range(nco):
                for b in range(nco):
                    key = f"{k}_eq[{i1 % L1},{i2 % L2},{a},{b}]"
                    v = r[key]
                    if bg:
                        v -= r[f"{bg}[{a}]"] * r[f"{bg}[{b}]"]
                    arr[a, b, iR] = v
        out[f"{k}_eqR"] = arr
    return out, listr


if __name__ == "__main__":
    L1 = int(sys.argv[1])
    kw = {k: (int(v) if k == "nmax" else float(v)) for k, v in (a.split("=") for a in sys.argv[2:])}
    r = compute(L1, 1, **kw)
    print(json.dumps({k: v for k, v in r.items() if "_eq[" not in k}, indent=1))
