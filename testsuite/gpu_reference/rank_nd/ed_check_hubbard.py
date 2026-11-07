#!/usr/bin/env python3
"""Sanity check of the shared conventions: 4-site Hubbard ring, U=4, beta=2 -> upstream value -1.47261997
(testsuite/test_vs_ed/test_specs.yaml hubbard_finT_chain)."""
import scipy.sparse as sp
from edcommon import Fock, thermal

L, t, U, beta = 4, 1.0, 4.0, 2.0
F = Fock(L)
H = sp.csr_matrix((F.dim, F.dim))
for i in range(L):
    H = H - t * F.hopsym(i, (i + 1) % L)
    # ALF Predefined_Int_U_SUN: alpha = -1/2, weight U/N_SUN:  (U/N_SUN) (sum_s (n_s - 1/2))^2 = (U/2)(n-1)^2
    #   = U (n_up-1/2)(n_dn-1/2) + U/4   -> the energy ALF/upstream quote contains the +U/4 per site
    m = F.n(i) - F.identity()
    H = H + (U / 2.0) * (m @ m)
e = thermal(H, {"E": H}, beta, F.sector)["E"]
print(f"E = {e:.8f}   (reference -1.47261997, diff {e + 1.47261997:.2e})")
assert abs(e + 1.47261997) < 2e-8
