"""Shared tools for the independent exact-diagonalisation references (rank_nd/ed_*.py).

Only numpy/scipy.  Conventions (all derived from Prog/ sources, see README.md):
  * modes: m = orb + Norb*spin, spin 0 = up, 1 = down; Fock basis of 2*Norb modes, bit m of the state integer.
  * Jordan-Wigner order = mode index.
  * H is the operator ALF samples (exp(-beta H), grand canonical, Z = Tr exp(-beta H)); constant offsets that
    the Fortran Obser adds to the *measured* energies are mirrored in the individual ed_*.py files.
  * thermal() diagonalises block by block in the conserved (N_up, N_dn) sectors of the fermions.
"""
import numpy as np
import scipy.sparse as sp


def popcount(a):
    a = np.asarray(a, dtype=np.uint64)
    out = np.zeros(a.shape, dtype=np.int64)
    for k in range(64):
        out += ((a >> np.uint64(k)) & np.uint64(1)).astype(np.int64)
    return out


class Fock:
    """Spinful fermions on norb orbitals (modes = 2*norb)."""

    def __init__(self, norb):
        self.norb = norb
        self.nm = 2 * norb
        self.dim = 1 << self.nm
        self.states = np.arange(self.dim, dtype=np.int64)
        self._c = {}
        up_mask = (1 << norb) - 1
        self.nup = popcount(self.states & up_mask)
        self.ndn = popcount(self.states >> norb)
        self.sector = self.nup * (norb + 1) + self.ndn

    def mode(self, orb, spin):
        return orb + self.norb * spin

    def c(self, m):
        """annihilation operator of mode m (sparse, dim x dim)"""
        if m not in self._c:
            s = self.states
            has = ((s >> m) & 1) == 1
            src = s[has]
            dst = src ^ (1 << m)
            sign = 1 - 2 * (popcount(src & ((1 << m) - 1)) % 2)
            self._c[m] = sp.csr_matrix((sign.astype(float), (dst, src)), shape=(self.dim, self.dim))
        return self._c[m]

    def cd(self, m):
        return self.c(m).T.tocsr()

    def cdc(self, m1, m2):
        """c^dagger_{m1} c_{m2}"""
        return (self.cd(m1) @ self.c(m2)).tocsr()

    def hop(self, i, j):
        """sum_sigma c^dagger_{i,sigma} c_{j,sigma}  (orbital indices); NOT hermitian-completed"""
        return sum(self.cdc(self.mode(i, s), self.mode(j, s)) for s in (0, 1)).tocsr()

    def hopsym(self, i, j):
        """sum_sigma (c^dagger_i c_j + h.c.)"""
        h = self.hop(i, j)
        return (h + h.T).tocsr()

    def nspin(self, i, s):
        return self.cdc(self.mode(i, s), self.mode(i, s))

    def n(self, i):
        return (self.nspin(i, 0) + self.nspin(i, 1)).tocsr()

    def identity(self):
        return sp.identity(self.dim, format="csr")

    def sz2(self, i):
        """n_up - n_dn  (= 2 S^z)"""
        return (self.nspin(i, 0) - self.nspin(i, 1)).tocsr()


def thermal(H, ops, beta, sector, return_spectrum=False):
    """Z-normalised thermal expectation values.
    H: sparse hermitian (real), ops: dict name -> sparse matrix (must conserve `sector`), sector: int label per basis state.
    Returns dict name -> float (plus 'Z_rel' = Z exp(beta*Emin))."""
    H = sp.csr_matrix(H)
    labels = np.unique(sector)
    parts = []
    emin = np.inf
    for lab in labels:
        idx = np.nonzero(sector == lab)[0]
        Hs = H[idx][:, idx].toarray()
        # hermiticity check
        assert np.abs(Hs - Hs.T).max() < 1e-10 * max(1.0, np.abs(Hs).max()), "H not symmetric in sector"
        E, V = np.linalg.eigh(Hs)
        parts.append((idx, E, V))
        emin = min(emin, E[0])
    Z = 0.0
    acc = {k: 0.0 for k in ops}
    for idx, E, V in parts:
        w = np.exp(-beta * (E - emin))
        Z += w.sum()
        for k, O in ops.items():
            Os = sp.csr_matrix(O)[idx][:, idx]
            # sector conservation check (leakage out of the sector)
            full = sp.csr_matrix(O)[:, idx]
            leak = abs(full).sum() - abs(Os).sum()
            assert leak < 1e-9 * max(1.0, abs(full).sum()), f"operator {k} leaves the sector"
            OV = Os @ V
            acc[k] += float((w * np.einsum("ij,ij->j", V.conj(), OV).real).sum())
    out = {k: v / Z for k, v in acc.items()}
    out["_Z_rel"] = Z
    out["_Emin"] = emin
    return out


def lattice_sites(L1, L2):
    """cells (n1, n2) with n1 in 0..L1-1, n2 in 0..L2-1 (ALF: Latt%list); returns list and index map"""
    cells = [(n1, n2) for n2 in range(L2) for n1 in range(L1)]
    return cells, {c: i for i, c in enumerate(cells)}


def shift(cell, d, L1, L2):
    return ((cell[0] + d[0]) % L1, (cell[1] + d[1]) % L2)


def latt_displacements(L1, L2):
    return [(dx, dy) for dy in range(L2) for dx in range(L1)]


def eq_lattice_ops(F, L1, L2, norb_cell, orb_index, kinds=("Green", "SpinZ", "Den")):
    """Operators for the ALF equal-time lattice observables (Predefined_Obs_eq_*_measure, N_SUN=2, N_FL=1):
        Obs_Latt(R,a,b) = (1/Ns) sum_I  O(I,a ; J,b),   R = cell(I)-cell(J) (Lattice%imj(I,J) is I-J)
        Green : sum_sigma <c^dag_{I a s} c_{J b s}>
        SpinZ : <(n_up-n_dn)_{I a} (n_up-n_dn)_{J b}>      (= 4 <Sz Sz>)
        Den   : <n_{I a} n_{J b}>,    Latt0 = <n_{I a}>
    orb_index(cell, a) -> orbital number in F.  Returns dict name -> sparse operator and the key list."""
    cells, _ = lattice_sites(L1, L2)
    Ns = len(cells)
    ops = {}
    for (dx, dy) in latt_displacements(L1, L2):
        for a in range(norb_cell):
            for b in range(norb_cell):
                acc = {k: None for k in kinds}
                for I in cells:
                    J = ((I[0] - dx) % L1, (I[1] - dy) % L2)
                    i, j = orb_index(I, a), orb_index(J, b)
                    terms = {}
                    if "Green" in kinds:
                        terms["Green"] = F.hop(i, j)
                    if "SpinZ" in kinds:
                        terms["SpinZ"] = (F.sz2(i) @ F.sz2(j)).tocsr()
                    if "Den" in kinds:
                        terms["Den"] = (F.n(i) @ F.n(j)).tocsr()
                    for k, t in terms.items():
                        acc[k] = t if acc[k] is None else acc[k] + t
                for k in kinds:
                    ops[f"{k}_eq[{dx},{dy},{a},{b}]"] = (acc[k] / Ns).tocsr()
    if "Den" in kinds:
        for a in range(norb_cell):
            s = None
            for I in cells:
                t = F.n(orb_index(I, a))
                s = t if s is None else s + t
            ops[f"Den_eq0[{a}]"] = (s / Ns).tocsr()
    return ops


def alf_lattice_listr(L1, L2):
    """R list in the order of ALF's Lattice%list (the pyALF Lattice class reproduces it); needs py_alf."""
    from py_alf.lattice import Lattice
    lat = Lattice({"L1": [L1, 0.0], "L2": [0.0, L2], "a1": [1.0, 0.0], "a2": [0.0, 1.0]})
    return [tuple(int(x) for x in r) for r in lat.listr]


def to_alf_arrays(r, L1, L2, norb_cell, kinds=("Green", "SpinZ", "Den")):
    """r: result of thermal() with keys 'Kind_eq[dx,dy,a,b]' (+ 'Den_eq0[a]').
    Returns dict 'Kind_eqR' -> array [a][b][iR] in ALF listr order; Den is made connected like the pyALF analysis
    (N_Back = 1 / ana: subtract <n_a><n_b>), Green/SpinZ have no background."""
    import numpy as np
    listr = alf_lattice_listr(L1, L2)
    out = {}
    for k in kinds:
        arr = np.zeros((norb_cell, norb_cell, len(listr)))
        for iR, (i1, i2) in enumerate(listr):
            for a in range(norb_cell):
                for b in range(norb_cell):
                    v = r[f"{k}_eq[{i1 % L1},{i2 % L2},{a},{b}]"]
                    if k == "Den":
                        v = v - r[f"Den_eq0[{a}]"] * r[f"Den_eq0[{b}]"]
                    arr[a, b, iR] = v
        out[f"{k}_eqR"] = arr
    return out, listr
