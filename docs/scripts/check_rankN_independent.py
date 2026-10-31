#!/usr/bin/env python3
"""Independent numerical check of rank-r fast / delayed / submatrix updates in ALF conventions.

Conventions (per flavour, single flavour here):
  G = (1 + A)^{-1},  A = cyclic product of slice matrices cut at the current position.
  Wraps are similarity transforms G <- L G R (R = L^{-1}), L acting on rows P, R on columns P.
  In the eigenframe of vertex n (after type-1 wrap up / type-2 wrap down) a field change
  phi -> phi' multiplies A from the left by F = 1 + sum_{m<=r} Delta_m e_{P(m)} e_{P(m)}^T,
  Delta_m = exp(g (phi'-phi) E_m) - 1.

Schemes compared against the exact (1+A)^{-1}:
  std   : literal transcription of ALF Upgrade2 (sequential Sherman-Morrison x_v, y_v).
  dly   : delayed  G = G0 - X Y^T   (stored as X [N x k], Yt = Y^T [k x N]).
  sub   : submatrix  G = Gb - Gb[:,S] M (1-Gb)[S,:],  Gb = G0 - X Yt (delayed base),
          lazy reads of Gb, pending-site padding on partial overlap with non-diagonal wraps,
          conversion of (S,M) into delayed columns at hopping wraps (or when |S| > kmax).
Run with: python3 -I check_indep.py
"""
import numpy as np

rng = np.random.default_rng(20261031)
N = 10
dtau = 0.1


def herm(n, scale=1.0):
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return scale * (a + a.conj().T) / 2


def expmh(H, c):
    """exp(c*H) for Hermitian H, complex scalar c."""
    w, v = np.linalg.eigh(H)
    return (v * np.exp(c * w)) @ v.conj().T


# ----------------------------------------------------------------------------------------
# Operators (emulating Op_set): O Hermitian on P; non-diagonal -> eigen-decomposition with
# non-zero eigenvalues first, column 1 of U divided by det(U); diagonal -> U = 1, r = N_op.
# ----------------------------------------------------------------------------------------
class Op:
    def __init__(self, P, O, g):
        self.P = np.array(P)
        self.O = np.array(O, dtype=complex)
        self.g = g
        n = len(P)
        self.diag = bool(np.allclose(self.O, np.diag(np.diag(self.O)), atol=0.0))
        if self.diag:
            self.U = np.eye(n, dtype=complex)
            self.E = np.real(np.diag(self.O)).copy()
            self.r = n
        else:
            w, v = np.linalg.eigh(self.O)
            nz = [i for i in range(n) if abs(w[i]) > 1e-10]
            z = [i for i in range(n) if abs(w[i]) <= 1e-10]
            # ALF fills zero eigenvalues from the back; order among them irrelevant here
            order = nz + z[::-1]
            self.U = v[:, order].astype(complex)
            self.E = np.array([w[i] if abs(w[i]) > 1e-10 else 0.0 for i in order])
            self.r = len(nz)
            self.U[:, 0] /= np.linalg.det(self.U)
        assert np.allclose(self.U @ np.diag(self.E) @ self.U.conj().T, self.O)

    def V(self, phi):
        """embedded exp(g phi O) (N x N)."""
        M = np.eye(N, dtype=complex)
        M[np.ix_(self.P, self.P)] = self.U @ np.diag(np.exp(self.g * phi * self.E)) @ self.U.conj().T
        return M

    def wrap(self, kind, phi):
        """(L, R) blocks for the four ALF wraps, read off Op_Wrapup / Op_Wrapdo + ZSLGEMM."""
        U, Ud = self.U, self.U.conj().T
        e_p = np.diag(np.exp(self.g * phi * self.E))
        e_m = np.diag(np.exp(-self.g * phi * self.E))
        if kind == 'up1':   # Mat = e^{g phi E} U^dag Mat U e^{-g phi E}
            return e_p @ Ud, U @ e_m
        if kind == 'up2':   # Mat = U Mat U^dag
            return U, Ud
        if kind == 'do2':   # Mat = U^dag Mat U
            return Ud, U
        if kind == 'do1':   # Mat = U e^{-g phi E} Mat e^{g phi E} U^dag
            return U @ e_m, e_p @ Ud
        raise ValueError(kind)


def apply_LR(Mat, P, L, R):
    Mat[P, :] = L @ Mat[P, :]
    Mat[:, P] = Mat[:, P] @ R


# random rank-2 operator on 3 sites (zero eigenvalue -> r < N_op)
v1 = rng.normal(size=3) + 1j * rng.normal(size=3)
v2 = rng.normal(size=3) + 1j * rng.normal(size=3)
O_r2 = 0.6 * np.outer(v1, v1.conj()) - 0.4 * np.outer(v2, v2.conj())
t45 = 0.8 * np.exp(0.7j)
ops = [
    Op([3], [[1.0]], 0.30 + 0.10j),                                   # diag single site
    Op([3, 4], np.diag([0.7, -1.2]), 0.25 - 0.05j),                   # diag 2-site, overlaps site 3
    Op([4, 5], [[0, t45], [np.conj(t45), 0]], 0.35 + 0.08j),         # non-diag bond, r = 2
    Op([5, 6, 7], O_r2, 0.30 - 0.10j),                                # non-diag, r = 2 < 3
    Op([1, 2], [[1, 1], [1, 1]], 0.40 + 0.0j),                        # non-diag, eig 2,0 -> r = 1
    Op([4, 5], [[0.3, -0.5j], [0.5j, -0.2]], 0.30 + 0.12j),          # same bond again (repeated sites)
    Op([7, 2, 8, 3], herm(4, 0.5), 0.20 + 0.05j),                     # 4-site, unsorted P
    Op([3], [[-0.8]], 0.30 - 0.20j),                                  # site 3 again
    Op([0, 9], [[0.1, 0.6], [0.6, -0.3]], 0.30 + 0.0j),               # disjoint bond
    Op([6, 1], [[0.0, 0.9j], [-0.9j, 0.0]], 0.25 + 0.1j),             # bond touching padded sites
]
Nv = len(ops)
Ntau = 2
print("operators (P, diag, r, N_op):", [([int(x) for x in o.P], o.diag, o.r, len(o.P)) for o in ops])

T = herm(N, 1.0)
B_T = expmh(T, -dtau)
B_Ti = expmh(T, dtau)
Rest = expmh(herm(N, 1.0), -0.3) @ expmh(herm(N, 0.5), 0.2)

phi = rng.normal(size=(Ntau, Nv))           # continuous fields (type-3 like)
# fixed proposals and accept mask for each visit (up-sweep, down-sweep)
prop = {(d, t, n): rng.normal() for d in ('up', 'do') for t in range(Ntau) for n in range(Nv)}
mask = {(d, t, n): bool(rng.random() < 0.75) for d in ('up', 'do') for t in range(Ntau) for n in range(Nv)}


def slice_mat(t):
    """S_t = V_Nv ... V_1 B_T  (rightmost factor is wrapped first in the up-sweep)."""
    M = B_T.copy()
    for n in range(Nv):
        M = ops[n].V(phi[t, n]) @ M
    return M


def scratch_A(stage):
    """A from scratch at slice boundaries."""
    S = [slice_mat(t) for t in range(Ntau)]
    if stage == 'start':          # before up sweep of slice 0
        A = Rest
        for t in reversed(range(Ntau)):
            A = A @ S[t]
        return A
    if stage == 'after_up':       # S_{Ntau-1} ... S_0 Rest
        A = np.eye(N, dtype=complex)
        for t in reversed(range(Ntau)):
            A = A @ S[t]
        return A @ Rest
    raise ValueError


# ----------------------------------------------------------------------------------------
# Scheme std: literal ALF Upgrade2 (single flavour, ratio without alpha / gamma factors)
# ----------------------------------------------------------------------------------------
def upgrade2_ratio_mat(G, op, Delta):
    r = op.r
    P = op.P[:r]
    Mat = np.zeros((r, r), dtype=complex)
    for m in range(r):
        for n in range(r):
            Mat[n, m] = -Delta[m] * G[P[n], P[m]]
        Mat[m, m] = (1 + Delta[m]) + Mat[m, m]
    return Mat


def upgrade2_accept(G, op, Delta):
    """literal transcription of the accept branch of Upgrade2; returns new G."""
    r = op.r
    P = op.P[:r]
    u = np.zeros((N, r), dtype=complex)
    v = np.zeros((N, r), dtype=complex)
    for n in range(r):
        u[P[n], n] = Delta[n]
        v[:, n] = -G[P[n], :]
        v[P[n], n] = 1 - G[P[n], P[n]]
    x_v = np.zeros((N, r), dtype=complex)
    y_v = np.zeros((N, r), dtype=complex)
    i = P[0]
    x_v[i, 0] = u[i, 0] / (1 + v[i, 0] * u[i, 0])
    y_v[:, 0] = v[:, 0]
    for n in range(1, r):
        x_v[:, n] = u[:, n]
        y_v[:, n] = v[:, n]
        Z = 1 + u[P[n], n] * v[P[n], n]
        syu = -(y_v[:, :n].T @ u[:, n])       # zgemv 'T' : plain transpose
        sxv = -(x_v[:, :n].T @ v[:, n])
        x_v[:, n] += x_v[:, :n] @ syu
        y_v[:, n] += y_v[:, :n] @ sxv
        Z = Z - np.sum(syu * sxv)
        x_v[:, n] *= 1 / Z
    # identity check:  x_v y_v^T == u Mat^{-1} v^T  (Mat as in the ratio)
    Mat = upgrade2_ratio_mat(G, op, Delta)
    assert np.allclose(x_v @ y_v.T, u @ np.linalg.solve(Mat, v.T), atol=1e-11)
    xp = G[:, P] @ x_v[P, :]
    return G - xp @ y_v.T


# ----------------------------------------------------------------------------------------
# Scheme dly: delayed  G = G0 - X Yt
# ----------------------------------------------------------------------------------------
class Delayed:
    def __init__(self, G, nd):
        self.G0 = G.copy()
        self.X = np.zeros((N, 0), dtype=complex)
        self.Yt = np.zeros((0, N), dtype=complex)
        self.nd = nd
        self.flushes = 0

    def G(self):
        return self.G0 - self.X @ self.Yt

    def read(self, rows, cols):
        return self.G0[np.ix_(rows, cols)] - self.X[rows, :] @ self.Yt[:, cols]

    def wrap(self, P, L, R, diag=False):
        apply_LR(self.G0, P, L, R)            # what ALF does to GR
        self.X[P, :] = L @ self.X[P, :]       # X  <- L X
        self.Yt[:, P] = self.Yt[:, P] @ R     # Y^T <- Y^T R   (Y <- R^T Y)

    def hop(self, Lf, Rf):                    # G <- Lf G Rf, dense
        self.G0 = Lf @ self.G0 @ Rf
        self.X = Lf @ self.X
        self.Yt = self.Yt @ Rf

    def ratio_mat(self, Q, Delta):
        GQQ = self.read(Q, Q)
        return np.eye(len(Q)) + (np.eye(len(Q)) - GQQ) @ np.diag(Delta)

    def accept(self, Q, Delta, Mat):
        allc = np.arange(N)
        Gcol = self.read(allc, Q)                         # G[:,Q]
        Grow = self.read(Q, allc)                         # G[Q,:]
        EQt = np.zeros((len(Q), N), dtype=complex)
        EQt[np.arange(len(Q)), Q] = 1
        xnew = Gcol @ np.diag(Delta) @ np.linalg.inv(Mat)  # N x r
        ytnew = EQt - Grow                                 # (1-G)[Q,:]
        self.X = np.hstack([self.X, xnew])
        self.Yt = np.vstack([self.Yt, ytnew])
        if self.X.shape[1] >= self.nd:
            self.flush()

    def append_cols(self, Xn, Ytn):
        self.X = np.hstack([self.X, Xn])
        self.Yt = np.vstack([self.Yt, Ytn])
        if self.X.shape[1] >= self.nd:
            self.flush()

    def flush(self):
        self.G0 = self.G()
        self.X = np.zeros((N, 0), dtype=complex)
        self.Yt = np.zeros((0, N), dtype=complex)
        self.flushes += 1


# ----------------------------------------------------------------------------------------
# Scheme sub: submatrix on top of a delayed base
# ----------------------------------------------------------------------------------------
class Submatrix:
    def __init__(self, G, kmax, nd, flush_at_hop=False):
        self.base = Delayed(G, nd)
        self.S = []                                      # pending sites (ordered list)
        self.M = np.zeros((0, 0), dtype=complex)
        self.kmax = kmax
        self.flush_at_hop = flush_at_hop
        self.pads = 0
        self.conversions = 0

    def G(self):
        Gb = self.base.G()
        if not self.S:
            return Gb
        S = np.array(self.S)
        k = len(S)
        ESt = np.zeros((k, N), dtype=complex)
        ESt[np.arange(k), S] = 1
        return Gb - Gb[:, S] @ self.M @ (ESt - Gb[S, :])

    def pad(self, sites):
        new = [s for s in sites if s not in self.S]
        if new:
            k = len(self.S)
            Mn = np.zeros((k + len(new), k + len(new)), dtype=complex)
            Mn[:k, :k] = self.M
            self.M = Mn
            self.S = self.S + new
        return [self.S.index(s) for s in sites]

    def wrap(self, P, L, R, diag):
        self.base.wrap(P, L, R)
        if not self.S:
            return
        Pl = list(P)
        inter = [p for p in Pl if p in self.S]
        if not inter:
            return                                           # disjoint: M unchanged
        if diag:                                             # L, R diagonal: rescale only
            for j, p in enumerate(Pl):
                if p in self.S:
                    i = self.S.index(p)
                    self.M[i, :] *= L[j, j]
                    self.M[:, i] *= R[j, j]
            return
        if len(inter) < len(Pl):                             # partial overlap -> pad
            self.pads += 1
        pos = self.pad(Pl)
        self.M[pos, :] = L @ self.M[pos, :]                  # M <- L_S M R_S
        self.M[:, pos] = self.M[:, pos] @ R

    def ratio_and_cache(self, Q, Delta):
        Q = list(Q)
        r = len(Q)
        S = self.S
        if S:
            GbQS = self.base.read(Q, S)
            GbSQ = self.base.read(S, Q)
            ISQ = np.array([[1.0 if s == q else 0.0 for q in Q] for s in S])
            t = self.M @ (ISQ - GbSQ)                        # k x r
            GQQ = self.base.read(Q, Q) - GbQS @ t
        else:
            GQQ = self.base.read(Q, Q)
        Mat = np.eye(r) + (np.eye(r) - GQQ) @ np.diag(Delta)
        return Mat

    def accept(self, Q, Delta, Mat):
        Q = list(Q)
        r = len(Q)
        posQ = self.pad(Q)                                   # S' = S u Q, M padded with zeros
        S = self.S
        k = len(S)
        JQ = np.zeros((k, r))
        JQ[posQ, np.arange(r)] = 1
        GbSQ = self.base.read(S, Q)
        GbQS = self.base.read(Q, S)
        a = JQ - self.M @ (JQ - GbSQ)                        # G[:,Q]  = Gb[:,S'] a
        bt = JQ.T + GbQS @ self.M                            # (1-G)[Q,:] = bt (1-Gb)[S',:]
        self.M = self.M + a @ np.diag(Delta) @ np.linalg.solve(Mat, bt)
        if len(self.S) > self.kmax:
            self.convert()

    def convert(self):
        """(S, M) -> delayed columns of the base:  X += Gb[:,S] M,  Y^T += (1-Gb)[S,:]."""
        if not self.S:
            return
        S = np.array(self.S)
        k = len(S)
        allc = np.arange(N)
        ESt = np.zeros((k, N), dtype=complex)
        ESt[np.arange(k), S] = 1
        Xn = self.base.read(allc, list(S)) @ self.M
        Ytn = ESt - self.base.read(list(S), allc)
        self.S, self.M = [], np.zeros((0, 0), dtype=complex)
        self.conversions += 1
        self.base.append_cols(Xn, Ytn)

    def hop(self, Lf, Rf):
        self.convert()
        if self.flush_at_hop:
            self.base.flush()
        self.base.hop(Lf, Rf)


# ----------------------------------------------------------------------------------------
# Driver: up sweep over Ntau slices then down sweep, following WRAPGRUP / WRAPGRDO
# ----------------------------------------------------------------------------------------
def run(quiet=False):
    phi_save = phi.copy()
    A = scratch_A('start')
    Gex = np.linalg.inv(np.eye(N) + A)
    Gstd = Gex.copy()
    schemes = {
        'dly(nd=5)': Delayed(Gex, 5),
        'dly(nd=inf)': Delayed(Gex, 10 ** 9),
        'sub(kmax=6,nd=8)': Submatrix(Gex, 6, 8),
        'sub(kmax=inf,nd=inf)': Submatrix(Gex, 10 ** 9, 10 ** 9),
        'sub(flush@hop)': Submatrix(Gex, 10 ** 9, 10 ** 9, flush_at_hop=True),
    }
    err = {k: 0.0 for k in list(schemes) + ['std']}
    rerr = {k: 0.0 for k in list(schemes) + ['std']}
    nacc = 0

    def all_wrap(op, kind, ph):
        nonlocal A, Gstd
        L, R = op.wrap(kind, ph)
        Lf = np.eye(N, dtype=complex); Lf[np.ix_(op.P, op.P)] = L
        Rf = np.eye(N, dtype=complex); Rf[np.ix_(op.P, op.P)] = R
        assert np.allclose(Lf @ Rf, np.eye(N))
        A = Lf @ A @ Rf
        apply_LR(Gstd, op.P, L, R)
        for s in schemes.values():
            s.wrap(op.P, L, R, op.diag) if isinstance(s, Submatrix) else s.wrap(op.P, L, R)

    def all_hop(Lf, Rf):
        nonlocal A, Gstd
        A = Lf @ A @ Rf
        Gstd = Lf @ Gstd @ Rf
        for s in schemes.values():
            s.hop(Lf, Rf)

    def vertex(t, n, d):
        nonlocal A, Gstd, nacc
        op = ops[n]
        r = op.r
        Q = op.P[:r]
        phn = prop[(d, t, n)]
        Delta = np.exp(op.g * (phn - phi[t, n]) * op.E[:r]) - 1
        F = np.eye(N, dtype=complex)
        F[Q, Q] += Delta
        ratio_ex = np.linalg.det(np.eye(N) + F @ A) / np.linalg.det(np.eye(N) + A)
        Mat_std = upgrade2_ratio_mat(Gstd, op, Delta)
        rat = {'std': np.linalg.det(Mat_std)}
        mats = {}
        for k, s in schemes.items():
            Mk = s.ratio_and_cache(Q, Delta) if isinstance(s, Submatrix) else s.ratio_mat(Q, Delta)
            mats[k] = Mk
            rat[k] = np.linalg.det(Mk)
        for k in rat:
            rerr[k] = max(rerr[k], abs(rat[k] / ratio_ex - 1))
        if mask[(d, t, n)]:
            nacc += 1
            A = F @ A
            phi[t, n] = phn
            Gstd = upgrade2_accept(Gstd, op, Delta)
            for k, s in schemes.items():
                s.accept(Q, Delta, mats[k])
        Gref = np.linalg.inv(np.eye(N) + A)
        err['std'] = max(err['std'], np.abs(Gstd - Gref).max())
        for k, s in schemes.items():
            err[k] = max(err[k], np.abs(s.G() - Gref).max())

    # ---- up sweep (WRAPGRUP per slice) ----
    for t in range(Ntau):
        all_hop(B_T, B_Ti)                                  # HOP_MOD_mmthr / mmthl_m1
        for n in range(Nv):
            all_wrap(ops[n], 'up1', phi[t, n])              # type 1, old field
            vertex(t, n, 'up')
            all_wrap(ops[n], 'up2', phi[t, n])              # type 2
    Ascr = scratch_A('after_up')
    frame_err_up = np.abs(A - Ascr).max() / np.abs(Ascr).max()
    # ---- down sweep (WRAPGRDO per slice) ----
    for t in reversed(range(Ntau)):
        for n in reversed(range(Nv)):
            all_wrap(ops[n], 'do2', phi[t, n])              # type 2
            vertex(t, n, 'do')
            all_wrap(ops[n], 'do1', phi[t, n])              # type 1, new field
        all_hop(B_Ti, B_T)                                  # Hop_mod_mmthl / mmthr_m1
    Ascr = scratch_A('start')
    frame_err_do = np.abs(A - Ascr).max() / np.abs(Ascr).max()
    Gfinal = np.linalg.inv(np.eye(N) + Ascr)

    if quiet:
        phi[:] = phi_save
        return None, err
    print(f"vertex visits: {2 * Ntau * Nv}, accepted: {nacc}")
    print(f"frame check (tracked A vs product from scratch): after up {frame_err_up:.2e}, after down {frame_err_do:.2e}")
    print(f"{'scheme':24s} {'max|G-G_exact|':>16s} {'max|ratio/exact-1|':>20s} {'final |G-G_scratch|':>20s}")
    print(f"{'std (Upgrade2)':24s} {err['std']:16.2e} {rerr['std']:20.2e} {np.abs(Gstd - Gfinal).max():20.2e}")
    for k, s in schemes.items():
        extra = ''
        if isinstance(s, Submatrix):
            extra = f"  [pads={s.pads}, conversions={s.conversions}, base flushes={s.base.flushes}]"
        else:
            extra = f"  [flushes={s.flushes}]"
        print(f"{k:24s} {err[k]:16.2e} {rerr[k]:20.2e} {np.abs(s.G() - Gfinal).max():20.2e}{extra}")
    ok = max(list(err.values()) + list(rerr.values())) < 1e-9 and max(frame_err_up, frame_err_do) < 1e-12
    if not quiet:
        print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    phi[:] = phi_save
    return ok, err


def negative_controls():
    """Deliberately wrong variants must fail (shows the check is sensitive)."""
    print("negative controls (max|G-G_exact| of the broken scheme; must be O(1e-3) or larger):")
    # 1) Y transformed with R^dagger instead of R^T  (conjugate-transpose confusion)
    orig = Delayed.wrap
    def bad_wrap(self, P, L, R, diag=False):
        apply_LR(self.G0, P, L, R)
        self.X[P, :] = L @ self.X[P, :]
        self.Yt[:, P] = (R.conj().T @ self.Yt[:, P].T).T
    Delayed.wrap = bad_wrap
    _, e = run(quiet=True)
    Delayed.wrap = orig
    print(f"  Y <- R^dag Y instead of R^T Y      : dly(nd=inf) {e['dly(nd=inf)']:.2e}")
    # 2) partial overlap handled without padding (only the intersection block transformed)
    origw = Submatrix.wrap
    def nopad_wrap(self, P, L, R, diag):
        self.base.wrap(P, L, R)
        Pl = list(P)
        idx = [j for j, p in enumerate(Pl) if p in self.S]
        if not idx:
            return
        pos = [self.S.index(Pl[j]) for j in idx]
        self.M[pos, :] = L[np.ix_(idx, idx)] @ self.M[pos, :]
        self.M[:, pos] = self.M[:, pos] @ R[np.ix_(idx, idx)]
    Submatrix.wrap = nopad_wrap
    _, e = run(quiet=True)
    Submatrix.wrap = origw
    print(f"  no padding on partial overlap      : sub(kmax=inf,nd=inf) {e['sub(kmax=inf,nd=inf)']:.2e}")
    # 3) hopping wrap: X <- B X but Y left untouched
    origh = Delayed.hop
    def bad_hop(self, Lf, Rf):
        self.G0 = Lf @ self.G0 @ Rf
        self.X = Lf @ self.X
    Delayed.hop = bad_hop
    _, e = run(quiet=True)
    Delayed.hop = origh
    print(f"  hop: Y^T not multiplied by B^-1    : dly(nd=inf) {e['dly(nd=inf)']:.2e}")


if __name__ == '__main__':
    run()
    negative_controls()
