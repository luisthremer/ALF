#!/usr/bin/env python3
"""Numerical check of the rank-N delayed and generalised submatrix updates (docs/algorithm-notes.md, sections 9-12).

One or two time slices of vertices are processed exactly as ALF's WRAPGRUP / WRAPGRDO do it, vertex by vertex:
  up  : type-1 wrap (old field)  -> propose/accept in the eigenframe -> type-2 wrap
  down: type-2 wrap              -> propose/accept in the eigenframe -> type-1 wrap (new field)
with diagonal and non-diagonal multi-site vertices (overlapping supports, repeated sites, zero eigenvalues,
complex couplings). Four G's are carried:
  * exact    : G = (1+A)^{-1}, A transformed as A <- L A R (wraps), A <- F A (accepts)
  * standard : Upgrade2 (Woodbury on the full G)
  * delayed  : G = G0 - X Y^T, hooks X[P] <- L X[P], Y[P] <- R^T Y[P]
  * submatrix: G = G0 - G0[:,S] Ginv Q (E_S^T - G0[S,:]) with dense Q, padding, paths B and W
and optionally the delayed factors are carried across a hopping wrap between two slices (cross-slice).
The accept/reject mask is fixed; every scheme must give the same ratios and the same G.
Submatrix: the adopted M-form (class SubmatrixM, state S and M = Gamma~^{-1} Q) runs by default; the rejected
Gamma~^{-1} + Q form with paths B/W (class Submatrix) is compared with --both-forms (it loses accuracy, up to
~1e-7, in the "mixed" model with many path-W updates; notes Sec. 12).
Usage: python3 -I check_rankN.py [--both-forms]   (prints max relative errors, exit 1 on failure)
"""
import sys
import numpy as np

TOL = 1e-10
rng = np.random.default_rng(20261031)


def herm(n, diag):
    if diag:
        return np.diag(rng.choice([-1.0, 1.0, 0.5, 2.0], size=n)).astype(complex)
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return (a + a.conj().T) / 2


def make_vertex(Ndim, sites, diag, zero_eig=False):
    n = len(sites)
    O = herm(n, diag)
    if zero_eig:  # rank-deficient operator: one zero eigenvalue
        w, V = np.linalg.eigh(O)
        w[0] = 0.0
        O = V @ np.diag(w) @ V.conj().T
        if diag:
            O = np.diag(np.diag(O))
    if diag:
        E = np.real(np.diag(O)).copy()
        U = np.eye(n, dtype=complex)
    else:
        E, U = np.linalg.eigh(O)
    order = np.argsort(np.abs(E) <= 1e-12, kind="stable")  # non-zero eigenvalues first (Op_set)
    E, U = E[order], U[:, order]
    r = int(np.sum(np.abs(E) > 1e-12))
    g = 0.3 * np.exp(1j * rng.uniform(0, 0.6))
    return dict(P=np.array(sites), E=E, U=U, r=r, g=g, diag=diag)


def wrap_mats(v, phi, kind, direction):
    """(L, R) on P for ALF's Op_Wrapup/Op_Wrapdo; G[P-block] <- L G R (L on rows P, R on cols P)."""
    e = np.exp(v["g"] * phi * v["E"])
    U = v["U"]
    if direction == "up":
        if kind == 1:
            return np.diag(e) @ U.conj().T, U @ np.diag(1 / e)
        return U, U.conj().T
    if kind == 2:
        return U.conj().T, U
    return U @ np.diag(1 / e), np.diag(e) @ U.conj().T


def embed(Ndim, P, Ml, Mr):
    L = np.eye(Ndim, dtype=complex)
    R = np.eye(Ndim, dtype=complex)
    L[np.ix_(P, P)] = Ml
    R[np.ix_(P, P)] = Mr
    return L, R


class Delayed:
    def __init__(self, G0):
        self.G0 = G0.copy()
        self.X = np.zeros((G0.shape[0], 0), complex)
        self.Y = np.zeros((G0.shape[0], 0), complex)

    def G(self):
        return self.G0 - self.X @ self.Y.T

    def wrap(self, P, Ml, Mr):  # ALF wraps G0; the hook transforms the pending factors
        L, R = embed(self.G0.shape[0], P, Ml, Mr)
        self.G0 = L @ self.G0 @ R
        self.X[P, :] = Ml @ self.X[P, :]
        self.Y[P, :] = Mr.T @ self.Y[P, :]

    def gpp(self, Pa):
        return self.G0[np.ix_(Pa, Pa)] - self.X[Pa, :] @ self.Y[Pa, :].T

    def accept(self, Pa, Delta, Mat):
        Gc = self.G0[:, Pa] - self.X @ self.Y[Pa, :].T
        Gr = self.G0[Pa, :] - self.X[Pa, :] @ self.Y.T
        Xn = Gc @ np.diag(Delta) @ np.linalg.inv(Mat)
        Yt = -Gr
        Yt[np.arange(len(Pa)), Pa] += 1.0
        self.X = np.hstack([self.X, Xn])
        self.Y = np.hstack([self.Y, Yt.T])

    def flush(self):
        self.G0 = self.G()
        self.X = self.X[:, :0]
        self.Y = self.Y[:, :0]

    def hop(self, B, Binv):  # cross-slice: G <- B G B^{-1}
        self.G0 = B @ self.G0 @ Binv
        self.X = B @ self.X
        self.Y = Binv.T @ self.Y


class Submatrix:
    """G = G0 - G0[:,S] Ginv Q (E_S^T - G0[S,:]),  Ginv = (I + Q (I-G0)_SS)^{-1}."""

    def __init__(self, G0):
        self.G0 = G0.copy()
        self.S = []
        self.Q = np.zeros((0, 0), complex)
        self.Ginv = np.zeros((0, 0), complex)
        self.npad = 0
        self.npathB = 0
        self.npathW = 0

    def ES(self):
        E = np.zeros((self.G0.shape[0], len(self.S)))
        E[self.S, np.arange(len(self.S))] = 1.0
        return E

    def G(self):
        if not self.S:
            return self.G0.copy()
        S = self.S
        return self.G0 - self.G0[:, S] @ self.Ginv @ self.Q @ (self.ES().T - self.G0[S, :])

    def pad(self, sites):
        new = [s for s in sites if s not in self.S]
        if not new or not self.S:
            if new and not self.S:
                pass
            if not new:
                return
        k, m = len(self.S), len(new)
        TR = self.Ginv @ self.Q @ self.G0[np.ix_(self.S, new)] if k else np.zeros((0, m), complex)
        Gi = np.zeros((k + m, k + m), complex)
        Gi[:k, :k] = self.Ginv
        Gi[:k, k:] = TR
        Gi[k:, k:] = np.eye(m)
        Q = np.zeros((k + m, k + m), complex)
        Q[:k, :k] = self.Q
        self.Ginv, self.Q, self.S = Gi, Q, self.S + new
        self.npad += m

    def prewrap(self, P, diag):
        inter = [p for p in P if p in self.S]
        if not diag and inter and len(inter) < len(P):
            self.pad(list(P))

    def wrap(self, P, Ml, Mr):
        L, R = embed(self.G0.shape[0], P, Ml, Mr)
        self.G0 = L @ self.G0 @ R
        idx = [i for i, p in enumerate(P) if p in self.S]
        if not idx:
            return
        sl = [self.S.index(P[i]) for i in idx]
        k = len(self.S)
        LS = np.eye(k, dtype=complex)
        RS = np.eye(k, dtype=complex)
        LS[np.ix_(sl, sl)] = Ml[np.ix_(idx, idx)]
        RS[np.ix_(sl, sl)] = Mr[np.ix_(idx, idx)]
        if len(idx) < len(P):  # partial overlap: allowed only for diagonal L, R
            assert np.allclose(Ml, np.diag(np.diag(Ml))) and np.allclose(Mr, np.diag(np.diag(Mr)))
        self.Ginv = LS @ self.Ginv @ RS
        self.Q = LS @ self.Q @ RS

    def gpp(self, Pa):
        if not self.S:
            return self.G0[np.ix_(Pa, Pa)].copy()
        S = self.S
        ESP = self.ES()[Pa, :].T  # k x r
        return self.G0[np.ix_(Pa, Pa)] - self.G0[np.ix_(Pa, S)] @ self.Ginv @ self.Q @ (ESP - self.G0[np.ix_(S, Pa)])

    def accept(self, Pa, Delta, Mat):
        r = len(Pa)
        D = np.diag(Delta)
        if not any(p in self.S for p in Pa):  # path B: border
            self.npathB += 1
            k = len(self.S)
            if k:
                S = self.S
                a = self.Ginv @ self.Q @ self.G0[np.ix_(S, Pa)]
                bT = self.G0[np.ix_(Pa, S)] @ self.Ginv
                Gpp = self.G0[np.ix_(Pa, Pa)] + self.G0[np.ix_(Pa, S)] @ a
            else:
                a = np.zeros((0, r), complex)
                bT = np.zeros((r, 0), complex)
                Gpp = self.G0[np.ix_(Pa, Pa)]
            Sig = np.eye(r) + D @ (np.eye(r) - Gpp)
            Si = np.linalg.inv(Sig)
            Gi = np.zeros((k + r, k + r), complex)
            Gi[:k, :k] = self.Ginv + a @ Si @ D @ bT
            Gi[:k, k:] = a @ Si
            Gi[k:, :k] = Si @ D @ bT
            Gi[k:, k:] = Si
            Q = np.zeros((k + r, k + r), complex)
            Q[:k, :k] = self.Q
            Q[k:, k:] = D
            self.Ginv, self.Q, self.S = Gi, Q, self.S + list(Pa)
            return
        # path W: pad the missing active sites, then Woodbury on Ginv
        self.npathW += 1
        self.pad(list(Pa))
        S = self.S
        k = len(S)
        sl = [S.index(p) for p in Pa]
        ES = self.ES()
        GkPS = self.G0[np.ix_(Pa, S)] - (self.G0[np.ix_(Pa, S)] @ self.Ginv @ self.Q) @ (np.eye(k) - self.G0[np.ix_(S, S)])
        Wr = ES[Pa, :] - GkPS  # (I - G_k)(P,S)
        MatW = np.eye(r) + Wr[:, sl] @ D
        assert np.allclose(MatW, Mat, atol=1e-11), "path W: Mat from (I-G_k)(P,S) differs from the proposal"
        self.Ginv = self.Ginv - (self.Ginv[:, sl] @ D) @ np.linalg.inv(Mat) @ Wr
        Qs = self.Q[sl, :].copy()
        Qs[np.arange(r), sl] += 1.0
        self.Q[sl, :] += D @ Qs

    def flush(self):
        self.G0 = self.G()
        self.S, self.Q, self.Ginv = [], np.zeros((0, 0), complex), np.zeros((0, 0), complex)

    def to_delayed(self):
        d = Delayed(self.G0)
        if self.S:
            d.X = self.G0[:, self.S] @ self.Ginv @ self.Q
            d.Y = (self.ES().T - self.G0[self.S, :]).T
        return d


class SubmatrixM:
    """M-form (adopted, D3): G = G0 - G0[:,S] M (E_S^T - G0[S,:]); one accept formula for any overlap."""

    def __init__(self, G0):
        self.G0 = G0.copy()
        self.S = []
        self.M = np.zeros((0, 0), complex)
        self.npad = 0
        self.npathB = 0
        self.npathW = 0

    def J(self, Pa):
        Jm = np.zeros((len(self.S), len(Pa)))
        for m, p in enumerate(Pa):
            Jm[self.S.index(p), m] = 1.0
        return Jm

    def G(self):
        if not self.S:
            return self.G0.copy()
        E = np.zeros((self.G0.shape[0], len(self.S)))
        E[self.S, np.arange(len(self.S))] = 1.0
        return self.G0 - self.G0[:, self.S] @ self.M @ (E.T - self.G0[self.S, :])

    def pad(self, sites):
        new = [s for s in sites if s not in self.S]
        if not new:
            return
        k, m = len(self.S), len(new)
        Mn = np.zeros((k + m, k + m), complex)
        Mn[:k, :k] = self.M
        self.M, self.S = Mn, self.S + new
        self.npad += m

    def prewrap(self, P, diag):
        inter = [p for p in P if p in self.S]
        if not diag and inter and len(inter) < len(P):
            self.pad(list(P))

    def wrap(self, P, Ml, Mr):
        L, R = embed(self.G0.shape[0], P, Ml, Mr)
        self.G0 = L @ self.G0 @ R
        idx = [i for i, p in enumerate(P) if p in self.S]
        if not idx:
            return
        sl = [self.S.index(P[i]) for i in idx]
        k = len(self.S)
        LS = np.eye(k, dtype=complex)
        RS = np.eye(k, dtype=complex)
        LS[np.ix_(sl, sl)] = Ml[np.ix_(idx, idx)]
        RS[np.ix_(sl, sl)] = Mr[np.ix_(idx, idx)]
        if len(idx) < len(P):
            assert np.allclose(Ml, np.diag(np.diag(Ml))) and np.allclose(Mr, np.diag(np.diag(Mr)))
        self.M = LS @ self.M @ RS

    def gpp(self, Pa):
        if not self.S:
            return self.G0[np.ix_(Pa, Pa)].copy()
        Jm = np.zeros((len(self.S), len(Pa)))
        for m, p in enumerate(Pa):
            if p in self.S:
                Jm[self.S.index(p), m] = 1.0
        return self.G0[np.ix_(Pa, Pa)] - self.G0[np.ix_(Pa, self.S)] @ self.M @ (Jm - self.G0[np.ix_(self.S, Pa)])

    def accept(self, Pa, Delta, Mat):
        if any(p in self.S for p in Pa):
            self.npathW += 1
        else:
            self.npathB += 1
        self.pad(list(Pa))
        S = self.S
        Jm = self.J(Pa)
        a = Jm - self.M @ (Jm - self.G0[np.ix_(S, Pa)])
        bT = Jm.T + self.G0[np.ix_(Pa, S)] @ self.M
        self.M = self.M + a @ np.diag(Delta) @ np.linalg.inv(Mat) @ bT

    def flush(self):
        self.G0 = self.G()
        self.S, self.M = [], np.zeros((0, 0), complex)

    def to_delayed(self):
        d = Delayed(self.G0)
        if self.S:
            E = np.zeros((self.G0.shape[0], len(self.S)))
            E[self.S, np.arange(len(self.S))] = 1.0
            d.X = self.G0[:, self.S] @ self.M
            d.Y = (E.T - self.G0[self.S, :]).T
        return d


def run(Ndim, vertices, direction, accept_mask, nd_cap=10**9, cross_slice=False, seed_fields=None, sub_cls=None):
    A = rng.normal(size=(Ndim, Ndim)) + 1j * rng.normal(size=(Ndim, Ndim))
    A = 0.6 * A / np.sqrt(Ndim) + np.eye(Ndim) * 0.4
    G = np.linalg.inv(np.eye(Ndim) + A)
    std = G.copy()
    dl = Delayed(G)
    sm = (sub_cls or Submatrix)(G)
    fields = [rng.choice([-1.0, 1.0, -2.0, 2.0]) for _ in vertices]
    nslices = 2 if cross_slice else 1
    T = rng.normal(size=(Ndim, Ndim)) * 0.1
    T = (T + T.T) / 2
    w, V = np.linalg.eigh(T)
    B = V @ np.diag(np.exp(-w)) @ V.T
    Binv = V @ np.diag(np.exp(w)) @ V.T
    errs = dict(ratio_dl=0.0, ratio_sm=0.0, G_dl=0.0, G_sm=0.0, G_std=0.0, cross=0.0)
    iacc = 0
    for sl in range(nslices):
        if sl > 0:  # hopping wrap between slices
            A = B @ A @ Binv
            std = B @ std @ Binv
            dl.hop(B, Binv)
            sm = sm.to_delayed()  # submatrix state converted to delayed form, then continued as delayed
            sm.hop(B, Binv)
        order = range(len(vertices)) if direction == "up" else range(len(vertices) - 1, -1, -1)
        for n in order:
            v = vertices[n]
            P = v["P"]
            first = 1 if direction == "up" else 2
            second = 2 if direction == "up" else 1
            phi_old = fields[n]
            if hasattr(sm, 'prewrap'):
                sm.prewrap(list(P), v["diag"])
            Ml, Mr = wrap_mats(v, phi_old, first, direction)
            L, R = embed(Ndim, P, Ml, Mr)
            A = L @ A @ R
            std = L @ std @ R
            dl.wrap(P, Ml, Mr)
            sm.wrap(P, Ml, Mr)
            # proposal in the eigenframe of vertex n
            phi_new = -phi_old if abs(phi_old) == 1 else phi_old / 2.0
            r = v["r"]
            Pa = P[:r]
            Delta = np.exp(v["g"] * (phi_new - phi_old) * v["E"][:r]) - 1.0
            D = np.diag(Delta)
            if r > 0:
                Gpp_std = std[np.ix_(Pa, Pa)]
                Mat_std = np.eye(r) + (np.eye(r) - Gpp_std) @ D
                R_std = np.linalg.det(Mat_std)
                Mat_dl = np.eye(r) + (np.eye(r) - dl.gpp(Pa)) @ D
                Mat_sm = np.eye(r) + (np.eye(r) - sm.gpp(Pa)) @ D
                errs["ratio_dl"] = max(errs["ratio_dl"], abs(np.linalg.det(Mat_dl) / R_std - 1))
                errs["ratio_sm"] = max(errs["ratio_sm"], abs(np.linalg.det(Mat_sm) / R_std - 1))
                acc = accept_mask[iacc % len(accept_mask)]
                iacc += 1
                if acc:
                    F = np.eye(Ndim, dtype=complex)
                    F[Pa, Pa] += Delta
                    A = F @ A
                    std = std - std[:, Pa] @ D @ np.linalg.inv(Mat_std) @ (np.eye(Ndim)[Pa, :] - std[Pa, :])
                    dl.accept(Pa, Delta, Mat_dl)
                    sm.accept(list(Pa), Delta, Mat_sm)
                    fields[n] = phi_new
                    if dl.X.shape[1] >= nd_cap:
                        dl.flush()
                    if hasattr(sm, 'prewrap') and len(sm.S) >= nd_cap:
                        sm.flush()
            if hasattr(sm, 'prewrap'):
                sm.prewrap(list(P), v["diag"])
            Ml, Mr = wrap_mats(v, fields[n], second, direction)
            L, R = embed(Ndim, P, Ml, Mr)
            A = L @ A @ R
            std = L @ std @ R
            dl.wrap(P, Ml, Mr)
            sm.wrap(P, Ml, Mr)
    Gex = np.linalg.inv(np.eye(Ndim) + A)
    scale = np.abs(Gex).max()
    errs["G_std"] = np.abs(std - Gex).max() / scale
    errs["G_dl"] = np.abs(dl.G() - std).max() / scale
    errs["G_sm"] = np.abs(sm.G() - std).max() / scale
    info = dict(npad=getattr(sm, "npad", -1), pathB=getattr(sm, "npathB", -1), pathW=getattr(sm, "npathW", -1))
    return errs, info


def model(Ndim, kind):
    V = []
    if kind == "hubbard":  # rank-1 diagonal, each site once
        for i in rng.permutation(Ndim):
            V.append(make_vertex(Ndim, [i], True))
    elif kind == "jz":  # diagonal rank-2 bonds (overlapping supports)
        for i in range(Ndim):
            V.append(make_vertex(Ndim, [i, (i + 1) % Ndim], True))
    elif kind == "sasb_zz":  # diagonal rank-4, overlapping
        for i in range(0, Ndim, 2):
            V.append(make_vertex(Ndim, [i, (i + 1) % Ndim, (i + 2) % Ndim, (i + 3) % Ndim], True))
    elif kind == "tv":  # non-diagonal bonds, overlapping
        for i in range(Ndim):
            V.append(make_vertex(Ndim, [i, (i + 1) % Ndim], False))
    elif kind == "tv_families":  # disjoint bond families (even bonds, then odd bonds)
        for st in (0, 1):
            for i in range(st, Ndim - 1, 2):
                V.append(make_vertex(Ndim, [i, i + 1], False))
    elif kind == "kondo":  # on-site rank-1 + non-diagonal bonds with a zero eigenvalue
        for i in range(Ndim):
            V.append(make_vertex(Ndim, [i], True))
        for i in range(0, Ndim - 2, 2):
            V.append(make_vertex(Ndim, [i, i + 1, i + 2], False, zero_eig=True))
    elif kind == "mixed":  # everything, repeated sites
        for i in range(Ndim):
            V.append(make_vertex(Ndim, [i], True))
            V.append(make_vertex(Ndim, [i, (i + 2) % Ndim], False))
            V.append(make_vertex(Ndim, [(i + 1) % Ndim, i], True))
            V.append(make_vertex(Ndim, [i], True))
    return V


def main():
    Ndim = 10
    fail = False
    masks = {"all": [True], "half": [True, False], "rare": [False, False, True]}
    print(f"{'model':12s} {'dir':4s} {'mask':5s} {'nd':>5s} {'fx':3s} | ratio_dl  ratio_sm  G_dl      G_sm      G_std     | pad/B/W")
    for kind in ["hubbard", "jz", "sasb_zz", "tv", "tv_families", "kondo", "mixed"]:
        for direction in ["up", "down"]:
            for mname, mask in masks.items():
                for nd in [10**9, 3]:
                  for form in (["G", "M"] if "--both-forms" in sys.argv else ["M"]):
                    for xs in [False, True]:
                        V = model(Ndim, kind)
                        e, info = run(Ndim, V, direction, mask, nd_cap=nd, cross_slice=xs,
                                      sub_cls=(Submatrix if form == "G" else SubmatrixM))
                        bad = max(e.values()) > TOL
                        fail |= bad
                        print(f"{kind:12s} {direction:4s} {mname:5s} {('inf' if nd > 10**6 else nd):>5} {form}{'x' if xs else '-':2s} | "
                              f"{e['ratio_dl']:.1e}  {e['ratio_sm']:.1e}  {e['G_dl']:.1e}  {e['G_sm']:.1e}  {e['G_std']:.1e}  | "
                              f"{info['npad']}/{info['pathB']}/{info['pathW']}" + ("   FAIL" if bad else ""))
    print("RESULT:", "FAIL" if fail else "PASS", f"(tolerance {TOL})")
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
