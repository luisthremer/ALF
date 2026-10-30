"""Numerical check of rank-1 / delayed / submatrix update formulas in ALF conventions.

ALF: G = (1 + A)^{-1}, G_ij = <c_i c_j^dag>. A vertex n acts on site p with
e^{V_n} = 1 + (e^{g phi E} - 1) e_p e_p^T (diagonal, N=1, rank 1).
Flip phi -> phi': e^{V'} = (1 + Delta e_p e_p^T) e^{V}, Delta = e^{g(phi'-phi)E} - 1,
so A -> (1 + Delta P) A. Ratio (fermion part, one flavour) R = 1 + Delta (1 - G_pp).

Up-sweep (WRAPGRUP):   for n: wrap G <- e^{V_n} G e^{-V_n} (old field); propose/update at n.
Down-sweep (WRAPGRDO): for n (descending): propose/update at n; unwrap G <- e^{-V_n'} G e^{V_n'} (new field).
"""
import numpy as np

rng = np.random.default_rng(7)
N = 24
g = 0.4 + 0.3j          # complex coupling to exercise phases
phis = rng.choice([-2.0, -1.0, 1.0, 2.0], size=N)
sites = rng.permutation(N)          # vertex n acts on site sites[n]
accept = rng.random(N) < 0.6        # fixed accept mask
new_phi = -phis


def rand_A():
    M = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    return 0.3 * M


def run(direction, scheme, nd=1000):
    A = rand_A()
    G = np.linalg.inv(np.eye(N) + A)
    order = range(N) if direction == "up" else range(N - 1, -1, -1)
    # low-rank state
    G0 = G.copy()
    X = np.zeros((N, 0), complex); Y = np.zeros((N, 0), complex)   # delayed: G = G0 - X Y^T
    S = []; Gi = np.zeros((0, 0), complex)                          # submatrix: G = G0 - G0[:,S] Gam^{-1} (1-G0)[S,:]
    ratios = []
    phi = phis.copy()

    def cur_pp(p):
        if scheme == "standard":
            return G[p, p]
        if scheme == "delayed":
            return G0[p, p] - X[p, :] @ Y[p, :]
        # submatrix, p not in S
        assert p not in S
        return G0[p, p] + G0[p, S] @ Gi @ G0[S, p]

    def wrap(p, f, sign):  # G <- T G T^{-1}, T = diag with t at p
        nonlocal G, G0, X, Y, Gi
        t = np.exp(sign * g * f)
        if scheme == "standard":
            G[p, :] *= t; G[:, p] /= t
        elif scheme == "delayed":
            G0[p, :] *= t; G0[:, p] /= t
            X[p, :] *= t; Y[p, :] /= t
        else:
            G0[p, :] *= t; G0[:, p] /= t
            if p in S:   # frame change of Gamma^{-1}: Gam -> T_S Gam T_S^{-1}
                k = S.index(p)
                Gi[k, :] *= t; Gi[:, k] /= t

    def flush():
        nonlocal G0, X, Y, S, Gi
        if scheme == "delayed" and X.shape[1]:
            G0 = G0 - X @ Y.T
            X = np.zeros((N, 0), complex); Y = np.zeros((N, 0), complex)
        if scheme == "submatrix" and S:
            G0 = G0 - G0[:, S] @ Gi @ (np.eye(N)[S, :] - G0[S, :])
            S = []; Gi = np.zeros((0, 0), complex)

    for n in order:
        p = sites[n]
        if direction == "up":
            wrap(p, phi[n], +1)
            A = (np.eye(N) + (np.exp(g * phi[n]) - 1) * np.outer(np.eye(N)[p], np.eye(N)[p])) @ A @ \
                (np.eye(N) + (np.exp(-g * phi[n]) - 1) * np.outer(np.eye(N)[p], np.eye(N)[p]))
        Delta = np.exp(g * (new_phi[n] - phi[n])) - 1
        Gpp = cur_pp(p)
        R = 1 + Delta * (1 - Gpp)
        ratios.append(R)
        if accept[n]:
            ep = np.eye(N)[p]
            A = (np.eye(N) + Delta * np.outer(ep, ep)) @ A
            if scheme == "standard":
                G = G - np.outer(G[:, p], (ep - G[p, :])) * Delta / R
            elif scheme == "delayed":
                col = G0[:, p] - X @ Y[p, :]           # current G(:,p)
                row = G0[p, :] - X[p, :] @ Y.T          # current G(p,:)
                X = np.column_stack([X, col * Delta / R])
                Y = np.column_stack([Y, ep - row])
                if X.shape[1] >= nd:
                    flush()
            else:
                # border Gamma = D^{-1} + (1 - G0)_{SS}; new row/col for p
                w = -G0[S, p]                 # (1-G0)[S,p]
                q = -G0[p, S]                 # (1-G0)[p,S]
                sig = 1 / Delta + 1 - G0[p, p]
                # Schur complement: s = sig - q Gi w ; note R = Delta * s (det-lemma check)
                s = sig - q @ Gi @ w
                assert np.isclose(Delta * s, R), (Delta * s, R)
                Giw = Gi @ w; qGi = q @ Gi
                k = len(S)
                new = np.zeros((k + 1, k + 1), complex)
                new[:k, :k] = Gi + np.outer(Giw, qGi) / s
                new[:k, k] = -Giw / s
                new[k, :k] = -qGi / s
                new[k, k] = 1 / s
                Gi = new; S.append(p)
                if len(S) >= nd:
                    flush()
            phi[n] = new_phi[n]
        if direction == "down":
            wrap(p, phi[n], -1)
            A = (np.eye(N) + (np.exp(-g * phi[n]) - 1) * np.outer(np.eye(N)[p], np.eye(N)[p])) @ A @ \
                (np.eye(N) + (np.exp(g * phi[n]) - 1) * np.outer(np.eye(N)[p], np.eye(N)[p]))
    flush()
    Gfin = G if scheme == "standard" else G0
    Gex = np.linalg.inv(np.eye(N) + A)
    return np.array(ratios), Gfin, Gex


for direction in ("up", "down"):
    for nd in (1000, 5):
        rs, Gs, Gex = run(direction, "standard")
        rng = np.random.default_rng(7); rng.choice([-2.0, -1.0, 1.0, 2.0], size=N); rng.permutation(N); rng.random(N)
        for scheme in ("delayed", "submatrix"):
            rng = np.random.default_rng(7); rng.choice([-2.0, -1.0, 1.0, 2.0], size=N); rng.permutation(N); rng.random(N)
            r, Gf, _ = run(direction, scheme, nd)
            print(f"{direction:4s} nd={nd:4d} {scheme:9s}: max|dR/R|={np.max(np.abs(r - rs) / np.abs(rs)):.2e} "
                  f"max|dG|/max|G|={np.max(np.abs(Gf - Gs)) / np.max(np.abs(Gs)):.2e}")
        print(f"{direction:4s} standard vs exact: max|dG|/max|G|={np.max(np.abs(Gs - Gex)) / np.max(np.abs(Gex)):.2e}")
        rng = np.random.default_rng(7); rng.choice([-2.0, -1.0, 1.0, 2.0], size=N); rng.permutation(N); rng.random(N)
